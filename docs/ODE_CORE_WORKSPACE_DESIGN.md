# Step 3 design (DESIGN ONLY, v2): AD-safe private workspace for the remaining per-call RHS temporaries

Status: **v2; main approved increment 1 (private primal-only)**, implemented (results: `docs/ODE_CALLBACK_OPTIMIZATIONS.md`, Step 3 increment 1); increments 2–3 not started. v1 is kept outside git as
`A/chunk18/step3/ODE_CORE_WORKSPACE_DESIGN_v1_superseded.md`. Its §3 marker-tag claim was **wrong** and is corrected below. Implementation waits for
main's design gate. Everything proposed is **private and additive**: no exported symbol, public struct, StaticArrays, new dependency, global or
implicit/lazy loading. A public export would be a separate question for Marco.

## 1. What still allocates (current tree: F patch + Steps 1–2)

Per warm `recombination_rhs!` call (Float64): 12-state **1440 B / 20 objects**, 7-state **1040 B / 14 objects**. Full prediction: 17.4 M allocation events,
1.92 GB, GC median ≈ 156 ms of ≈ 2.05 s (Step 2 benchmark). Sites, from source (an exact `Profile.Allocs` inventory at sample rate 1 on single calls —
12/7-state, feedback, Float64 and inner-Jacobian Dual — is the **first task** of increment 1, before any code):

| site | object(s) |
|---|---|
| `ode_unpack` (`RecombinationODE.jl`) | `X` (15) |
| `recombination_rhs!` | `g` (15) |
| `fcn_effective!` | `dXH` (5); `dXHe` (4), allocated even when `flag_He = false` |
| `hydrogen_rhs!(dX, table, Tg, rho, …)` → `get_rates(table, Tg, rho*Tg)` (`RHS.jl:121`, `RateTable.jl:215–217`) | `A` (n), `B` (n), `R` (n × neq) in a `RateResult` |
| `helium_base_rhs!(dX, table, Tg, …)` → `get_helium_rates(table, Tg, c)` (`HeRHS.jl:130,201`) | `A`, `B` (n), `R` (n × n) |
| H-I absorber | to be measured (the `HIAbsorption.jl:106–107,172` arrays look like setup) |

The physics primitives already take explicit rate arrays: `hydrogen_rhs!(dX, Tg, Xe, Xp, NH, Hz, X, A, B, R, …)` and
`helium_base_rhs!(dX, Tg, Xe, NH, Hz, X, fHe, A, B, R, …)`. In-place `get_rates!` exists. A buffered path can therefore call the same primitives,
with **no duplicated physics**.

## 2. Element types come from the actual dependencies (corrected table)

The source dependencies are:
- `X = ode_unpack(y, fHe)` has type **promote(eltype(y), typeof(fHe)) only**.
- `bg = ode_background(rm, z; hscale, nbscale)`: `Tg = TCMB(z)` depends on z only; `NH` on z and nbscale; `Hz` on z and hscale.
- `g` has type promote(X, Tg, NH, Hz, z, `feedback_eltype(rm.diffusion)`) (`RecombinationODE.jl`, already in the source).
- H rates `get_rates(table, Tg, Te = rho·Tg)`: **A** depends on Tg and Te, so promote(Tg, rho); **B, R** depend on **Tg only**.
- He rates `get_helium_rates(table, Tg)`: A, B and R depend on Tg only.
- hscale/nbscale enter NH and Hz, hence g, but **not** the rates, unless the actual Tg or Te carries them (it does not: Tg = TCMB(z)).

| caller | y | z | X | g, dXH, dXHe | H A | H B, R | He A, B, R |
|---|---|---|---|---|---|---|---|
| primal (SciML `f`) | T | T | T | T | T | T | T |
| state Jacobian (`PreparedJac5`, Dual{J}) | D_J | T | D_J | D_J | D_J (via rho) | **T** | **T** |
| time gradient (Dual{t} in z) | T | D_t | **T** | D_t | D_t (via Tg) | D_t | D_t |
| hscale/nbscale Dual (D_h) | T or D_h (promoted state in a Dual solve) | T | as y | D_h | as y and Tg | T | T |
| outer ForwardDiff solve (4-parameter Jacobian) | every T above is Dual{Outer,Float64,k}; Jacobian/tgrad Duals nest inside | | | | | | |

The workspace therefore keeps **separate type parameters `TX, TG, TA, TBR, THe`** derived from these actual dependencies. It preserves the
`RateResult{TA,TB,TR}` split **from the start**. Promoting B/R to A's type is **not** assumed bitwise-safe. Each workspace instance is
built per solve from the actual (y, z, hscale, nbscale, feedback) types, and every call checks compatibility. A mismatch goes to the unchanged
allocating public path (counted), never a conversion. **Baseline capture before implementation:** the currently valid mixed and nested paths
(the rows above, plus BigFloat state) are recorded as bitwise reference outputs from the current allocating path.

## 3. ForwardDiff tag construction (corrected; proven on a toy, not yet on CosmoRec)

**Source fact** (installed ForwardDiff 1.4.6, `src/config.jl:34–37`): `checktag(::Type{Tag{F,V}}, f::F, x::AbstractArray{V}) = true` and any other
`Tag{FT,VT}` throws `InvalidTagException`. The default check (in `jacobian!`, via the config) therefore requires the callable `f isa F` of the tag.

- **v1's empty-marker idea is wrong:** `Tag(Marker(), Float64)` with a `Buffered{W}` functor throws `InvalidTagException` (reproduced, item 1 below).
- **The default concrete tag `Tag{Buffered{W},V}` is recursive** when `W` holds `Dual{Tag{Buffered{W},V}}` buffers.
- **Family tag (main's construction):** `Tag{Buffered,V}()`, with `Buffered` the UnionAll functor family, satisfies `f isa Buffered` under the **default**
  check. It needs no `Val(false)`, no custom global method and no Tuple bypass, and it breaks the recursion because the tag no longer contains `W`.
- **Ordering:** `ForwardDiff.Tag(f, V)` triggers the internal `@generated tagcount` at construction. A family tag must do the same, by calling
  `ForwardDiff.tagcount(Tag{Buffered,V})` once when the cache is built. That is an internal (unexported) ForwardDiff function, but it is the exact
  mechanism `Tag(f, V)` uses. The tag is built per solve, inside any outer ForwardDiff context, so the inner count exceeds the outer one.

**Toy proof** (`A/chunk18/step3_tag_probe.jl` → `step3/tag_probe.txt`, exit 0; ForwardDiff only, no CosmoRec). A nonlinear 4-component RHS, a functor
writing through a Dual workspace, the default check, vs a closure-based oracle (default tags):

| check | result |
|---|---|
| 1) marker tag with `Buffered` functor | `InvalidTagException` |
| 2) family tag, primal J | bitwise |
| 3) one cache with alternating live s/z; two caches interleaved | both equal to the oracle |
| 4a) d J/d s, outer Dual over a parameter (nested) | bitwise |
| 4b) d J/d z | bitwise |
| 4c) d J/d u (Jacobian of the Jacobian) | bitwise |
| 5) d²J/dz² (two nested outer levels) | bitwise |
| 6a) outer variable entering both the state and the parameter | bitwise |
| 6b) classic perturbation-confusion form x · Σ J(u; s = x) | bitwise |
| 7) ordering | outer ≺ inner = true, inner ≺ outer = false |

- In 4–6 the **inner** level is family-tagged and the outer levels use default closure tags.
- **Not tested:** two family-tagged caches with the *same V* nested inside one computation (both levels family-tagged at equal V). A shared tag type would then risk perturbation confusion. Our call graph never nests a prepared Jacobian inside another at the same element type: the inner V always contains the outer Dual. The concrete default tags already used by `jac5!`/`PreparedJac5` have the same per-type sharing.
- Probe attempts 1–2 failed on harness bugs (a constructor that overwrote the default constructor; an unpromoted oracle state); they are kept as `step3/tag_probe_attempt{1,2}*`.
- **Before the core plan is approved**, the same checks must pass on CosmoRec's real RHS in increment 2: the outer-Dual (4-parameter) solve, d/dz of the state Jacobian, two caches and composition.

## 4. Scope and API (private entry point; public allocating path stays the oracle)

- **New private entry point**, not a new method of an exported generic: e.g. `CosmoRec._recombination_rhs_ws!(f, z, y, rm, ws; flag_He, hscale, nbscale)`
  and `CosmoRec._RHSWorkspace` (both unexported, underscore-named). Exported `recombination_rhs!`/`recombination_rhs`/`fcn_effective!` are **not
  changed**: they keep their current allocating implementation and serve as the oracle. They do **not** build a full workspace per call, which could
  regress AD, types or allocations, especially with helium off.
- **Shared math, no duplication:** the private path calls the existing primitives with workspace buffers. The only new helpers are private in-place
  wrappers: `_ode_unpack!`, a private `_get_helium_rates!` mirroring `get_helium_rates`, and buffer-taking private variants of the
  `fcn_effective!` assembly (`dXH`, `dXHe`, rates from the workspace; `dXHe` untouched when `flag_He = false`).
- **Only the opt-in reference callback** (test helpers) uses the workspace initially. No public struct, StaticArrays or dependency decision is
  needed for private additive internals. Exporting any of it would be a question for Marco.

## 5. Increments (each behind its own gate)

1. **Primal-only buffered path:** private workspace (`TX, TG, TA, TBR, THe`), private entry point, primal callback functor in the test helpers (per
   solve, counted fallback). The Jacobian (`PreparedJac5`) and tgrad stay unchanged (they keep using the allocating RHS).
2. **Jacobian workspace** inside the Jacobian cache, with the **family tag** of §3. This needs a separate gate with the CosmoRec-level tag checks.
3. **tgrad workspace**, later and separately. No all-storage rewrite.

## 6. Gates per increment

- **Numerics:** private path bitwise vs the public allocating oracle at all 5a rows (12/7-state, off-table DP fallback z = 1500, inactive absorber, feedback model); public X_e/T_m byte-identical (fiducial + 4 cosmologies); identical per-solve stats; 10a's known failures unchanged.
- **AD:** the captured mixed/nested baselines (§2) bitwise; ForwardDiff through the callback paths; 10b 4-parameter Jacobian bitwise; 5e unchanged; prepared DI-Mooncake (5a VJP, 5f) within existing tolerances.
- **Reverse-mode aliasing:** test workspace reuse **within one Mooncake-taped function** (two RHS evaluations sharing buffers), because reverse mode must not read buffers overwritten later.
- **Allocation:** `@allocated` per call, primal path → near 0 B for increment 1; full-call allocations and bytes.
- **Timing:** matched, thermally logged full-call benchmark (warm, in-process cool, pinned CPU 8, unchanged rule), incremental vs the previous step.
- **Full suite:** once at the end of the accepted batch; the known 6 failures with identical values.

## 7. Risks

- **Internal ForwardDiff API:** `tagcount` is unexported; a ForwardDiff upgrade could change it. Pin it with a regression test: the tag-ordering assertion plus nested checks.
- **Reverse-mode aliasing** (§6).
- **Compile-time growth.**
- **Scope:** the H-I absorber internals and PDE stages are separate targets.

## 8. CosmoRec-level proof for increment 2 (external, no production change) — proof gate FAILED on one check; counterfactual: pre-existing limitation

**Proof** (`A/chunk19/inc2_jacws_proof.jl`). The design under test is `BufRHS` (a functor writing through a Dual-typed private `_rhs_workspace`) and
`JacWS` (`JacobianConfig(…, Tag{BufRHS,V}())`, default tag check, `tagcount` triggered explicitly). Attempt history:
- **attempt 1** (`inc2_proof_20261003T174133`): harness bug, a Float64 `J` received a Dual Jacobian;
- **attempt 2** (`inc2_proof_retry_20261003T185714`): harness `ParseError` (misplaced `end`). My "PARSE OK" pre-check was invalid: `Meta.parseall` returns `:error` nodes and does not throw. It was replaced by a tree walk, verified to flag attempt 2;
- **attempt 3** (`inc2_proof_retry_20261003T191335`): exit **1**, `PROOF RESULT: 1 FAILURE(S)`. Every mismatch, error or fallback counts.

| check (vs the unprepared oracle `jac5!`) | 12-state z = 2500 | 12-state z = 1500 (off-table DP fallback) | 7-state z = 400 |
|---|---|---|---|
| primal J | bitwise | bitwise | bitwise |
| dJ/du, dJ/dz, dJ/dhscale (outer Dual over the family-tagged inner Jacobian) | bitwise | bitwise | bitwise |
| d²J/dz² (two nested outer levels) | bitwise | bitwise | bitwise |
| two caches, alternating live (y, p, z) | = oracle, 0 fallbacks | = oracle, 0 fallbacks | = oracle, 0 fallbacks |
| prebuilt cache, prepared Mooncake, one tape, two Jacobians, vs ForwardDiff | 4.132e-14 | **4.450e-01 (FAIL)** | 5.191e-16 |

- **Also passed:**
  - tag ordering: outer ≺ inner true, inner ≺ outer false;
  - whole public call with the Jacobian-workspace callback: 0 fallbacks; X_e/T_m byte-identical to the primal baseline; all per-solve stats identical;
  - 10b 18×4 ForwardDiff Jacobian (Dual solves) byte-identical, 0 fallbacks;
  - full-call allocations 983759392 B / 4491547 → 501099568 B / 2801905 (one sample);
  - Jacobian kernel 0 B (vs 5776/3888 B), indicative.

**Bounded counterfactual** (`A/chunk19/inc2_counterfactual_z1500.jl` → `inc2_counterfactual_20261003T192241/counterfactual.txt`, exit 0). Same state, W1/W2 and
perturbation; prepared DI-Mooncake vs the ForwardDiff gradient through the oracle:

| 12-state z = 1500 | A oracle `jac5!` | B `PreparedJac5` | C family tag + allocating RHS | D1 fresh Dual workspace per call | D2 one reused Dual workspace |
|---|---|---|---|---|---|
| single-J ⟨W1, J(u)⟩ | 5.117e-01 | 5.117e-01 | 5.117e-01 | 5.117e-01 | 5.117e-01 |
| two-J ⟨W1, J(u)⟩ + ⟨W2, J(u·cs)⟩ | 4.450e-01 | 4.450e-01 | 4.450e-01 | 4.450e-01 | 4.450e-01 |

At the in-table control z = 2500, all variants give 2.832e-14 (single-J) and 4.132e-14 (two-J). In every row the repeated gradient is identical, the primal objective
equals the reference, and there are 0 fallbacks.

**Conclusion.** The disagreement is **pre-existing**. The unprepared oracle, the current `PreparedJac5`, the family tag with the allocating RHS, and the
buffered design (fresh or reused) all give the same Mooncake gradient of J at the off-table DP-fallback state. It is therefore **not** caused by the family tag,
the Dual workspace or workspace reuse/aliasing; it is not a second-call effect either, since the single-J objective fails as well. It is a disagreement between
Mooncake and ForwardDiff for **derivatives of the state Jacobian** (second order) through the off-table DP-fallback path. Which one is correct is **not
established** here: an independent arbiter, e.g. a high-precision central difference of ⟨W, J(u)⟩, would decide. Earlier Mooncake-through-Jacobian coverage
(Step-2 probe, test 5f) used only z = 2500 and z = 400, so this case was untested before.

Status: **increment 2 NOT approved**; no production Jacobian change. Main decides the next gate (e.g. arbitration of the z = 1500 second-order derivative,
and whether the gate should require Mooncake-through-J agreement on the off-table branch, given that the baseline shows the same behaviour).

### 8.1 Follow-ups (main): componentwise parity, then bounded arbitration (external; no production, tolerance or clamp change)

**Componentwise parity** (`A/chunk19/inc2_parity_z1500.jl` → `inc2_parity_20261003T195301/parity.txt`, vectors in `parity.txt.vectors.txt`; exit **0**).
- For both objectives (single-J, two-J) at 12-state z = 1500 and the z = 2500 control, the full 12-component Mooncake gradients of B (`PreparedJac5`), C (family tag + allocating RHS), D1 and D2 (buffered, fresh / reused Dual workspace) are **bitwise identical** to the legacy oracle A.
- Primal objectives equal the reference, repeated gradients are identical, and there are 0 fallbacks.
- This is a **no-new-regression** result; it is not evidence that Mooncake is correct.

**Arbitration** (`A/chunk19/inc2_arbitration_z1500.jl` → `inc2_arbitration_20261003T195402/arbitration.txt`, exit 0, 5 min). Single case: 12-state z = 1500, s(u) = ⟨W1, J(u)⟩, J = `jac5!`.

| step | result |
|---|---|
| fallback calls captured in one primal RHS | singlet orders (4,4,5,4,4,4,5) at Tg 4090.98, η 6.35e7, τS 2.44e7, pd 1.69e-3; triplet orders (4,4,4,2,2,2,2), τS 2.49, pd 1 |
| frozen-order diagnostic model | primal RHS and primal J **bitwise** equal to the adaptive model |
| ForwardDiff vs Mooncake, adaptive orders | 5.117e-01 (max-norm relative to max\|g_F\|); worst component k = 10 |
| the same with **frozen** orders | 5.117e-01; ForwardDiff adaptive == frozen, Mooncake adaptive == frozen, so the adaptive Patterson branch is **not** the cause |
| `hi_absorption = false` diagnostic model | **2.090e-16**: the disagreement is confined to the H-I absorber path |
| component k = 10 (He state, y₁₀ ≈ 6.4e-17), frozen | ForwardDiff **1.4193e15**, Mooncake **−5.0463e15** |
| BigFloat(256) ForwardDiff-over-ForwardDiff, frozen | **4.2321119038513456e14** (97 s) |
| BigFloat(256) central differences of s along e₁₀, h = 6.4e-19 … 6.4e-25, and Richardson | **4.2321119038513456e14** at every step |
| Float64 central differences, same steps | exactly 0: the perturbation is absorbed by Float64 rounding, so they are unusable as an arbiter |

- **Reading (factual):** BigFloat ForwardDiff and BigFloat fixed-step FD agree at the **17 digits displayed** (values were printed after a Float64 cast). Agreement
  at full 256-bit precision was **not** shown in this run; see §8.2.
- **Neither Float64 method reproduces that value.**
  - The global difference normalized by max|g_F(Float64)| is 7.883e-2 for ForwardDiff and 4.329e-1 for Mooncake.
  - For the **individual component k = 10**, Float64 ForwardDiff 1.4193e15 vs BigFloat 4.2321e14 is a **235 %** relative error, and Float64 Mooncake −5.0463e15 has the **wrong sign**.
- **Demonstrated:** a precision-sensitive Float64 second-derivative path through the absorber at this state, for this component.
- **Not proven:** the mathematical conditioning or the culprit operation; other components or states (only k = 10 at z = 1500 was arbitrated, within the cap); any mechanism inside Mooncake. "Truth" here means the 256-bit evaluation of the same code.
- This is consistent with the parity result: the behaviour is a property of the existing absorber code under Float64 second-order AD, unchanged by the cache, tag or workspace.

**Implication for the gate (main's decision).**
- An off-table Float64 correctness gate cannot use a reference that is itself ~8 % (globally normalized) inaccurate.
- The **regression gate** at this state is therefore:
  - bitwise legacy-vs-new Mooncake vector parity;
  - forward and nested ForwardDiff value parity;
  - high-precision new-vs-old parity (§8.2).
- The original Float64 Mooncake-vs-ForwardDiff **FAIL is retained** in the evidence, neither labelled correct nor hidden.
- The in-table and control AD accuracy gates are **unchanged**.
- This is a scope-specific non-regression gate. It is **not** a claim of accurate Float64 off-table Hessians, nor of a correct true reverse through Rodas.

### 8.2 Final bounded gate: high-precision new-vs-legacy parity at the frozen z = 1500 point — PASSED

`A/chunk19/inc2_bigparity_z1500.jl` → `inc2_bigparity_20261003T202459/bigparity.txt` (full-precision strings), exit **0**.
- Attempt 1 (`inc2_bigparity_20261003T202237`) failed on a harness bug, kept on record: a 3-argument outer constructor overwrote the 3-field default constructor.
- Setup: frozen orders captured from the primal fallback calls (singlet (4,4,5,4,4,4,5), triplet (4,4,4,2,2,2,2)); the frozen model's primal RHS is bitwise equal to the adaptive one.
- Measured quantity: ForwardDiff-over-ForwardDiff gradient of s(u) = ⟨W1, J(u)⟩ in BigFloat, with J from the legacy oracle `jac5!` vs the NEW `JacWS` (Dual private workspace, family tag, default check).

| precision | new vs legacy, all 12 components | fallbacks | time legacy / new |
|---|---|---|---|
| 256 bit | **bitwise identical** | 0 | 91.2 s / 111.9 s |
| 128 bit | **bitwise identical** | 0 | 0.2 s / 0.2 s |

Component k = 10 at full precision:
- 256 bit: 4.232111903851345627861555251935371480329198577529797195630165000496252092867308e+14 (legacy = new);
- 128 bit: 4.232111903851345627861471556457063414186e+14 (legacy = new).

Central-FD ladder of s along e₁₀ (relative error vs the ForwardDiff value at the same precision; legacy and new identical):

| h | 256 bit | 128 bit |
|---|---|---|
| 6.435e-19 | 1.236e-35 | 5.872e-18 |
| 6.435e-21 | 1.236e-39 | 1.123e-16 |
| 6.435e-23 | 1.236e-43 | 4.094e-14 |
| 6.435e-25 | 1.237e-47 | 3.830e-13 |

- At 256 bits the FD error falls by 1e-4 per factor-100 step, i.e. O(h²) truncation, so the 256-bit ForwardDiff value is confirmed to ~47 digits.
- At 128 bits the ladder is roundoff-limited for small h.
- The 128- and 256-bit ForwardDiff values of this component agree to ~22 significant digits, out of ~38 for 128 bits. This is consistent with a strongly precision-amplifying evaluation path for this second derivative, and hence with the Float64 results (§8.1) losing the component entirely. Consistent with, not a proven mechanism.

**Gate status:** the regression gate of §8.1 holds:
- bitwise legacy-vs-new Mooncake vectors;
- bitwise forward and nested ForwardDiff values;
- bitwise high-precision new-vs-legacy parity at 128 and 256 bits.

The Float64 off-table Mooncake-vs-ForwardDiff FAIL remains recorded as a pre-existing, precision-sensitive behaviour of the existing absorber path. It is **not** fixed and **not** claimed accurate. Increment-2 implementation awaits main's review of this gate and of the plan in §9.

## 9. Increment-2 implementation plan (main approved the proof gate; IMPLEMENTED as planned — results in `docs/ODE_CALLBACK_OPTIMIZATIONS.md`, Step 3 increment 2)

**Minimal callback patch** (test helpers only, `test/chunk5_helpers.jl`; no `src` change is needed, since the private `_rhs_workspace` /
`_recombination_rhs_ws!` from increment 1 already accept Dual element types):
- `mutable struct BufRHS5{W,P,Z}`: Dual private workspace plus LIVE `p`/`z` and a fallback counter; calls `_recombination_rhs_ws!` and counts a `false` return.
- Family tag `Tag{BufRHS5,V}()`, with `ForwardDiff.tagcount` triggered once at construction (explicit internal dependency, documented).
- `PreparedJac5` gets the buffered path: built per solve from the actual `u0`/`p`/`z0` types. On a type mismatch, the existing unprepared oracle `jac5!` runs and is counted, never converted. The default tag check stays on; no `Val(false)`, piracy or Tuple bypass; no globals.
- `tgrad5!` unchanged (later, separate gate). `ODEFUN5` and the allocating public RHS remain the oracles.

**Normal unconditional regression tests** (new file `test/chunk5h_jacobian_workspace.jl`, added to `runtests.jl`):
1. **Primal J** bitwise vs `jac5!` at the 5a rows (12-state in-table z = 2500, off-table z = 1500, 7-state z = 400).
2. **ForwardDiff dJ/du, dJ/dz, dJ/dhscale, d²J/dz²** bitwise vs the oracle (cache built from the Dual types).
3. **Two caches** alternating live (y, p, z); no stale parameters; 0 fallbacks.
4. **Prepared DI-Mooncake, one tape, two Jacobian evaluations** through one prebuilt cache:
   - in-table states (z = 2500, z = 400): agreement with ForwardDiff within the existing 1e-12 accuracy gate (unchanged);
   - off-table z = 1500: **bitwise parity with the legacy oracle's Mooncake vector** (the no-new-regression gate of §8.1). It is labelled as such and is not an accuracy claim. The known Float64 Mooncake-vs-ForwardDiff difference there is recorded in this doc, not asserted as correct.
5. **Tag ordering** assertion (pins the internal `tagcount` behaviour).
6. **5b:** the solver's Jacobian is the buffered cache; a real block is bitwise vs the allocating callback; 0 fallbacks.
7. **Allocation:** buffered Jacobian per call below the current `PreparedJac5`.

The high-precision new-vs-old check (§8.2) stays an external gate artifact: about 2 × 97 s at 256 bits is too slow for the normal suite.

**Validation sequence after the patch:** the focused native + AD suites, the 5 histories (bitwise), the 18×4 10b Jacobian (bitwise), then the
matched benchmark (warm then in-process cool, pinned CPU 8, unchanged rule; incremental vs increment 1). The root `Pkg.test` runs once at the end of the accepted batch.
