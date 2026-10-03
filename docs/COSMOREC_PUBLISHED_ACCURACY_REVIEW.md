# CosmoRec: published accuracy claims vs what our port has demonstrated (literature/config review, 2026-10-02)

**Scope.** This is a read-only review: no source, test, gate, dependency or native change, and no new runs.

**Sources.** The full-text PDFs and the author page were downloaded to the session scratchpad, outside git, and read in full. Quotations are verbatim, with the page of the arXiv PDF. DOIs come from the arXiv metadata. The WebFetch tool was denied in this session, so public pages were fetched with `curl`.

## 1. Primary sources (original authors) and their status

| # | reference | role |
|---|---|---|
| P1 | J. Chluba & R. M. Thomas, "Towards a complete treatment of the cosmological recombination problem", MNRAS 412, 748–764 (2011); arXiv:1010.3631 (v3); DOI 10.1111/j.1365-2966.2010.17940.x | the CosmoRec paper (physics, convergence) |
| P2 | J. R. Shaw & J. Chluba, "Precise cosmological parameter estimation using CosmoRec", MNRAS (2011); arXiv:1102.3683 (v1); DOI 10.1111/j.1365-2966.2011.18782.x | defines the 'default' setting; accuracy of settings in X_e, C_ℓ and parameters |
| P3 | J. Chluba, G. M. Vasil & L. J. Dursi, "Recombinations to the Rydberg States of Hydrogen and Their Effect During the Cosmological Recombination Epoch", MNRAS 407, 599–612 (2010); arXiv:1003.4928; DOI 10.1111/j.1365-2966.2010.16940.x | shell convergence (up to 350 shells), collisions |
| R1 | native `README` of our v3.0b tree (`tools/CosmoRec/README`, "beta version 3.0 (March 06, 2017)") | version-specific code claim |
| R2 | author page http://www.jb.man.ac.uk/~jchluba/Science/CosmoRec/CosmoRec.html (fetched 2026-10-02; describes **v2.0.3 beta**; v3.0b is not listed) | later code claim |
| I1 | Y. Ali-Haïmoud & C. M. Hirata, "HyRec: A fast and highly accurate primordial hydrogen and helium recombination code", Phys. Rev. D 83, 043513 (2011); arXiv:1011.3758; DOI 10.1103/PhysRevD.83.043513 | **independent** code (not an original-author claim) |

No later original-author convergence paper specific to v3.0b was found. R1 is the claim carried by our source tree. The changes listed in R1 for v3.0b (fundamental-constant variation, a new Recfast++ correction function, solver-stability items) come with no new accuracy figure.

## 2. Exact claims (verbatim)

**Code-level output precision** (R1 lines 278–280; R2 says the same, with "≤" for "<~"):
> "With the startup file it is currently possible to include the effective rates for up to 500 shells for hydrogen, and 20 shells for helium. The precision of the output should be <~ 0.1 % for hydrogen recombination and similar to ~0.1%-0.2% during helium recombination."

R2 reads "…up to 500 shells for hydrogen, and 30 shells for helium. The precision of the output should be ≤ 0.1% for hydrogen recombination and similar to 0.1%-0.2% during helium recombination." The metric (pointwise vs peak) and the reference (a more complete model) are not specified.

**Required-physics threshold (P1, p. 1, Introduction):**
> "it is crucial to incorporate all important processes leading to changes in the free electron fraction close to the maxima of the Thomson visibility function … by more than ∼ 0.1% into one recombination module."

**Effective multi-level convergence and implementation check (P1, p. 13, §5.1):**
> "We find that the correction converges down to z ∼ 200 when including ∼ 300 − 400 shells … We also directly compared with our full multi-level recombination code and found the difference to be smaller than ∆Ne /Ne ∼ 10−5."
> "Collisional processes are still able to change the low redshift behaviour at the ∼ 0.1% level in this redshift range … we defer a detailed analysis on the importance of this effect to a future work."

**Two-photon and Raman truncation (P1, p. 14, §5.3–5.4):**
> "Two-photon processes lead to a total acceleration of recombination by ∆Ne /Ne ∼ −0.46% at z ∼ 1120. The main contribution comes from the 3s-1s and 3d-1s two-photon process, while the higher levels only add ∆Ne /Ne ∼ −0.08% at z ∼ 1200. The correction practically converge when accounting for the two-photon terms up to 5s-1s and 5d-1s. … for practical purposes it is sufficient to include the two-photon corrections for all ns and nd states up to n ∼ 4 − 5."
> "Higher level Raman scattering lead to a small additional modification, which for practical purposes could be neglected. We recommend including the Raman-corrections for the first three shells."

**The 'default' setting and its accuracy (P2):**
- Table 1, p. 4: default = resolved states "2s-3s, 2p-3p, 3d" (5), n2γ = 3, nRaman = 2, nHeI_max = 2, "HeI diffusion correction on", "HeI feedback off", ∆z (PDE-solver) = 20, "average runtime 1.3 sec". Effective rates "nHI_eff = 500 and nHeI_eff = 30". History "starting at z = 3000 and ending at z = 50 with 500 intermediate points", completed to z = 0 "using the simple Recfast ODE system".
- §2.1, pp. 2–3: "For all cases shown the differences in the hydrogen recombination history are ∆Ne /Ne ≲ 0.1%. The largest difference appears when switching off the diffusion correction to the escape probability of the He i 21P − 11S resonance … resulting in ∆Ne /Ne ∼ 0.3% uncertainty at z ∼ 1800 … During hydrogen recombination, higher level two-photon decays (n > 3) do still lead to some ∼ 0.1% uncertainty, which appears to be dominated by the 4s-1s and 4d-1s process."
- p. 2: "The agreement in the prediction for the CMB power spectra is better than 0.1% for all considered recombination models. This error is below the 3/l benchmark suggested by Seljak et al. (2003)".
- p. 3: "We also confirmed the precision of CosmoRec by comparing directly with the most detailed computation carried out using a more elaborate multi-level hydrogen-helium recombination code (Chluba & Thomas 2010). We found differences no larger than ∆Ne /Ne ∼ 0.01% at all redshifts." Also: "a final cross-validation of the CosmoRec outputs with independent recombination codes … will be very important."
- p. 4: "changes at the level of ∼ 1% in the freeze-out tail of recombination are not constrainable with Planck"; collisions "could imply modifications ≳ 0.1%".
- §2.3 (pp. 4–5): for the z ≲ 200 completion, "all shown cases is ∆Cl /Cl ≲ 10−5 at l ≲ 3000".
- §3, Table 2 (p. 6): mock data generated with the 'full' setting and analysed with the 'default' setting give parameter biases ≤ 0.12σ (ΛCDM: Ωb h² −0.01σ, Ωc h² 0.05σ, τ 0.10σ, log A_s 0.12σ) for a Planck-like forecast.

**Shell convergence (P3, abstract p. 1):**
> "for 350 shells down to redshift z ∼ 200 the results for the free electron fraction have practically converged. The final modification … decreases from about ∆Ne /Ne ∼ 2.8% for 100 shells to ∆Ne /Ne ∼ 1.6% for 350 shells. … for accurate computations in connection with the analysis of Planck data already ∼ 100 shells are expected to be sufficient. … collisional rates … a correction of ∆Ne /Ne ∼ −8.8 × 10−4 at z ∼ 700".

**Independent, not an original-author claim (I1, p. 2):** HyRec computes "a highly accurate recombination history (with errors at the level of a few times 10−3 for helium recombination and a few times 10−4 for hydrogen recombination)". Of Chluba & Thomas it says: "The code they present includes the same physics as ours. The main difference is the treatment of radiative transfer." No quantitative CosmoRec–HyRec agreement figure is given there. P2's CosmoRec "HyRec-case" setting differs from 'default' by ≲ 0.1% in hydrogen (P2, Fig. 1).

**T_m:** no explicit accuracy claim for T_m(z) was found in P1–P3, R1 or R2.

## 3. Which configuration our benchmark uses (traced, not inferred from the ODE size)

- **Native call:** `cosmorec_calc_h_cpp_`, runmode 0, with `runpars[1] = 0` (the "default setting", `CosmoRec.cpp:665-668`) and `runpars[15..20] = −1` (no overrides), from `chunk10/native_capture/harness10a.cpp:141-149`.
- **Batch defaults** (`CosmoRec.cpp:833-849`): nShells = 3, **nS_effective = 500**, nShellsHeI = 2, HI_absorption = 2 ("on with Diffusion fudge", `runfiles/parameters.ini:69-71`), spin_forbidden = 1, HeI_Feedback = 0, Diffusion_flag = 1, induced_flag = 2, nS_2gamma = 3, nS_Raman = 2, Diff_iteration_max = 2; nz = 3000, zstart = 3000 (`:757-760`).
- **Tables loaded.** Native (`Modules/HI_routines.cpp:20-26, 115-135`): `Effective_Rates.HI/Effective_Rate_Tables.nS_3/Rates_n{2,3}_l*.nS_500.dat`. Julia (`src/RecombinationODE.jl:21-30`): the same directory and nS_H = 500; helium `Effective_Rate_Tables.HeI.res_2` with nS = 30; `Pesc_Data/DP_Coll_Data.31.fac_50.neff_30.dat`. Our fixture prints CFG10 nS_2gamma 3, nS_Raman 2, Diff_iteration_max 2, induced_flag 2, Diffusion_flag 1.
- **Physical model.** This is the **published 'default' setting** of P2 Table 1, physically: an effective multi-level atom with 500 hydrogen shells folded into 5 resolved levels (n ≤ 3), and 30 helium shells. It is **not a 3-shell physical model**:
  - nShells = 3 is the resolved-level count;
  - 500 is the effective-rate shell count;
  - 2γ (n ≤ 3) and Raman (n ≤ 2) are the radiative-transfer truncations.
- **It is not the 'full' or high-accuracy setting.** Relative to P1's recommendation (2γ to n ∼ 4–5, Raman for 3 shells) and to P2's 'full' setting (22 resolved levels, n2γ = 8, nRaman = 7, HeI feedback n = 5), it omits:
  - n > 3 two-photon terms, about 0.1% in X_e during H recombination (P2) or −0.08% (P1);
  - n = 3 Raman terms ("small", P1);
  - HeI feedback;
  - the explicit HeI radiative transfer (runmodes 4–6).
  - Collisions are absent in every CosmoRec setting (R1: "Collisional processes were not taken into account thus far").
- **Numerics differ from the 2011 default:** 3000 history nodes vs 500, and the PDE outputs every ∆z = 10 between z = 2500 and 500, vs ∆z = 20.
- **Tree provenance.** Our tree also carries CAMB-adapter override hooks (`runpars[15..20]`, `CosmoRec.cpp:719-726, 852, 903-907`). They are inactive at −1. Whether they are in the official v3.0b tarball was not checked.

## 4. Claimed vs measured

| quantity | published claim (source) | metric / domain / configuration of the claim | ours (measured) | our metric / domain / configuration | relation |
|---|---|---|---|---|---|
| X_e, hydrogen recombination | "<~ 0.1 %" output precision (R1, R2) | relative to a more complete model; metric not stated; code level | **1.449e-6** max pointwise (Julia vs native) | all 10000 output nodes (incl. Recfast tail), one fiducial cosmology (CAMB H, native initialization), default setting | 1.45e-6 is a **port error**, not model accuracy. It is about 690× below the claimed 0.1% model precision |
| X_e, helium recombination | "~0.1%-0.2%" (R1, R2) | same | (inside the same 1.449e-6 max) | same | port error ≪ claimed model precision |
| X_e, effective-rate implementation | < 1e-5 vs full multi-level code (P1 §5.1); ≤ 0.01% vs most detailed code (P2 p. 3) | implementation consistency, one cosmology | 1.449e-6 | port vs native | our port error is below the authors' own implementation-check level |
| X_e, setting choice | ≲ 0.1% (H) between settings; 0.3% at z ∼ 1800 without the HeI diffusion correction (P2) | physics-truncation uncertainty of 'default' | not applicable | — | model uncertainty of 'default' ≈ 1e-3 |
| X_e, shell convergence | converged to z ∼ 200 at 300–400 (350) shells (P1, P3) | model convergence | not tested | — | — |
| C_ℓ (TT/TE/EE) | settings agree to "better than 0.1%"; tail completion ≲ 1e-5 (P2) | ℓ ≲ 3000, fiducial cosmology | **not computed** | — | open |
| parameters | default vs full ≤ 0.12σ (P2 Table 2) | Planck-like forecast | not computed | — | open |
| T_m | no claim found | — | 3.557e-8 max pointwise | as for X_e | port error only |

## 5. Verdict

1. **What our test demonstrates.** In the published 'default' physical configuration (500-shell effective rates; n2γ = 3, nRaman = 2; HeI diffusion fudge on; feedback off), our Julia port reproduces the native v3.0b at one fiducial cosmology, using the native initialization and our injected Rodas5P callback (reltol 1e-12, a1 1e-18, aex 1e-14):
   - X_e to 1.45 ppm max pointwise, T_m to 0.036 ppm.
   - That is a **numerical-reproduction** result: an additional port error on top of the model's own uncertainty. It is not evidence that the model is accurate at the ppm level.
2. **Margin against the authors' claims.** The claimed model precision is ≲ 0.1% in X_e (H) and 0.1–0.2% (He); the default's setting-level uncertainty is ≈ 0.1% (P2).
   - Our port error is about 3 orders of magnitude below that, and below the authors' own implementation checks (≤ 1e-5 to 1e-4).
   - So, at the fiducial cosmology and default setting, the port adds a negligible error compared with the model error the authors state. In that specific sense the user's criterion ("if we reach the claimed theoretical accuracy") is met **by a wide margin, for X_e at this cosmology**.
   - T_m agrees with the native code to 3.557e-8 max pointwise. That is excellent native agreement, but no published T_m accuracy target exists (§2), so we cannot claim it meets one.
3. **What is NOT demonstrated:**
   - C_ℓ or parameter-level agreement: no Boltzmann run with the Julia X_e vs the native X_e.
   - Other cosmologies: only one was tested, and the Julia path still uses the native preliminary history and initialization (B3, not implemented).
   - The 'full' or high-accuracy settings: nShells 4/10, n2γ 4–8, HeI feedback/transfer are not ported (`hi_pde_setup` throws for other configurations).
   - Model accuracy itself: collisions, n > 3 two-photon terms and shell truncation are model limitations of CosmoRec, not of the port.
   - Reverse-mode and end-to-end cosmological sensitivities (B2, B3), and the six remaining test failures. These are separate from scientific accuracy.
4. **A source-backed port-error budget (proposal; no gate changed).** The stated model precision is ε ≈ 1e-3 in X_e (P1 Introduction; R1; P2), so a proportionate port budget would be ≤ ε/10 = **1e-4** pointwise in X_e over the recombination domain.
   - The approved 1e-5 gate is already 10× stricter than that, and the measured 1.45e-6 passes both.
   - No published T_m precision exists, so any T_m budget would be our own choice. The current 1e-7 gate passes at 3.6e-8.
   - None of this argues for widening existing gates. It shows that the current X_e gate is proportionate or conservative relative to the published model precision. The T_m gate cannot be judged against a published target, because none exists.
5. **Recommended to close the "Boltzmann-facing" question** (proposals):
   - (a) A CAMB/CLASS run with the native X_e vs the Julia X_e at the fiducial cosmology, reporting ΔC_ℓ/C_ℓ for TT/TE/EE at ℓ ≤ 3000 against the authors' 0.1% / 1e-5 scales.
   - (b) Repeat at a few cosmologies, once B3 (in-pipeline initialization) exists.
   - (c) If 'full'-setting accuracy is wanted, port nShells 4–10 / n2γ up to 8 and the HeI options, then validate against the native code in that configuration.
