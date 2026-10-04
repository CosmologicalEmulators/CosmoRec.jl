# Chunk 5a: the first-pass recombination ODE right-hand side dy/dz of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3):
#   Modules/ODEdef_CosmoRec.cpp:87-135 (copy_LI_to_ysol / copy_ysol_to_LI, compute_fractions :31-47), :898-935 (fcn_effective(int*, ...): f = g / dz_dt, dz_dt = -H (1+z)),
#   :329-380 (fcn_effective(z, LI), with flag_He = 0 skipping fcn_HeI_effective), Modules/main.CosmoRec.cpp:84-262 (Xe_frac_effective_rates).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39) and
# Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK5A_RESULTS.md NOTICE).
#
# STATE (native `ysol`, 1-based): flag_He = 1 (12 entries): [rho = Te/Tg, X(H 1s), X(H 2s, 2p, 3s, 3p, 3d), X(He 1^1S), X(He 2^1S, 2^1P, 2^3S, 2^3P1)]; flag_He = 0 (7 entries): the hydrogen part only.
# It is mapped to the 15-entry `Level_I.X` of Chunk 3d: X[1] = Xe (recomputed, native copy_ysol_to_LI), X[2:7] = H, X[8] = He 1s, X[9, 10, 11, 13] = resolved He, X[12], X[14] = the two
# UNRESOLVED triplet levels (placeholders: verified not to enter g), X[15] = rho. All data are explicit (tables, background accessors, DP fallback); nothing is loaded by default.

"""Resolved hydrogen levels of the production 3-shell atom (native `res_state_list.dat`, as printed by the Chunk 2 harness): `n, l, gw, nuion [Hz], mu_red`."""
const NATIVE_HYDROGEN_RESOLVED = ResolvedStates([2, 2, 3, 3, 3], [0, 1, 0, 1, 2], [2.0, 6.0, 2.0, 6.0, 10.0],
    [822012807922942.25, 822012807922942.25, 365339025743529.88, 365339025743529.88, 365339025743529.88], fill(0.99945567942448077, 5))

"""
    load_native_effective_model(rec_database, fcorr_path; nS_H = 500, nS_He = 30, spin_forbidden = true) -> EffectiveModel

Explicit loader of the production tables from a native `Rec_database` directory (read-only; never bundled): HI `Effective_Rates.HI/Effective_Rate_Tables.nS_3/Rates_n*_l*.nS_500.dat`,
HeI `Effective_Rates.HeI/Effective_Rate_Tables.HeI.res_2/HeI_Rates_*.nS_30.dat` (also the Bitot column), `Pesc_Data/DP_Coll_Data.31.fac_50.neff_30.dat`, and the `f.corr.dat` of `fcorr_path`.
"""
function load_native_effective_model(rec_database::AbstractString, fcorr_path::AbstractString; nS_H::Int = 500, nS_He::Int = 30, spin_forbidden::Bool = true)
    for p in (rec_database, fcorr_path)
        ispath(p) || throw(ArgumentError("native data path does not exist: $p"))
    end
    hdir = joinpath(rec_database, "Effective_Rates.HI", "Effective_Rate_Tables.nS_3")
    hedir = joinpath(rec_database, "Effective_Rates.HeI", "Effective_Rate_Tables.HeI.res_2")
    htable = load_rate_table(hdir, nS_H, NATIVE_HYDROGEN_RESOLVED)
    hetable = load_helium_rate_table(hedir; nS = nS_He)
    dp = read_native_dp_table(joinpath(rec_database, "Pesc_Data", "DP_Coll_Data.31.fac_50.neff_$(nS_He).dat"))
    bitot = read_native_helium_bitot(joinpath(hedir, "HeI_Rates_n2_l1_S0_J1.nS_$(nS_He).dat"))
    fc = read_native_fcorr(fcorr_path)
    return EffectiveModel(htable, hetable, dp, bitot, fc; spin_forbidden = spin_forbidden)
end

"""Everything the ODE right-hand side needs: the 3d `EffectiveModel` (tables), the 4b `CosmosAccessors` (background and loaded Hubble), the native DP fallback
and, in a diffusion iteration, the `HIDiffusionFeedback` of the previous PDE stage (`nothing` = first pass)."""
struct RecombinationModel{M<:EffectiveModel,A<:CosmosAccessors,F,D}
    eff::M
    cosmos::A
    dp_fallback::F
    diffusion::D
end
RecombinationModel(eff, cosmos) = RecombinationModel(eff, cosmos, tr -> dpesc_fallback(NATIVE_DPESC, tr))
RecombinationModel(eff, cosmos, dp_fallback) = RecombinationModel(eff, cosmos, dp_fallback, nothing)
"""The same model with the diffusion feedback `fb` switched on (native `Diffusion_correction_is_on = 1` after `setup_DF_interpol_data`)."""
with_diffusion(rm::RecombinationModel, fb) = RecombinationModel(rm.eff, rm.cosmos, rm.dp_fallback, fb)

const ODE_NRES_H = 5
const ODE_NRES_HE = 4
ode_nstate(flag_He::Bool) = flag_He ? 2 + ODE_NRES_H + 1 + ODE_NRES_HE : 2 + ODE_NRES_H

"""Explicit background at `z` from the 4b accessors: `Tg = TCMB`, `NH`, `H`, `fHe`. `hscale`, `nbscale` (default 1) multiply `H(z)` and `NH(z)` (Chunk 5e differentiable parameters)."""
ode_background(rm::RecombinationModel, z; hscale = 1.0, nbscale = 1.0) =
    EffectiveBackground(cosmos_TCMB(rm.cosmos, z), nbscale * cosmos_NH(rm.cosmos, z), hscale * cosmos_H(rm.cosmos, z), cosmos_fHe(rm.cosmos))

"""Native `copy_ysol_to_LI`: packed `y` -> 15-entry `X` (unresolved He levels and, for `flag_He = 0`, all He slots set to `floor`; `X[1]` = `Xe` from the ground states)."""
const ODE_UNPACK_FLOOR = 1.0e-300
function ode_unpack(y::AbstractVector, fHe, flag_He::Bool; floor = ODE_UNPACK_FLOOR)
    length(y) == ode_nstate(flag_He) || throw(DimensionMismatch("state length $(length(y)) != $(ode_nstate(flag_He)) for flag_He = $flag_He"))
    T = promote_type(eltype(y), typeof(fHe))
    return _ode_unpack!(fill(T(floor), EFF_NEQ), y, fHe, flag_He)
end
# private in-place core (X pre-filled with the floor by the caller; the workspace path refills it); same assignments as before
function _ode_unpack!(X, y, fHe, flag_He::Bool)
    X[15] = y[1]
    for k in 1:(1 + ODE_NRES_H)
        X[1 + k] = y[1 + k]            # X[2] = H 1s, X[3:7] = resolved H
    end
    if flag_He
        X[8] = y[8]
        X[9] = y[9]; X[10] = y[10]; X[11] = y[11]; X[13] = y[12]
        X[1] = (1 - X[2]) + (fHe - X[8])
    else
        X[1] = 1 - X[2]
    end
    return X
end

"""
    recombination_rhs!(f, z, y, rm; flag_He = true) -> f

Native `fcn_effective(int*, ...)` (col < 0): `f = dy/dz = g / dz_dt`, `dz_dt = -H (1+z)`, with `g` the 3d `fcn_effective!` of the unpacked state. `f` has `ode_nstate(flag_He)` entries.
"""
function recombination_rhs!(f::AbstractVector, z, y::AbstractVector, rm::RecombinationModel; flag_He::Bool = true, hscale = 1.0, nbscale = 1.0)
    bg = ode_background(rm, z; hscale = hscale, nbscale = nbscale)
    X = ode_unpack(y, bg.fHe, flag_He)
    T = promote_type(eltype(X), typeof(bg.Tg), typeof(bg.NH), typeof(bg.Hz), typeof(z), feedback_eltype(rm.diffusion))
    g = Vector{T}(undef, EFF_NEQ)
    fcn_effective!(g, z, X, bg, rm.eff; dp_fallback = rm.dp_fallback, flag_He = flag_He, diffusion = rm.diffusion)
    return _pack_rhs!(f, g, bg, z, flag_He)
end
# private: dy/dz from g (shared by the public RHS and the private workspace RHS)
function _pack_rhs!(f, g, bg, z, flag_He::Bool)
    dz_dt = -bg.Hz * (1.0 + z)
    f[1] = g[15] / dz_dt
    for k in 1:(1 + ODE_NRES_H)
        f[1 + k] = g[1 + k] / dz_dt
    end
    if flag_He
        f[8] = g[8] / dz_dt
        f[9] = g[9] / dz_dt; f[10] = g[10] / dz_dt; f[11] = g[11] / dz_dt; f[12] = g[13] / dz_dt
    end
    return f
end

"""Allocating form of [`recombination_rhs!`](@ref)."""
function recombination_rhs(z, y::AbstractVector, rm::RecombinationModel; flag_He::Bool = true, hscale = 1.0, nbscale = 1.0)
    T = promote_type(eltype(y), typeof(z), typeof(hscale), typeof(nbscale), feedback_eltype(rm.diffusion), Float64)
    return recombination_rhs!(Vector{T}(undef, length(y)), z, y, rm; flag_He = flag_He, hscale = hscale, nbscale = nbscale)
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 5b: one native pass (`Xe_frac_effective_rates`, main.CosmoRec.cpp:84-262) on the node grid with the sampled-HeI switch
# ----------------------------------------------------------------------------------------------------------------------------------------------------
"""Native `init_xarr(x0, xm, xarr, npts, 0)` (linear, accumulated like the C++ loop: `xarr[i] = xarr[i-1] + dx`); the last node is NOT exactly `xm` (e.g. 50.000000000188)."""
function init_xarr_linear(x0::Float64, xm::Float64, npts::Int)
    xarr = zeros(npts); dx = (xm - x0) / (npts - 1); xarr[1] = x0
    for i in 2:npts
        xarr[i] = xarr[i - 1] + dx
    end
    return xarr
end

"""ODE parameters passed through SciML: `rm` the explicit model, `flag` the helium flag (`flag_He`). `recombination_ode!` is a plain function (no captured state)."""
struct RecombinationODEParams{R,T}
    rm::R
    flag::Bool
    hscale::T
    nbscale::T
end
RecombinationODEParams(rm, flag::Bool) = RecombinationODEParams(rm, flag, 1.0, 1.0)
recombination_ode!(du, u, p::RecombinationODEParams, z) = recombination_rhs!(du, z, u, p.rm; flag_He = p.flag, hscale = p.hscale, nbscale = p.nbscale)

"""Initial packed 12-state at `zstart` from the 4a Saha initializer fed by the 4b accessors (explicit; native `Set_*_Levels_to_Saha` + `copy_LI_to_ysol`)."""
function ode_initial_state(rm::RecombinationModel, zstart)
    X = saha_initial_state(saha_inputs_at(rm.cosmos, zstart), NATIVE_SAHA_HYDROGEN, NATIVE_SAHA_HELIUM)
    return pack_ysol(X, 6, 7)
end

"""Per-node observables (native `pass_on_the_Solution_CosmoRec` row without z): `[Xe, X1s, 2s, 2p, 3s, 3p, 3d, rho]` from a packed state (`Xe = xp + XHeII`, `XHeII = fHe - He 1s` only while `flag_He = 1`)."""
function ode_observables(y::AbstractVector, fHe, flag_He::Bool)
    xe = flag_He ? (1 - y[2]) + (fHe - y[8]) : 1 - y[2]
    return [xe, y[2], y[3], y[4], y[5], y[6], y[7], y[1]]
end

"""
    recombination_pass(rm, solve_ode; zstart = 3000.0, zend = 50.0, nz = 3000, k_switch = nothing, block = 50) -> NamedTuple

One native pass without diffusion: node grid `init_xarr_linear(zstart, zend, nz)`, initial state ode_initial_state, stiff solves of the 12-state system with the helium test
`(z < 200 || fHe - He1s <= 1e-7)` at every node (4d decision on primal values), the 4d reset at the first node that satisfies it and a restart with the 7-state system.
`solve_ode(f!, u0, z0, znodes, p) -> n_state x length(znodes)` is the caller-supplied stiff solve (`znodes[1] = z0`, descending); `p` is a [`RecombinationODEParams`](@ref).
`k_switch` (1-based node index of the switch, from a previous primal run) freezes the discrete branch (differentiable evaluation, no event detection). Returns the nodes, the 12-state
history up to the switch node, the 7-state history after it, the observables per node (`Xe, X1s, 2s, 2p, 3s, 3p, 3d, rho`) and `Te = rho TCMB`.
The switch node is a DISCRETE quantity: no derivative of its location is defined.
"""
function recombination_pass(rm::RecombinationModel, solve_ode; zstart::Float64 = 3000.0, zend::Float64 = 50.0, nz::Int = 3000, k_switch = nothing, block::Int = 50,
                            hscale = 1.0, nbscale = 1.0, nodes = nothing)
    zfull = init_xarr_linear(zstart, zend, nz)
    if nodes !== nothing                 # Chunk 5e: solve/report only the selected nodes (requires a frozen branch); the switch node and the end nodes are always included
        k_switch === nothing && throw(ArgumentError("recombination_pass: `nodes` requires a frozen `k_switch`"))
        idx = sort(unique(vcat(1, collect(nodes), k_switch, nz)))
        k_sub = findfirst(==(k_switch), idx)
        zn = zfull[idx]; nz = length(idx); k_switch = k_sub; block = nz
    else
        zn = zfull; idx = collect(1:nz)
    end
    fHe = cosmos_fHe(rm.cosmos)
    y0 = ode_initial_state(rm, zstart)          # frozen: depends on the preliminary (4b/4c) history, not differentiated w.r.t. hscale/nbscale
    T = promote_type(eltype(y0), typeof(hscale), typeof(nbscale))
    obs = Vector{Vector{T}}(undef, nz)
    obs[1] = ode_observables(y0, fHe, true)
    ys = Vector{Vector{T}}(undef, nz)              # full packed state at each node (12 entries before the switch node inclusive, 7 after)
    ys[1] = y0
    kswitch = 0
    # ---- 12-state phase ----
    i = 1                                           # current node (1-based), state ys[i] known and its observables stored
    p1 = RecombinationODEParams(rm, true, hscale, nbscale)
    while true
        sw = k_switch === nothing ? helium_switch_condition(zn[i], fHe, ys[i][8]) : (i == k_switch)
        if sw
            kswitch = i
            break
        end
        i == nz && break
        j = min(i + block, nz)
        if k_switch !== nothing
            j = min(j, k_switch)
        end
        sol = solve_ode(recombination_ode!, ys[i], zn[i], zn[i:j], p1)
        for m in (i + 1):j
            ys[m] = collect(sol[:, m - i + 1]); obs[m] = ode_observables(ys[m], fHe, true)
        end
        # the switch test is applied at every node in order (native: start of each interval)
        found = false
        for m in (i + 1):j
            sw = k_switch === nothing ? helium_switch_condition(zn[m], fHe, ys[m][8]) : (m == k_switch)
            if sw
                kswitch = m; found = true; i = m; break
            end
        end
        found && break
        i = j
    end
    # ---- helium reset and 7-state phase ----
    if kswitch > 0 && kswitch < nz
        X = ode_unpack(ys[kswitch], fHe, true)
        y7 = pack_ysol(helium_switch_reset(X, 6, 7, fHe), 6, 7; helium = false)
        ys7 = Vector{Vector{eltype(y7)}}(undef, nz)
        ys7[kswitch] = y7
        p0 = RecombinationODEParams(rm, false, hscale, nbscale)
        i = kswitch
        while i < nz
            j = min(i + block * 4, nz)
            sol = solve_ode(recombination_ode!, ys7[i], zn[i], zn[i:j], p0)
            for m in (i + 1):j
                ys7[m] = collect(sol[:, m - i + 1])
                obs[m] = ode_observables(ys7[m], fHe, false)
            end
            i = j
        end
        # observables are promoted to a common element type
        T2 = promote_type(T, eltype(y7))
        obsT = [Vector{T2}(o) for o in obs]
        ysout = Vector{Vector{T2}}(undef, nz)
        for m in 1:kswitch; ysout[m] = Vector{T2}(ys[m]); end
        for m in (kswitch + 1):nz; ysout[m] = Vector{T2}(ys7[m]); end
        obs = obsT; ys = ysout
    else
        kswitch == 0 && (kswitch = nz + 1)
    end
    Xe = [o[1] for o in obs]; rho = [o[8] for o in obs]
    Te = [rho[m] * cosmos_TCMB(rm.cosmos, zn[m]) for m in 1:nz]
    return (z = zn, node_index = idx, k_switch = kswitch, states = ys, observables = obs, Xe = Xe, rho = rho, Te = Te, zswitch = kswitch <= nz ? zn[kswitch] : NaN)
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 5c: the Recfast tail (CosmoRec.cpp:278-324 compute_Recfast_part; Recfast++.cpp:336-441 Xe_frac_rescaled; Cosmos.cpp:795-840 recombine_using_Recfast_system)
# ----------------------------------------------------------------------------------------------------------------------------------------------------
"""Native `init_xarr(x0, xm, xarr, npts, 1)` (logarithmic, accumulated: `xarr[i] = xarr[i-1] dx`, `dx = (xm/x0)^(1/(npts-1))`)."""
function init_xarr_log(x0::Float64, xm::Float64, npts::Int)
    xarr = zeros(npts); dx = (xm / x0)^(1.0 / (npts - 1)); xarr[1] = x0
    for i in 2:npts
        xarr[i] = xarr[i - 1] * dx
    end
    return xarr
end

"""ODE parameters of the rescaled Recfast system: `θ` (4c parameter vector), constants, Hubble function and the rescaling factor `ff` of the hydrogen equation (native `fcn_rescaled`)."""
struct RecfastTailParams{T,C,H}
    θ::Vector{T}
    c::C
    Hfun::H
    ff::T
end
"""Plain 3-state right-hand side `(xHep, xp, TM)` of native `fcn_rescaled`: the Recfast system with `f[1] *= ff` (`f[1]` is the `xp` equation)."""
function recfast_tail_rhs!(du, u, p, z)
    f = recfast_rhs(p.θ, z, u, p.Hfun(z), p.c)
    du[1] = f[1]; du[2] = f[2] * p.ff; du[3] = f[3]
    return du
end

"""
    recfast_tail(θ, c, Hfun, zi, Xe_Hi, Xe_Hei, Xei, TMi, dXei, solve_ode; nz = 200, ze = 0.001) -> NamedTuple(z, Xe_H, Xe_He, Xe, TM, ff)

Native `Xe_frac_rescaled` as called by `compute_Recfast_part`: nodes `init_xarr_log(zi, max(1e-5, ze), nz)`, initial state `(Xe_He = Xe_Hei, Xe_H = Xe_Hi, TM = TMi)`, rescaling factor
`ff = dXei / f_xp(zi, y0)` (the UNscaled Recfast `xp` equation), stiff solve of the rescaled system via the caller's `solve_ode(f!, u0, z0, znodes, p)`. `Xe = Xe_H + Xe_He`. Row 1 is the initial state.
"""
function recfast_tail(θ, c::RecfastConstants, Hfun, zi, Xe_Hi, Xe_Hei, Xei, TMi, dXei, solve_ode; nz::Int = 200, ze::Float64 = 0.001)
    zarr = init_xarr_log(Float64(zi), max(1.0e-5, ze), nz)
    y0 = [Xe_Hei, Xe_Hi, TMi]
    f0 = recfast_rhs(θ, zarr[1], (y0[1], y0[2], y0[3], 0.0), Hfun(zarr[1]), c)
    ff = dXei / f0[2]
    Tp = promote_type(eltype(θ), typeof(ff))
    p = RecfastTailParams(Vector{Tp}(θ), c, Hfun, Tp(ff))
    sol = solve_ode(recfast_tail_rhs!, y0, zarr[1], zarr, p)
    Xe_He = vcat(Xe_Hei, [sol[1, k] for k in 2:nz]); Xe_H = vcat(Xe_Hi, [sol[2, k] for k in 2:nz]); TM = vcat(TMi, [sol[3, k] for k in 2:nz])
    return (z = zarr, Xe_H = Xe_H, Xe_He = Xe_He, Xe = vcat(Xei, Xe_H[2:end] .+ Xe_He[2:end]), TM = TM, ff = ff)
end

"""Native inputs of `compute_Recfast_part` from the last ODE node (flag_He = 0): `Xe_Hi = 1 - X1s`, `Xe_Hei = 1e-30` (`max(1e-30, 0)`), `Xei = Xe = xp`, `TMi = TCMB rho`, `dXei = -dX1s/dz` (native `-g[iHI]/dz_dt`)."""
function recfast_tail_inputs(rm::RecombinationModel, zi, y7::AbstractVector; hscale = 1.0, nbscale = 1.0)
    f = recombination_rhs(zi, y7, rm; flag_He = false, hscale = hscale, nbscale = nbscale)
    xp = 1 - y7[2]
    return (Xe_Hi = xp, Xe_Hei = 1.0e-30, Xei = xp, TMi = cosmos_TCMB(rm.cosmos, zi) * y7[1], dXei = -f[2])
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 5d: output assembly (CosmoRec.cpp:278-324 output rows, :545-597 return_solution_to_calling_program) and the full single pass
# ----------------------------------------------------------------------------------------------------------------------------------------------------
"""Native `output_CosmoRec` rows (z, Xe, Te [K]): the ODE nodes (3000) followed by the Recfast tail rows 2..200 (the initial tail row is not stored natively)."""
function recombination_output_rows(pass, tail)
    z = vcat(pass.z, tail.z[2:end])
    Xe = vcat(pass.Xe, tail.Xe[2:end])
    Te = vcat(pass.Te, tail.TM[2:end])
    return (z = z, Xe = Xe, Te = Te)
end

"""
    return_solution_to_grid(rm, rows, zgrid) -> (Xe, Te)

Native `return_solution_to_calling_program`: natural cubic splines (GSL, 4b) of `ln Xe` and `ln Te` against ascending z through the stored rows, evaluated with the native end nudge for
`zmin < z < zmax` of the stored rows; at or beyond the ends the Recfast preliminary history values `Xe_Seager(z)` and `Te(z)` (4b accessors) are returned instead.
"""
function return_solution_to_grid(rm::RecombinationModel, rows, zgrid::AbstractVector)
    za = reverse(rows.z); n = length(za)
    sXe = natural_cubic_spline(za, log.(reverse(rows.Xe)))
    sTe = natural_cubic_spline(za, log.(reverse(rows.Te)))
    zmin = min(za[1], za[n]); zmax = max(za[1], za[n])
    T = promote_type(eltype(rows.Xe), eltype(rows.Te), Float64)
    Xe = Vector{T}(undef, length(zgrid)); Te = Vector{T}(undef, length(zgrid))
    for (i, zq) in enumerate(zgrid)
        if zq >= zmax || zq <= zmin
            Xe[i] = cosmos_Xe_Seager(rm.cosmos, zq); Te[i] = cosmos_Te(rm.cosmos, zq)
        else
            Xe[i] = exp(spline_eval_native(sXe, zq)); Te[i] = exp(spline_eval_native(sTe, zq))
        end
    end
    return Xe, Te
end

"""
    recombination_history(rm, θ, solve_pass, solve_tail, zgrid; kwargs...) -> NamedTuple

The complete native single pass (iteration 0, no diffusion; native `runmode 1`): ODE pass (5b), Recfast tail (5c), output rows and the returned `Xe`, `Te` on `zgrid` (5d).
`θ` is the 4c Recfast parameter vector `[F, A2s1s, Yp, Omega_b, h100, T0]`; `hscale`, `nbscale` multiply `H(z)` and the baryon density `NH`/`Omega_b` consistently in the pass and the tail (Chunk 5e parameters; the initial state and the 4b preliminary history are frozen); `solve_pass` / `solve_tail` are the caller's stiff solves. Other keyword arguments go to [`recombination_pass`](@ref).
"""
function recombination_history(rm::RecombinationModel, θ, solve_pass, solve_tail, zgrid::AbstractVector; hscale = 1.0, nbscale = 1.0, kwargs...)
    pass = recombination_pass(rm, solve_pass; hscale = hscale, nbscale = nbscale, kwargs...)
    zi = pass.z[end]
    y7 = pass.states[end]
    inp = recfast_tail_inputs(rm, zi, y7; hscale = hscale, nbscale = nbscale)
    Hfun = z -> hscale * cosmos_H(rm.cosmos, z)
    θeff = [θ[1], θ[2], θ[3], θ[4] * nbscale, θ[5], θ[6]]       # Recfast NH is proportional to Omega_b
    tail = recfast_tail(θeff, NATIVE_RECFAST_CONSTANTS, Hfun, zi, inp.Xe_Hi, inp.Xe_Hei, inp.Xei, inp.TMi, inp.dXei, solve_tail)
    rows = recombination_output_rows(pass, tail)
    Xe, Te = return_solution_to_grid(rm, rows, zgrid)
    return (pass = pass, tail = tail, rows = rows, Xe = Xe, Te = Te, z = zgrid)
end
