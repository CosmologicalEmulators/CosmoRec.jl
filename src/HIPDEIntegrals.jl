# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 8: correction integrals over the HI PDE spectrum (fed back into the next recombination ODE pass).
# Sources (ORIGINAL CosmoRec v3.0b, read-only): PDE_Problem/Solve_PDEs_integrals.cpp (tau_S_function, calc_DF, compute_integral_over_resonances,
# compute_DF_2gamma, compute_DF_Raman, compute_DI1_2s_and_dump_it with DIFF_CORR_STOREII undefined, epsrel_HI = 1e-5), PDE_Problem/Solve_PDEs.cpp:840-966
# (outputs for zout < zs - 10, DF_step = 1), Development/Recombination/Sobolev.cpp (p_ij), Development/Integration/Patterson.cpp
# (Integrate_using_Patterson_adaptive; its refinement branch is unreachable, see `patterson_integrate`).
# Production outputs: DI1_2s (2s-1s two-photon), DF_2gamma for 3s and 3d, DF_Raman for 2s. Ly-n and nD-1s integrals are off in production.
# ----------------------------------------------------------------------------------------------------------------------------------------------------

const HI_EPSREL_DF = 1.0e-5

"""Native Sobolev escape probability `p_ij(tau_S)` (Sobolev.cpp:73-78)."""
function sobolev_p_ij(tau)
    tau <= 1.0e-10 && return 1.0 - 0.5 * tau + tau * tau / 6.0
    tau >= 5.0e2 && return 1.0 / tau
    return (1.0 - exp(-tau)) / tau
end

"""Native `tau_S_function(n, z, X1s)` = 3 A21 lambda21^3 X1s NH / (8 pi H) for the Ly-n line."""
hi_tau_S(m::HIPDEModel, n::Int, z, X1s) = 3.0 * m.atom.lyn_A21[n - 1] * m.atom.lyn_lambda21[n - 1]^3 * X1s * cosmos_NH(m.cosmos, z) / 8.0 / PDE_PI / cosmos_H(m.cosmos, z)

# calc_DF on a native cspline of Fwork
"""Record/replay of the Patterson stopping decisions of a PDE stage (one entry per sub-integral, in evaluation order): `mode = :record` stores the
native decisions, `mode = :replay` applies the stored rule without tests (frozen branch, for finite-difference verification of derivatives, which by
construction differentiate the frozen-decision map)."""
mutable struct PattersonLevels
    mode::Symbol
    levels::Vector{Int}
    pos::Int
end
PattersonLevels() = PattersonLevels(:record, Int[], 0)
replay(p::PattersonLevels) = PattersonLevels(:replay, p.levels, 0)

# returns (value, converged): `converged = false` when the native stopping rule |r_k - r_{k-1}| <= max(epsabs, |r_k| epsrel) was never met
# (the 255-point value is returned without refinement, as natively)
function _hi_calc_DF(s, a, b, epsabs, rec::Union{Nothing,PattersonLevels} = nothing)
    a >= b && return (zero(eltype(s.y)), true)
    f = t -> spline_eval_native(s, t)
    if rec !== nothing && rec.mode === :replay
        rec.pos += 1
        r, _ = patterson_integrate(f, a, b, HI_EPSREL_DF, epsabs; level = rec.levels[rec.pos])
        return (r, true)
    end
    r, lev = patterson_integrate(f, a, b, HI_EPSREL_DF, epsabs)
    rec === nothing || push!(rec.levels, lev)
    lev < 8 && return (r, true)
    r7, _ = patterson_integrate(f, a, b, HI_EPSREL_DF, epsabs; level = 7)
    return (r, abs(_primal(r) - _primal(r7)) <= max(epsabs, abs(_primal(r)) * HI_EPSREL_DF))
end

# native interval boundaries: xmin, xmin(1+0.005), then x_res(n)(1 -/+ dx) for n = nmin..nmax, xmax(1-dx), xmax
# (shared by the default path and the fixed-geometry quadrature plans of HIPDEQuadPlan.jl)
function _hi_resonance_edges(atom::HIPDEAtom, xmin, xmax, nmin::Int, nmax::Int)
    nu21 = atom.nu21
    dx_x = 1.0e-4
    edges = [xmin, xmin * (1.0 + 0.005)]
    for n in nmin:nmax
        xres = atom.lyn[n - 1].nu21 / nu21
        push!(edges, xres * (1.0 - dx_x), xres * (1.0 + dx_x))
    end
    push!(edges, xmax * (1.0 - dx_x), xmax)
    return edges
end

function _hi_integral_over_resonances(m::HIPDEModel, s, xmin, xmax, nmin::Int, nmax::Int, rec = nothing)
    edges = _hi_resonance_edges(m.atom, xmin, xmax, nmin, nmax)
    epsabs = 1.0e-60
    r = zero(eltype(s.y)); sabs = 0.0; nunc = 0
    for i in 1:(length(edges) - 1)
        r1, ok = _hi_calc_DF(s, edges[i], edges[i + 1], epsabs, rec)
        r += r1; sabs += abs(_primal(r1)); nunc += !ok
        epsabs = abs(_primal(r)) * HI_EPSREL_DF
    end
    return r, sabs, nunc   # sabs = sum of |sub-integrals|; nunc = number of sub-integrals that never met the native stopping rule
end

# integration limits of the production blocks (shared by the default path and the quadrature plans; expressions unchanged)
_hi_limits_di1(x, nnu) = (max(x[1], 0.5), min(x[nnu], 1.0))
_hi_limits_two3(x, nnu, a::HIPDEAtom) = (nui1 = a.Dnu_1s[3]; nu21 = a.nu21; (max(0.5 * nui1 / nu21, x[1]), min(nui1 / nu21, x[nnu - 1])))
_hi_limits_raman(x, k0, nnu, a::HIPDEAtom) = (nui1 = a.Dnu_1s[2]; nu21 = a.nu21; (max(nui1 / nu21, x[k0 + 1]), x[nnu + k0 - 1]))

"""
    hi_pde_integrals(m, z, y, st)

The production correction integrals at an output redshift `z` from the spectrum `y` (= nu21 Dn(x)) and the coefficient side products `st` evaluated at
`z` (native `PDE_funcs` after the step): `values = (DI1_2s, DF_2g_3s, DF_2g_3d, DF_R_2s)` in the native order of Solve_PDEs.cpp:868-924, and
`scale`, the source-derived residual scale of each (every output is a DIFFERENCE of an integral over the spectrum and a Sobolev/equilibrium term,
`A2s1s(I/A2s1s + DnL)` or `w(r - DRtot)`, which cancel near equilibrium: scale = the same combination with absolute values, where for the
resonance-split integrals `|r|` is replaced by the sum of the absolute sub-integrals, because the integrand changes sign and the native Patterson
tolerances (epsrel_HI = 1e-5) apply to the individual sub-integrals).
"""
function hi_pde_integrals(m::HIPDEModel, z, y::AbstractVector, st::HIPDEState; levels = nothing, quadplan = nothing)
    s = m.setup; a = m.atom; k = m.k; x = s.x; np = length(x)
    nu21 = a.nu21
    Tg = cosmos_TCMB(m.cosmos, z)
    X1s = hi_Xi(m.pops, z, 0)
    T = promote_type(eltype(y), typeof(z), eltype(st.exp_x))
    # DI1 2s-1s (compute_DI1_2s_and_dump_it)
    DI1 = let nnu = s.index_emission[1], nu2s = a.Dnu_1s[2]
        exp_i1 = exp(-k.h_kb * nu2s / Tg)
        Dnem = _pde_Dnem(m, z, 2, 0)
        F = Vector{T}(undef, nnu)
        for i in 1:nnu
            f_x = exp_i1 * (1.0 / st.exp_x[i] - 1.0)
            F[i] = (y[i] * f_x - nu2s * Dnem) * st.P2g_ns[1][i]
        end
        if quadplan === nothing
            sp = natural_cubic_spline(x[1:nnu], F)
            lim = _hi_limits_di1(x, nnu)
            r, ok = _hi_calc_DF(sp, lim[1], lim[2], 1.0e-60, levels)
        else
            r, ok = _hi_plan_single(quadplan.di1, F, levels)
        end
        r = r / PROF_A2S1S(m)
        DnLj = _pde_Xnl(m, z, 2, 0) / X1s - exp_i1
        (PROF_A2S1S(m) * (r + DnLj), PROF_A2S1S(m) * (abs(_primal(r)) + abs(_primal(DnLj))), Int(!ok))
    end
    # two-photon 3s / 3d (compute_DF_2gamma(3, li))
    DF2g = map(((0, st.P2g_ns[2]), (2, st.P2g_nd[2]))) do (li, P)
        ni = 3; nnu = s.index_emission[ni - 1]
        nui1 = a.Dnu_1s[ni]
        w = 2.0 * li + 1
        exp_i1 = exp(-k.h_kb * nui1 / Tg)
        DnLj = _pde_Xnl(m, z, ni, li) / X1s / w - exp_i1
        Dnem = _pde_Dnem(m, z, ni, li)
        DRtot = zero(T)
        for n in 2:(ni - 1)
            nuik = li == 0 ? a.nu_3s2p : a.nu_3d2p
            Aik = li == 0 ? a.A_3s2p : a.A_3d2p
            exp_n = exp(-k.h_kb * a.lyn[n - 1].nu21 / Tg)
            tau_S_n = hi_tau_S(m, n, z, X1s)
            pem = 1.0 - _pde_pd(m, z, n, 1)
            exp_ik = exp(-k.h_kb * nuik / Tg)
            PS = sobolev_p_ij(tau_S_n)
            Xkp = _pde_Xnl(m, z, n, 1)
            Atilde = Aik * pem / (1.0 - exp_ik)
            DnLkp = Xkp / X1s / 3.0 - exp_n
            Del = DnLkp - DnLj / exp_ik
            DRtot += Atilde * exp_ik * (Del - PS * DnLkp)
        end
        F = Vector{T}(undef, nnu)
        for i in 1:nnu
            f_x = exp_i1 * (1.0 / st.exp_x[i] - 1.0)
            F[i] = (y[i] * f_x - nu21 * Dnem) * P[i]
        end
        r, sabs, nunc = if quadplan === nothing
            sp = natural_cubic_spline(x[1:nnu], F)
            xmin, xmax = _hi_limits_two3(x, nnu, a)
            _hi_integral_over_resonances(m, sp, xmin, xmax, 2, ni - 1, levels)
        else
            _hi_plan_resonances(quadplan.two3, F, levels)
        end
        (w * (r - DRtot), w * (sabs + abs(_primal(DRtot))), nunc)
    end
    # Raman 2s (compute_DF_Raman(2, 0, nmax = nresmax = 3))
    DFR = let ni = 2, li = 0, nmax = s.nresmax, k0 = s.index_emission[1]
        nnu = np - k0
        nui1 = a.Dnu_1s[ni]
        w = 2.0 * li + 1.0
        exp_i1 = exp(-k.h_kb * nui1 / Tg)
        DnLj = _pde_Xnl(m, z, ni, li) / X1s / w - exp_i1
        Dnem = _pde_Dnem(m, z, ni, li)
        DRtot = zero(T)
        for n in (ni + 1):nmax
            nuik = a.nu_3p2s; Aik = a.A_3p2s
            exp_n = exp(-k.h_kb * a.lyn[n - 1].nu21 / Tg)
            tau_S_n = hi_tau_S(m, n, z, X1s)
            pem = 1.0 - _pde_pd(m, z, n, 1)
            exp_ik = exp(-k.h_kb * nuik / Tg)
            PS = sobolev_p_ij(tau_S_n)
            Xkp = _pde_Xnl(m, z, n, 1)
            Atilde = 3.0 / w * Aik * pem / (1.0 - exp_ik)
            DnLkp = Xkp / X1s / 3.0 - exp_n
            Del = DnLkp - DnLj * exp_ik
            DRtot += Atilde * (Del - PS * DnLkp)
        end
        F = Vector{T}(undef, nnu)
        for i in 1:nnu
            f_x = exp_i1 * (1.0 / st.exp_x[i + k0] - 1.0)
            F[i] = (y[i + k0] * f_x - nu21 * Dnem) * st.PR_ns[1][i + k0]
        end
        r, sabs, nunc = if quadplan === nothing
            sp = natural_cubic_spline(x[(k0 + 1):np], F)
            xmin, xmax = _hi_limits_raman(x, k0, nnu, a)
            _hi_integral_over_resonances(m, sp, xmin, xmax, ni + 1, nmax, levels)
        else
            _hi_plan_resonances(quadplan.raman, F, levels)
        end
        (w * (r - DRtot), w * (sabs + abs(_primal(DRtot))), nunc)
    end
    return (values = (DI1[1], DF2g[1][1], DF2g[2][1], DFR[1]), scale = (DI1[2], DF2g[1][2], DF2g[2][2], DFR[2]),
            unconverged = (DI1[3], DF2g[1][3], DF2g[2][3], DFR[3]))
end
PROF_A2S1S(m::HIPDEModel) = m.setup.prof.A2s1s

"""
    hi_pde_corrections(m; zs = 2500.0, ze = 500.0, kwargs...)

The production PDE stage of one diffusion iteration (`compute_DPesc_with_diffusion_equation_effective`): march the spectrum from `zs` to `ze` and evaluate
the correction integrals at every output (`zout < zs - 10`). Returns `(z, DI1_2s, DF_2g = [3s, 3d], DF_R = [2s], scale, y)` (`scale`: residual scales of\nthe four outputs, see `hi_pde_integrals`; `y`: final spectrum).
"""
function hi_pde_corrections(m::HIPDEModel; zs = 2500.0, ze = 500.0, T::Type = Float64, levels = nothing, quadplan = nothing, kwargs...)
    quadplan === nothing || _check_quadplan(quadplan, m.setup, m.atom)    # full-grid stale-plan guard, before the march
    zo = Float64[]; vals = [T[] for _ in 1:4]; scs = [Float64[] for _ in 1:4]; unc = [Int[] for _ in 1:4]
    cb = (k, zin, zout, y, S) -> begin
        if zout < zs - 10.0
            r = hi_pde_integrals(m, zout, y, S.prev; levels = levels, quadplan = quadplan)
            push!(zo, zout)
            for q in 1:4
                push!(vals[q], r.values[q]); push!(scs[q], _primal(r.scale[q])); push!(unc[q], r.unconverged[q])
            end
        end
    end
    y, _ = hi_pde_march(m; zs = zs, ze = ze, on_step = cb, T = T, kwargs...)
    return (z = zo, DI1_2s = vals[1], DF_2g = [vals[2], vals[3]], DF_R = [vals[4]], scale = Tuple(scs), unconverged = Tuple(unc), y = y)
end
