# Chunk 4c: preliminary Recfast++ history of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), production settings only:
#   Development/Cosmology/Recfast++/Recfast++.cpp:56-95 (SahaBoltz_HeIII/HeII, Te_QS), :119-335 (Xe_frac: Saha segments, grid, ODE start),
#   src/evalode.Recfast.cpp:40-260 (evaluate_Recfast_System and rate functions), src/cosmology.Recfast.cpp (NH, H pointer, fHe, t_cos_rad),
#   include/constants.Recfast.h (SI constants), src/Variation_constants.Recfast.cpp (factors == 1 for aS = mS = 1), src/ODE_solver.Recfast.cpp (native solver: variable-order Gear, rel tol 1e-8).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39) and
# Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK4C_RESULTS.md NOTICE).
#
# ACTIVE TERMS ONLY (verified at runtime from the native object, harness record RFV): `F` fudge, `A2s1s`, no DM annihilation/decay, no magnetic heating, no reionization terms,
# no Chluba-Thomas correction function (include_CF = 0), eval_ion_Tr = 0 (photoionization at T = Tm), varying-constant factors = 1 (aS = mS = 1, pS = 0).
# The package is solver-agnostic (no ODE dependency in `[deps]`): `recfast_history` takes the stiff solve as a callback; tests/benchmarks supply SciML (`Rodas5P`) solves.
# Everything in SI (as the native Recfast++): `NH` in m^-3, rates in m^3/s. The native Hubble function is the loaded Cosmos H(z) (a callback `Hfun(z)` in 1/s).
# Native column meaning of the history: Xe_H = n_e,H/n_H and Xe_He = n_e,He/n_H are IONIZED populations (not neutral fractions).

"""Native `RECFAST_physical_constants` / `RECFAST_atomic_constants` (constants.Recfast.h), SI."""
struct RecfastConstants
    cLight::Float64
    hPlanck::Float64
    kBoltz::Float64
    mElect::Float64
    mHatom::Float64
    G::Float64
    aRad::Float64
    sigmaT::Float64
    fac_mHemH::Float64
    Mpc::Float64
    Lam2s1sHe::Float64
    LyalphaH::Float64
    LyalphaHe::Float64
    EionH2s::Float64
    E2s1sH::Float64
    EionHe2s::Float64
    EionHeII::Float64
    E2s1sHe::Float64
    E2p2sHe::Float64
    pi::Float64
end

function RecfastConstants()
    cLight = 2.99792458E+08; hPlanck = 6.62606876E-34; kBoltz = 1.3806503E-23; mElect = 9.10938188E-31
    amu = 1.660538782e-24; mH_gr = 1.0078250321 * amu
    L_H_ion = 1.096787737e+7; L_H_alpha = 8.225916453e+6
    return RecfastConstants(cLight, hPlanck, kBoltz, mElect, mH_gr * 1.0e-3, 6.67428e-8 * 1.0e-3, 4.0 * 5.670400e-8 / cLight, 6.6524616e-29, 3.97152594 / 4.0, 3.08568025e+22,
                            51.3, 1.215670e-07, 5.843344e-08, (L_H_ion - L_H_alpha) * cLight * hPlanck, L_H_alpha * cLight * hPlanck, 6.363254e-19, 8.7186944e-18, 3.30301387e-18, 9.64908313e-20,
                            3.1415926535897932)
end
const NATIVE_RECFAST_CONSTANTS = RecfastConstants()

"""
Differentiable physical parameters `θ = [F, A2s1s, Yp, Omega_b, h100, T0]`: Recfast fudge factor, H 2s-1s rate [1/s] (native default 8.22458; the production run reads 8.2206),
helium mass fraction, baryon density, `h100` (used only by `NH`), CMB temperature [K].
"""
const RECFAST_THETA_NAMES = ("F", "A2s1s", "Yp", "Omega_b", "h100", "T0")

rf_fHe(c::RecfastConstants, Yp) = Yp / (4.0 * c.fac_mHemH * (1.0 - Yp))
rf_H0(c::RecfastConstants, h100) = h100 * 100.0 * 1.0e+3 / c.Mpc
rf_NH(c::RecfastConstants, θ, z) = 3.0 * rf_H0(c, θ[5])^2 * θ[4] / (8.0 * c.pi * c.G * c.mHatom * (1.0 / (1.0 - θ[3]))) * (1.0 + z)^3
rf_TCMB(θ, z) = θ[6] * (1.0 + z)

rf_boltzmann(gj, gi, E, T, c::RecfastConstants) = max(1.0e-300, (gj / gi) * exp(-E / (c.kBoltz * T)))
# native `1.0/Boltzmann(1,1,E,T)` = 1/max(1e-300, exp(-E/(kB T))) (<= 1e300); written as min(1e150, exp(+E/(kB T))) so that the clamped branch has a finite (zero) derivative
# (1/Dual(1e-300, 0) overflows to Inf*0 = NaN under AD); equal to the native value to rounding, and equal in the clamped regime up to 1e-16 (insensitive there: it only enters C_He).
rf_inv_boltzmann(E, T, c::RecfastConstants) = min(1.0e+150, exp(E / (c.kBoltz * T)))
function rf_sahaboltz(gi, gc, ge, E_ion, T, c::RecfastConstants)
    c1 = (c.hPlanck / (2.0 * c.pi * c.mElect)) * (c.hPlanck / c.kBoltz)
    return (ge * gi / (2.0 * gc) * c1^1.5 * T^(-1.5) * exp(min(680.0, E_ion / (c.kBoltz * T))))
end

# 1/rf_sahaboltz written without forming the huge Saha factor (exp(min(680, .)) ~ 1e295 at low Tm): the native `alpha/SahaBoltz` quotient is mathematically identical, but its
# second-order AD partials overflow (Inf/Inf = NaN inside the solver's nested Jacobian); equal to the native value to rounding.
function rf_inv_sahaboltz(gi, gc, ge, E_ion, T, c::RecfastConstants)
    c1 = (c.hPlanck / (2.0 * c.pi * c.mElect)) * (c.hPlanck / c.kBoltz)
    return (2.0 * gc / (ge * gi)) * c1^(-1.5) * T^1.5 * exp(-min(680.0, E_ion / (c.kBoltz * T)))
end

function rf_alphaH(F, TM)
    a1 = 4.309; a2 = -0.6166; a3 = 0.6703; a4 = 0.5300; t = TM / 1.0e+4
    return F * a1 * 1.0e-19 * t^a2 / (1.0 + a3 * t^a4)
end
function rf_alphaHe(TM)
    a1 = 10.0^(-16.744); a2 = 0.711; T0 = 3.0; T1 = 10.0^5.114
    sq0 = sqrt(TM / T0); sq1 = sqrt(TM / T1)
    return a1 / (sq0 * (1.0 + sq0)^(1.0 - a2) * (1.0 + sq1)^(1.0 + a2))
end

"""
    recfast_rhs(θ, z, y, Hz, c = NATIVE_RECFAST_CONSTANTS) -> (f1, f2, f3, f4)

Native `evaluate_Recfast_System` for the active production terms: `y = (xHep, xp, TM, t_cos)`, `Hz` = Hubble [1/s] (the native H pointer = Cosmos H). Independent variable `z`
(derivatives w.r.t. z, so integrate from large to small z). `f4 = -1/(Hz (1+z))` is the internal cosmic-time equation.
"""
function recfast_rhs(θ, z, y, Hz, c::RecfastConstants = NATIVE_RECFAST_CONSTANTS)
    F, A2s1s, Yp = θ[1], θ[2], θ[3]
    nHTot = rf_NH(c, θ, z)
    fHe = rf_fHe(c, Yp)
    xHep = y[1]; xp = y[2]; TM = y[3]
    xe = xp + xHep
    nH1 = (1.0 - xp) * nHTot
    nHe1 = nHTot * (fHe - xHep)
    TR = rf_TCMB(θ, z)
    Comp = 8.0 * c.sigmaT * c.aRad * TR^4 / (3.0 * Hz * (1.0 + z) * c.mElect * c.cLight)
    A2s1sH = A2s1s; A2s1sHe = c.Lam2s1sHe
    lambda21H = c.LyalphaH; lambda21He = c.LyalphaHe
    alphaH = rf_alphaH(F, TM)
    alphaHe = rf_alphaHe(TM)
    BH = rf_boltzmann(1.0, 1.0, c.E2s1sH, TM, c)
    betaH = rf_alphaH(F, TM) * rf_inv_sahaboltz(2.0, 1.0, 1.0, c.EionH2s, TM, c)
    BHe = rf_boltzmann(1.0, 1.0, c.E2s1sHe, TM, c)
    BHe2p2s = rf_inv_boltzmann(c.E2p2sHe, TM, c)
    betaHe = rf_alphaHe(TM) * rf_inv_sahaboltz(1.0, 2.0, 1.0, c.EionHe2s, TM, c)
    KHe = lambda21He^3 / (Hz * 8.0 * c.pi)
    KH = lambda21H^3 / (Hz * 8.0 * c.pi)
    CHe = (1.0 + KHe * A2s1sHe * nHe1 * BHe2p2s) / (1.0 + KHe * (A2s1sHe + betaHe) * nHe1 * BHe2p2s)
    CH = (1.0 + KH * A2s1sH * nH1) / (1.0 + KH * (A2s1sH + betaH) * nH1)
    f1 = (alphaHe * xe * xHep * nHTot - betaHe * (fHe - xHep) * BHe) * CHe / ((1.0 + z) * Hz)
    f2 = (alphaH * xe * xp * nHTot - betaH * (1.0 - xp) * BH) * CH / ((1.0 + z) * Hz)
    f3 = Comp * xe / (1.0 + xe + fHe) * (TM - TR) + 2.0 * TM / (1.0 + z)
    f4 = -1.0 / (Hz * (1.0 + z))
    return (f1, f2, f3, f4)
end

"""In-place 3-state SciML right-hand side `(xHep, xp, TM)` (`t_cos` is not used by the active terms); `p = (θ, c, Hfun)`; independent variable `z`."""
function recfast_rhs3!(du, u, p, z)
    θ, c, Hfun = p
    f = recfast_rhs(θ, z, u, Hfun(z), c)
    du[1] = f[1]; du[2] = f[2]; du[3] = f[3]
    return du
end

"""Native `Te_QS`: quasi-stationary matter temperature `TR (1 - 1/(1+kappa))`, `kappa = Comp fC Xe/(1+Xe+fHe)`, `Comp = 8 sigmaT aRad TR^4/(3 H me c)`."""
function rf_Te_QS(θ, z, Xe, Hz, c::RecfastConstants = NATIVE_RECFAST_CONSTANTS, fC = 1.0)
    TR = rf_TCMB(θ, z)
    Comp = 8.0 * c.sigmaT * c.aRad * TR^4 / (3.0 * Hz * c.mElect * c.cLight)
    kappa = Comp * fC * Xe / (1.0 + Xe + rf_fHe(c, θ[3]))
    return TR * (1.0 - 1.0 / (1.0 + kappa))
end

"""Native Recfast++ `SahaBoltz_HeIII(nH, fHe, T)` (nH in m^-3)."""
function rf_SahaBoltz_HeIII(nH, fHe, T, c::RecfastConstants = NATIVE_RECFAST_CONSTANTS)
    c1 = 2.0 * c.pi * c.mElect * c.kBoltz * T / (c.hPlanck * c.hPlanck)
    g = c1^1.5 * exp(-c.EionHeII / (c.kBoltz * T)) / nH; d = 1.0 + fHe - g
    return 0.5 * (d + sqrt(d * d + 4.0 * g * (1.0 + 2.0 * fHe)))
end
"""Native Recfast++ `SahaBoltz_HeII(nH, fHe, T)`. Uses `EionHeII` in the exponent exactly as the native code does (it is the native source; not altered)."""
function rf_SahaBoltz_HeII(nH, fHe, T, c::RecfastConstants = NATIVE_RECFAST_CONSTANTS)
    c1 = 2.0 * c.pi * c.mElect * c.kBoltz * T / (c.hPlanck * c.hPlanck)
    g = 4.0 * c1^1.5 * exp(-c.EionHeII / (c.kBoltz * T)) / nH; d = 1.0 - g
    return 0.5 * (d + sqrt(d * d + 4.0 * g * (1.0 + fHe)))
end

"""Node grid of one native `Xe_frac` call: `z` (descending, native order), the 1-based end index of each Saha segment, and the index `jode` of the first ODE node (initial state)."""
struct RecfastGrid
    z::Vector{Float64}
    seg_end::NTuple{3,Int}
    jode::Int
end

"""
    recfast_grid(θ, c, Hfun; npts = 6000, zstart = 2.5e4, zend = 0.0) -> RecfastGrid

Native `Xe_frac` grid logic (Recfast++.cpp:119-335) with the data-dependent segment breaks (`|Xe/Xe_Saha - 1| >= 5e-4` etc., decided on PRIMAL values: `θ` is converted to Float64) -
a discrete operator: the nodes depend on `θ` only through these decisions and `T0`.
"""
function recfast_grid(θ_in, c::RecfastConstants, Hfun; npts::Int = 6000, zstart::Float64 = 2.5e4, zend_rec::Float64 = 0.0)
    θ = Float64.(_primal.(θ_in))
    zarr = zeros(npts)
    fHe = rf_fHe(c, θ[3])
    Saha_eps_I = 5.0e-4; Saha_eps_II = 1.0e-5
    n_low_res = 100; n_high_res = (npts - 2 * n_low_res) ÷ 3
    zc_HeIII = 9000.0 * 2.726 / θ[6]; zc_HeIIIS = 4800.0 * 2.726 / θ[6]; zc_HeII = 3600.0 * 2.726 / θ[6]
    ends = [0, 0, 0]
    # 0-based j as in the native code; zarr[j + 1] in Julia
    j = 0
    z = zstart; zend = max(zc_HeIII, zend_rec); Dz = (z - zend) / n_low_res
    if z > zc_HeIII && z > zend
        while z >= zend + Dz
            zarr[j + 1] = z
            Xe = 1.0 + 2.0 * fHe
            TM = rf_Te_QS(θ, z, Xe, Hfun(z), c)
            Xe_Saha = rf_SahaBoltz_HeIII(rf_NH(c, θ, z), fHe, TM, c)
            if abs(Xe / Xe_Saha - 1.0) >= Saha_eps_I
                zend = z; break
            end
            z -= Dz; j += 1
        end
    end
    ends[1] = j + 1
    z = zend; zend = max(zc_HeIIIS, zend_rec); Dz = (z - zend) / n_high_res
    if z > zc_HeIIIS && z > zend
        while z >= zend + Dz
            zarr[j + 1] = z
            Xe = rf_SahaBoltz_HeIII(rf_NH(c, θ, z), fHe, rf_TCMB(θ, z), c)
            TM = rf_Te_QS(θ, z, Xe, Hfun(z), c)
            if abs(Xe / (1.0 + fHe) - 1.0) <= Saha_eps_II
                zend = z; break
            end
            z -= Dz; j += 1
        end
    end
    ends[2] = j + 1
    z = zend; zend = max(zc_HeII, zend_rec); Dz = (z - zend) / n_low_res
    if z > zc_HeII && z > zend
        while z >= zend + Dz
            zarr[j + 1] = z
            Xe = 1.0 + fHe
            TM = rf_Te_QS(θ, z, Xe, Hfun(z), c)
            Xe_Saha = rf_SahaBoltz_HeII(rf_NH(c, θ, z), fHe, TM, c)
            if abs(Xe / Xe_Saha - 1.0) >= Saha_eps_I
                zend = z; break
            end
            z -= Dz; j += 1
        end
    end
    z = zend
    zarr[j + 1] = z
    ends[3] = j + 1
    jode = j
    if z > zend_rec
        Dzl = (z - zend_rec) / (npts - j - 1)
        i = j + 1
        while i < npts && zarr[i] > 800.0          # native: i < nzpts && zarr[i-1] > 800 (0-based)
            zarr[i + 1] = zarr[i] - Dzl; i += 1
        end
        zfac = exp(log(zarr[i] / max(zend_rec, 1.0e-2)) / (npts - i))
        while i < npts && zarr[i] > zend_rec
            zarr[i + 1] = zarr[i] / zfac; i += 1
        end
        zarr[npts] = max(0.0, zend_rec)
    end
    return RecfastGrid(zarr, (ends[1], ends[2], ends[3]), jode + 1)
end

"""
    recfast_history(θ, c, Hfun, grid, solve_ode) -> NamedTuple(z, Xe_H, Xe_He, Xe, dXe, dXe_H, TM)

Native `Xe_frac` on a given [`RecfastGrid`](@ref): the three Saha segments are closed forms (`Xe_H = 1`, `Xe_He = fHe`, `dXe = dXe_H = 0`, `TM = Te_QS`), the ODE part integrates
`(xHep, xp, TM)` from the initial node `grid.jode` to the end with the callback `solve_ode(rhs!, u0, zs, znodes, p) -> 3 x length(znodes)` (states at `znodes`, descending,
`znodes[1] = zs`), with `rhs! = `[`recfast_rhs3!`](@ref) (a plain function: no captured state, as Mooncake requires) and `p = (θ, c, Hfun)`. The callback may return `(states, derivatives)`; then `dXe = d1 + d2`, `dXe_H = d2` use
the solver's derivative (native `Sz.dy` is the native solver's own derivative of its solution), otherwise the RHS at each node (ill-conditioned where the solution sits on the stiff equilibrium branch).
Native order: z descending. Columns are IONIZED populations, not neutral fractions.
"""
function recfast_history(θ, c::RecfastConstants, Hfun, grid::RecfastGrid, solve_ode)
    zarr = grid.z; n = length(zarr); j = grid.jode
    T = promote_type(eltype(θ), Float64)
    Xe_H = ones(T, n); Xe_He = fill(zero(T), n); Xe = zeros(T, n); dXe = zeros(T, n); dXH = zeros(T, n); TM = zeros(T, n)
    fHe = rf_fHe(c, θ[3])
    e1, e2, _ = grid.seg_end
    for i in 1:j
        Xe_He[i] = fHe
        if i <= e1 - 1
            Xe[i] = 1.0 + 2.0 * fHe
        elseif i <= e2 - 1
            Xe[i] = rf_SahaBoltz_HeIII(rf_NH(c, θ, zarr[i]), fHe, rf_TCMB(θ, zarr[i]), c)
        else
            Xe[i] = 1.0 + fHe
        end
        TM[i] = rf_Te_QS(θ, zarr[i], Xe[i], Hfun(zarr[i]), c)
    end
    # initial node of the ODE (native: segment-3 closed form at z = zend)
    Xe_He[j] = fHe; Xe[j] = 1.0 + fHe; TM[j] = rf_Te_QS(θ, zarr[j], Xe[j], Hfun(zarr[j]), c)
    if n > j
        u0 = [Xe_He[j], Xe_H[j], TM[j]]
        res = solve_ode(recfast_rhs3!, u0, zarr[j], zarr[j:n], (θ, c, Hfun))
        sol, dsol = res isa Tuple ? res : (res, nothing)
        for k in 2:(n - j + 1)
            i = j + k - 1
            Xe_He[i] = sol[1, k]; Xe_H[i] = sol[2, k]; TM[i] = sol[3, k]
            Xe[i] = Xe_H[i] + Xe_He[i]
            if dsol === nothing        # RHS at the node: ill-conditioned on the stiff equilibrium branch (see docs/CHUNK4C_RESULTS.md)
                f = recfast_rhs(θ, zarr[i], (Xe_He[i], Xe_H[i], TM[i], 0.0), Hfun(zarr[i]), c)
                dXe[i] = f[1] + f[2]; dXH[i] = f[2]
            else                       # derivative of the solver's dense output (native `Sz.dy` is the solver's own derivative, not the RHS)
                dXe[i] = dsol[1, k] + dsol[2, k]; dXH[i] = dsol[2, k]
            end
        end
    end
    return (z = zarr, Xe_H = Xe_H, Xe_He = Xe_He, Xe = Xe, dXe = dXe, dXe_H = dXH, TM = TM)
end
