# Chunk 4b: explicit Cosmos accessors and GSL-compatible splines of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3):
#   Development/Cosmology/Cosmos.cpp:167-180 (init: O_L closure), :393-480 (calc_coeff_X_spline), :504 (calc_spline_log), :567-578 (H), :650-781 (Saha-Boltzmann, X1s, XHeI1s, XHeII1s, NHe*),
#   Development/Cosmology/Cosmos.h:311-445 (TCMB, Xe_Seager, Xe_b, dXe_dz, Te_Tg, Te, kappa_cool, Nb, Ne, NH, Ntot, rho_g_1_cm3), :1141 (sigT),
#   Development/Cosmology/Cosmos.cpp:95-160 (Hubble::init / H / check_limits), Development/Simple_routines/routines.cpp:347-413 (GSL natural cubic splines + 1e-14 end nudge).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39) and
# Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK4B_RESULTS.md NOTICE).
#
# SCOPE: the accessors are functions of an EXPLICIT history (the 6000 Recfast++ nodes, native column meaning: Xe_H = n_e,H/n_H and Xe_He = n_e,He/n_H are IONIZED
# populations, NOT neutral X1s / XHeI1s), explicit cosmology constants and an optional explicit Hubble table. No Recfast++ integration, no cosmology closure (Omega_L, Nb0, rho_g are inputs),
# no hidden caches or loaders. Out-of-range spline queries throw (native: the GSL error handler aborts); nothing extrapolates.

"""Natural cubic spline with the exact coefficient algorithm of `gsl_interp_cspline` (tridiagonal LDL solve, `b,c,d` per interval)."""
struct NaturalCubicSpline{T}
    x::Vector{Float64}
    y::Vector{T}
    b::Vector{T}
    c::Vector{T}
    d::Vector{T}
end

function natural_cubic_spline(x::AbstractVector{<:Real}, y::AbstractVector{T}) where {T<:Real}
    n = length(x)
    (n >= 3 && length(y) == n) || throw(DimensionMismatch("natural_cubic_spline: need length(x) == length(y) >= 3"))
    all(i -> x[i + 1] > x[i], 1:(n - 1)) || throw(ArgumentError("natural_cubic_spline: x must be strictly increasing"))
    xs = Float64.(x)
    N = n - 2
    c = zeros(T, n)
    diag = Vector{Float64}(undef, N); off = Vector{Float64}(undef, N); gam = Vector{T}(undef, N)
    for i in 1:N
        h_i = xs[i + 1] - xs[i]; h_ip1 = xs[i + 2] - xs[i + 1]
        dy_i = y[i + 1] - y[i]; dy_ip1 = y[i + 2] - y[i + 1]
        g_i = h_i != 0.0 ? 1.0 / h_i : 0.0
        g_ip1 = h_ip1 != 0.0 ? 1.0 / h_ip1 : 0.0
        off[i] = h_ip1
        diag[i] = 2.0 * (h_ip1 + h_i)
        gam[i] = 3.0 * (dy_ip1 * g_ip1 - dy_i * g_i)
    end
    alpha = Vector{Float64}(undef, N); z = Vector{T}(undef, N)
    alpha[1] = diag[1]; z[1] = gam[1]
    for i in 2:N
        t = off[i - 1] / alpha[i - 1]
        alpha[i] = diag[i] - t * off[i - 1]
        z[i] = gam[i] - t * z[i - 1]
    end
    c[N + 1] = z[N] / alpha[N]
    for i in (N - 1):-1:1
        c[i + 1] = (z[i] - off[i] * c[i + 2]) / alpha[i]
    end
    b = Vector{T}(undef, n - 1); d = Vector{T}(undef, n - 1)
    for i in 1:(n - 1)
        dx = xs[i + 1] - xs[i]
        b[i] = (y[i + 1] - y[i]) / dx - dx * (c[i + 1] + 2.0 * c[i]) / 3.0
        d[i] = (c[i + 1] - c[i]) / (3.0 * dx)
    end
    return NaturalCubicSpline{T}(xs, collect(y), b, c, d)
end

"""Thrown where native `gsl_spline_eval` would call the GSL error handler (query outside `[x[1], x[end]]`)."""
struct SplineDomainError <: Exception
    x::Float64
    xmin::Float64
    xmax::Float64
end
Base.showerror(io::IO, e::SplineDomainError) = print(io, "SplineDomainError: x = ", e.x, " outside [", e.xmin, ", ", e.xmax, "] (native: GSL error handler aborts)")

"""`gsl_spline_eval`: interval = last node `<= x` (bisection, `x == x[end]` uses the last interval); throws `SplineDomainError` outside the node range."""
function spline_eval(s::NaturalCubicSpline, x)
    xv = _primal(x)
    n = length(s.x)
    (s.x[1] <= xv <= s.x[n]) || throw(SplineDomainError(xv, s.x[1], s.x[n]))
    lo = 1; hi = n
    while hi > lo + 1
        i = (hi + lo) >> 1
        if s.x[i] > xv
            hi = i
        else
            lo = i
        end
    end
    delta = x - s.x[lo]
    return s.y[lo] + delta * (s.b[lo] + delta * (s.c[lo] + delta * s.d[lo]))
end

"""Native `calc_spline_JC` (routines.cpp:392-413): `x` within 1e-14 (relative) of an end is nudged inside by 1e-14, then `gsl_spline_eval` (no extrapolation)."""
function spline_eval_native(s::NaturalCubicSpline, x)
    xmin = s.x[1]; xmax = s.x[end]
    xv = _primal(x)
    if xv <= xmin * (1.0 + 1.0e-14) && xv >= xmin * (1.0 - 1.0e-14)
        x = xmin * (1.0 + 1.0e-14)
    end
    xv = _primal(x)
    if xv >= xmax * (1.0 - 1.0e-14) && xv <= xmax * (1.0 + 1.0e-14)
        x = xmax * (1.0 - 1.0e-14)
    end
    return spline_eval(s, x)
end

"""
Explicit Cosmos constants (harness record COS2/CONST/CFG): `Y_p`, `T0`, `fac_mHemH`, `zsRe`, `Nb0`, analytic Hubble `H0, O_k, O_L, O_m, O_rel`, `rho_g_gr`, `sigT0`
(`sigT(z) = sigT0` for the default aS = mS = 1, p = 0), `cl`, Saha constants `hPlanck, kBoltz` [J], `mElect` [kg], `EH_inf_ergs`, `me_mp`, `me_malp`, `EionHeI` [J], `pi`.
`z_saha = 3500` (hard-coded in `Cosmos::X1s/XHeI1s/XHeII1s/NHeIII`).
"""
struct CosmosConstants
    Y_p::Float64
    T0::Float64
    fac_mHemH::Float64
    zsRe::Float64
    Nb0::Float64
    H0::Float64
    O_k::Float64
    O_L::Float64
    O_m::Float64
    O_rel::Float64
    rho_g_gr::Float64
    sigT0::Float64
    cl::Float64
    hPlanck::Float64
    kBoltz::Float64
    mElect::Float64
    EH_inf_ergs::Float64
    me_mp::Float64
    me_malp::Float64
    EionHeI::Float64
    pi::Float64
    z_saha::Float64
end

"""The six native splines of `Cosmos::calc_coeff_X_spline`: `ln Xe`, `dXe`, `ln Xe_H`, `dXe_H`, `rho = TM/TCMB`, `ln max(Xe_He, 1e-20)` on ascending z."""
struct RecfastSplines{A,B,C,D,E,F}
    lnXe::A
    dXe::B
    lnXeH::C
    dXeH::D
    rho::E
    lnXeHe::F
end

"""
    recfast_splines(c, z, Xe_H, Xe_He, Xe, dXe, dXe_H, TM)

Native `calc_coeff_X_spline`: the history arrays are in NATIVE order (z DESCENDING, as returned by `recombination_history`); they are reversed to ascending z.
`Xe_H = n_e,H/n_H` and `Xe_He = n_e,He/n_H` are ionized populations. `TM` in K; `rho_i = TM_i / TCMB(z_i)`.
"""
function recfast_splines(c::CosmosConstants, z, Xe_H, Xe_He, Xe, dXe, dXe_H, TM)
    M = length(z)
    all(length(a) == M for a in (Xe_H, Xe_He, Xe, dXe, dXe_H, TM)) || throw(DimensionMismatch("recfast_splines: history arrays differ in length"))
    zg = reverse(z)
    sp(y) = natural_cubic_spline(zg, y)
    return RecfastSplines(sp(log.(reverse(Xe))), sp(reverse(dXe)), sp(log.(reverse(Xe_H))), sp(reverse(dXe_H)),
                          sp(reverse([TM[i] / (c.T0 * (1.0 + z[i])) for i in 1:M])), sp(log.(max.(reverse(Xe_He), 1.0e-20))))
end

"""Native `Hubble` object with a loaded table: natural cubic spline of `ln H` vs `ln z` (see `hubble_table`) valid for `zmin < z < zmax`."""
struct HubbleTable
    spline::NaturalCubicSpline{Float64}
    zmin::Float64
    zmax::Float64
end

"""
    hubble_table(z, Hz)

Native `Hubble::init(const double *z, const double *Hz, int nz)` for the DESCENDING-z table of the CAMB adapter (`z[1] > z[2]`; the ascending branch is also implemented).
NATIVE BUG REPRODUCED: node 0 (`z = 1e-10` if the table ends at `z = 0`) keeps `ln H = 0` because the native loop fills `lHz` from `k = 1` only. The effect (H wrong for z <~ 1)
is confirmed against the native `cosmos.H(z)` outputs in `test/fixtures/native_cosmos_accessors.txt` (records ACC); it is not repaired.
"""
function hubble_table(z::AbstractVector{<:Real}, Hz::AbstractVector{<:Real})
    nz = length(z)
    length(Hz) == nz || throw(DimensionMismatch("hubble_table: z and Hz differ in length"))
    lz = zeros(nz); lHz = zeros(nz)
    if z[1] < z[2]
        lz[1] = z[1] == 0 ? log(1.0e-10) : log(z[1])
        zmin = exp(lz[1]); zmax = z[nz]
        for k in 2:nz
            lz[k] = log(z[k]); lHz[k] = log(Hz[k])
        end
    else
        lz[1] = z[nz] == 0 ? log(1.0e-10) : log(z[nz])
        zmin = exp(lz[1]); zmax = z[1]
        for k in 2:nz
            lz[k] = log(z[nz + 1 - k]); lHz[k] = log(Hz[nz + 1 - k])
        end
    end
    return HubbleTable(natural_cubic_spline(lz, lHz), zmin, zmax)
end

"""All inputs of the accessors: constants, the six history splines and an optional loaded Hubble table (`nothing`: analytic Seager form everywhere)."""
struct CosmosAccessors{SP,H}
    c::CosmosConstants
    sp::SP
    hubble::H
end

"""Native `Cosmos::H(z)`: loaded table for `zmin < z < zmax` (spline `exp(spline(log z))`), otherwise `H0 sqrt(O_L + (1+z)^2 (O_k + (1+z)(O_m + O_rel (1+z))))`."""
function cosmos_H(a::CosmosAccessors, z)
    h = a.hubble
    if h !== nothing && _primal(z) > h.zmin && _primal(z) < h.zmax
        return exp(spline_eval_native(h.spline, log(z)))
    end
    c = a.c
    zp1 = 1.0 + z
    return c.H0 * sqrt(c.O_L + zp1 * zp1 * (c.O_k + zp1 * (c.O_m + c.O_rel * zp1)))
end

cosmos_TCMB(a::CosmosAccessors, z) = a.c.T0 * (1.0 + z)
cosmos_Ttoz(a::CosmosAccessors, T) = T / a.c.T0 - 1.0
cosmos_Nb(a::CosmosAccessors, z) = a.c.Nb0 * (1.0 + z)^3
cosmos_NH(a::CosmosAccessors, z) = (1.0 - a.c.Y_p) * cosmos_Nb(a, z)
cosmos_fHe(a::CosmosAccessors) = a.c.Y_p / (4.0 * a.c.fac_mHemH) / (1.0 - a.c.Y_p)
cosmos_sigT(a::CosmosAccessors, z) = a.c.sigT0
cosmos_rho_g_1_cm3(a::CosmosAccessors, z) = a.c.rho_g_gr * (1.0 + z)^4

"""Native `Xe_Seager`: spline `exp(ln Xe)` below `zsRe`, closed form above."""
function cosmos_Xe_Seager(a::CosmosAccessors, z)
    c = a.c
    _primal(z) >= c.zsRe && return (1.0 - c.Y_p / 2.0 * (2.0 - 1.0 / c.fac_mHemH)) / (1.0 - c.Y_p)
    return exp(spline_eval_native(a.sp.lnXe, z))
end

"""Native `Xe_b`: 1 above `zsRe`, else `(1-Yp)/(1-Yp/2 (2 - 1/fac)) exp(ln Xe spline)`."""
function cosmos_Xe_b(a::CosmosAccessors, z)
    c = a.c
    _primal(z) >= c.zsRe && return one(z)
    return (1.0 - c.Y_p) / (1.0 - c.Y_p / 2.0 * (2.0 - 1.0 / c.fac_mHemH)) * exp(spline_eval_native(a.sp.lnXe, z))
end

"""Native `dXe_dz_Seager` (spline of dXe/dz, linear values; 0 above `zsRe`)."""
cosmos_dXe_dz(a::CosmosAccessors, z) = _primal(z) >= a.c.zsRe ? zero(z) : spline_eval_native(a.sp.dXe, z)

cosmos_Ne(a::CosmosAccessors, z) = cosmos_Xe_b(a, z) * (1.0 - a.c.Y_p / 2.0 * (2.0 - 1.0 / a.c.fac_mHemH)) * cosmos_Nb(a, z)
cosmos_Ntot(a::CosmosAccessors, z) = cosmos_NH(a, z) * (1.0 + cosmos_fHe(a) + cosmos_Xe_Seager(a, z))
cosmos_kappa_cool(a::CosmosAccessors, z) = 8.0 / 3.0 * cosmos_sigT(a, z) * a.c.cl / cosmos_H(a, z) * cosmos_rho_g_1_cm3(a, z) * cosmos_Ne(a, z) / cosmos_Ntot(a, z)

"""Native `Te_Tg`: `1/(1 + 1/kappa_cool)` above `zsRe`, else the `rho` spline."""
function cosmos_Te_Tg(a::CosmosAccessors, z)
    _primal(z) >= a.c.zsRe && return 1.0 / (1.0 + 1.0 / cosmos_kappa_cool(a, z))
    return spline_eval_native(a.sp.rho, z)
end
cosmos_Te(a::CosmosAccessors, z) = cosmos_TCMB(a, z) * cosmos_Te_Tg(a, z)

_saha_c1(a::CosmosAccessors, T) = 2.0 * a.c.pi * a.c.mElect * a.c.kBoltz * T / (a.c.hPlanck * a.c.hPlanck)
_EionHI(a) = a.c.EH_inf_ergs * 1.0e-7 / (1.0 + a.c.me_mp)
_EionHeII(a) = 4.0 * a.c.EH_inf_ergs * 1.0e-7 / (1.0 + a.c.me_malp)

"""Native `SahaBoltz_HeIII(T)`: Xe^Saha during HeIII -> HeII."""
function cosmos_SahaBoltz_HeIII(a::CosmosAccessors, T)
    c1 = _saha_c1(a, T); NHm3 = cosmos_NH(a, cosmos_Ttoz(a, T)) * 1.0e+6
    g = c1^1.5 * exp(-_EionHeII(a) / (a.c.kBoltz * T)) / NHm3; fHe = cosmos_fHe(a); d = 1.0 + fHe - g
    return 0.5 * (d + sqrt(d * d + 4.0 * g * (1.0 + 2.0 * fHe)))
end

"""Native `SahaBoltz_HeII(T)`: Xe^Saha during HeII -> HeI."""
function cosmos_SahaBoltz_HeII(a::CosmosAccessors, T)
    c1 = _saha_c1(a, T); NHm3 = cosmos_NH(a, cosmos_Ttoz(a, T)) * 1.0e+6
    g = 4.0 * c1^1.5 * exp(-a.c.EionHeI / (a.c.kBoltz * T)) / NHm3; d = 1.0 - g
    return 0.5 * (d + sqrt(d * d + 4.0 * g * (1.0 + cosmos_fHe(a))))
end

function _saha_g(a::CosmosAccessors, T)
    c1 = _saha_c1(a, T); NHm3 = cosmos_NH(a, cosmos_Ttoz(a, T)) * 1.0e+6
    gHeI = 4.0 * c1^1.5 * exp(-a.c.EionHeI / (a.c.kBoltz * T)) / NHm3
    gHeII = c1^1.5 * exp(-_EionHeII(a) / (a.c.kBoltz * T)) / NHm3
    return gHeI, gHeII, cosmos_Xe_Seager(a, cosmos_Ttoz(a, T))
end

"""Native `SahaBoltz_HeII1s(T)` (Xe taken from `Xe_Seager(Ttoz(T))`)."""
function cosmos_SahaBoltz_HeII1s(a::CosmosAccessors, T)
    gHeI, gHeII, Xe = _saha_g(a, T)
    return cosmos_fHe(a) * gHeI * Xe / (gHeI * (Xe + gHeII) + Xe * Xe)
end

"""Native `SahaBoltz_HeI1s(T)`."""
function cosmos_SahaBoltz_HeI1s(a::CosmosAccessors, T)
    gHeI, gHeII, Xe = _saha_g(a, T)
    return cosmos_fHe(a) * Xe * Xe / (gHeI * (Xe + gHeII) + Xe * Xe)
end

"""Native `SahaBoltz_HI1s(T)`."""
function cosmos_SahaBoltz_HI1s(a::CosmosAccessors, T)
    c1 = _saha_c1(a, T); NHm3 = cosmos_NH(a, cosmos_Ttoz(a, T)) * 1.0e+6
    g = c1^1.5 * exp(-_EionHI(a) / (a.c.kBoltz * T)) / NHm3
    Xe = cosmos_Xe_Seager(a, cosmos_Ttoz(a, T))
    return Xe / (Xe + g)
end

"""Native `X1s(z)`: Saha above `z_saha = 3500`, else `1 - exp(ln Xe_H spline)` (`Xe_H` is the IONIZED hydrogen population)."""
cosmos_X1s(a::CosmosAccessors, z) = _primal(z) >= a.c.z_saha ? cosmos_SahaBoltz_HI1s(a, cosmos_TCMB(a, z)) : 1.0 - exp(spline_eval_native(a.sp.lnXeH, z))
cosmos_Xp(a::CosmosAccessors, z) = 1.0 - cosmos_X1s(a, z)
"""Native `XHeII1s(z)` (`exp(spline)` is the IONIZED-electron helium population `Xe_He`, used as is by the native code)."""
cosmos_XHeII1s(a::CosmosAccessors, z) = _primal(z) >= a.c.z_saha ? cosmos_SahaBoltz_HeII1s(a, cosmos_TCMB(a, z)) : exp(spline_eval_native(a.sp.lnXeHe, z))
"""Native `XHeI1s(z)`: Saha above 3500, else `fHe - exp(ln Xe_He spline)`."""
cosmos_XHeI1s(a::CosmosAccessors, z) = _primal(z) >= a.c.z_saha ? cosmos_SahaBoltz_HeI1s(a, cosmos_TCMB(a, z)) : cosmos_fHe(a) - exp(spline_eval_native(a.sp.lnXeHe, z))
cosmos_NHeI(a::CosmosAccessors, z) = cosmos_XHeI1s(a, z) * cosmos_NH(a, z)
cosmos_NHeII(a::CosmosAccessors, z) = cosmos_XHeII1s(a, z) * cosmos_NH(a, z)
"""Native `NHeIII(z)`: `(fHe - XHeI1s - XHeII1s) NH` above 3500, else `1e-20 NH`."""
cosmos_NHeIII(a::CosmosAccessors, z) = _primal(z) >= a.c.z_saha ? (cosmos_fHe(a) - cosmos_XHeI1s(a, z) - cosmos_XHeII1s(a, z)) * cosmos_NH(a, z) : 1.0e-20 * cosmos_NH(a, z)

"""Explicit 4a inputs from the accessors at `z` (native `Set_*_Levels_to_Saha` reads exactly these): `SahaInputs(fHe, Xe_Seager, Xp, NH, Te, TCMB, NHeII/NH, NHeI/NH)`."""
function saha_inputs_at(a::CosmosAccessors, z)
    return SahaInputs(cosmos_fHe(a), cosmos_Xe_Seager(a, z), cosmos_Xp(a, z), cosmos_NH(a, z), cosmos_Te(a, z), cosmos_TCMB(a, z),
                      cosmos_NHeII(a, z) / cosmos_NH(a, z), cosmos_NHeI(a, z) / cosmos_NH(a, z))
end
