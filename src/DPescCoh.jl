# Pure-Julia port of the explicit coherent-scattering DPesc integral of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), used by the native
# H-I absorber whenever DP_interpol_S/T leave the pre-tabulated DP table:
#   Development/Recombination/ODEdef_Xi.HeI.cpp:292-323 (DPesc_coh), Modules/DPesc_HI_abs_tabulation.cpp:55-70 (call_DP_Singlet/Triplet),
#   Development/Recombination/Pesc.HI.cpp:37-121 (Inner_Int_appr, DPesc_appr_I_sym), Development/Integration/Patterson.cpp (Patterson nested rules, driver),
#   Development/Line_profiles/Voigtprofiles.{h,cpp} (Voigtprofile_Dawson), Development/Simple_routines/routines.cpp (Dawson_Int, erf_JC, DC_sumprod),
#   Development/Hydrogenic/Photoionization_cross_section.cpp (H 1s cross section shape).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39) and
# Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK3E_RESULTS.md NOTICE).
#
# Deviations from the native code (all documented in docs/CHUNK3E_RESULTS.md):
#  * Patterson nodes/weights are GENERATED here (BigFloat Stieltjes extension + interpolatory weights), not copied from the native 16-digit tables.
#  * The H 1s cross-section SHAPE sigma(nu)/sigma(nu_thr) uses the closed-form 1s Gaunt factor (equal to the native Storey-Hummer evaluation to ~2e-15).
#  * The `Patterson order` reached by each of the seven sub-intervals is a discrete, primal-valued branch (see `dpesc_orders`).

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Patterson nested quadrature rules, 1, 3, 7, 15, 31, 63, 127, 255 points on [-1, 1]
# ---------------------------------------------------------------------------------------------------------------------------------------------------
function _bigpolymul_lin(P::Vector{BigFloat}, root2::BigFloat)   # P(t) * (t - root2), ascending coefficients
    out = zeros(BigFloat, length(P) + 1)
    for (i, c) in enumerate(P)
        out[i + 1] += c
        out[i] -= c * root2
    end
    return out
end

function _bighorner(c::Vector{BigFloat}, x::BigFloat)
    r = c[end]
    for i in (length(c) - 1):-1:1
        r = r * x + c[i]
    end
    return r
end

# monic Stieltjes polynomial Q(t = x^2) of degree m with int x P(x^2) Q(x^2) x^k dx = 0 for odd k < n+... (P from the n current nodes)
function _stieltjes_even(pos::Vector{BigFloat})
    m = length(pos)                         # n = 2m - 1 current nodes, new polynomial has degree m in t
    P = BigFloat[1]
    for i in 2:m
        P = _bigpolymul_lin(P, pos[i]^2)
    end
    M = zeros(BigFloat, m, m + 1)
    for l in 0:(m - 1), j in 0:m
        s = zero(BigFloat)
        for (i, Pi) in enumerate(P)
            s += Pi * 2 / (2 * ((i - 1) + l + 1 + j) + 1)
        end
        M[l + 1, j + 1] = s
    end
    A = M[:, 1:m]
    b = -M[:, m + 1]
    for k in 1:m                            # Gaussian elimination with partial pivoting
        piv = argmax([abs(A[r, k]) for r in k:m]) + k - 1
        if piv != k
            A[[k, piv], :] = A[[piv, k], :]
            b[[k, piv]] = b[[piv, k]]
        end
        for r in (k + 1):m
            f = A[r, k] / A[k, k]
            for c in k:m
                A[r, c] -= f * A[k, c]
            end
            b[r] -= f * b[k]
        end
    end
    c = zeros(BigFloat, m)
    for k in m:-1:1
        s = b[k]
        for j in (k + 1):m
            s -= A[k, j] * c[j]
        end
        c[k] = s / A[k, k]
    end
    return vcat(c, BigFloat(1))
end

function _bigroot(Q::Vector{BigFloat}, lo::BigFloat, hi::BigFloat)
    f(x) = _bighorner(Q, x * x)
    flo = f(lo)
    for _ in 1:80
        mid = (lo + hi) / 2
        fm = f(mid)
        if sign(fm) == sign(flo)
            lo = mid; flo = fm
        else
            hi = mid
        end
    end
    x = (lo + hi) / 2
    dQ = [Q[i + 1] * i for i in 1:(length(Q) - 1)]
    for _ in 1:12
        x -= f(x) / (2 * x * _bighorner(dQ, x * x))
    end
    return x
end

function _interpolatory_weights(pos::Vector{BigFloat})
    nodes = vcat(-reverse(pos[2:end]), pos)         # all N nodes ascending
    N = length(nodes)
    R = BigFloat[1]
    for x in nodes                                  # R(x) = prod (x - x_j), ascending coefficients
        Rn = zeros(BigFloat, length(R) + 1)
        for (i, c) in enumerate(R)
            Rn[i + 1] += c
            Rn[i] -= c * x
        end
        R = Rn
    end
    w = BigFloat[]
    for xi in pos
        S = zeros(BigFloat, N)                      # R / (x - xi) by synthetic division (descending)
        carry = zero(BigFloat)
        for k in N:-1:1
            carry = R[k + 1] + carry * xi
            S[k] = carry
        end
        Sx = _bighorner(S, xi)
        integral = zero(BigFloat)
        for k in 1:2:N                              # even powers x^(k-1): 2/k
            integral += S[k] * 2 / k
        end
        push!(w, integral / Sx)
    end
    return w
end

function _generate_patterson()
    xs = Vector{Vector{Float64}}()
    ws = Vector{Vector{Float64}}()
    setprecision(BigFloat, 2048) do
        pos = BigFloat[0]
        push!(ws, Float64.(_interpolatory_weights(pos)))
        push!(xs, Float64[])
        for _ in 1:7
            Q = _stieltjes_even(pos)
            bounds = vcat(pos, BigFloat(1))
            newx = [_bigroot(Q, bounds[i], bounds[i + 1]) for i in 1:length(pos)]
            merged = BigFloat[]
            for i in eachindex(pos)
                push!(merged, pos[i]); push!(merged, newx[i])
            end
            pos = merged
            push!(xs, Float64.(newx))
            push!(ws, Float64.(_interpolatory_weights(pos)))
        end
    end
    return xs, ws
end

const PATTERSON_NPOINTS = (1, 3, 7, 15, 31, 63, 127, 255)
const PATTERSON_X, PATTERSON_W = _generate_patterson()   # PATTERSON_X[l]: new positive nodes of level l; PATTERSON_W[l]: weights of [0, interleaved positive nodes]

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# special functions (rational approximations of routines.cpp, coefficients transcribed mechanically)
# ---------------------------------------------------------------------------------------------------------------------------------------------------
const _DAWSON_A = (5.574810327568612837e-01, 1.480548760676749464e-01, 2.473412335012745505e-02, 2.885568838534784757e-03, 2.460310593071530085e-04, 1.551300911064284532e-05, 7.087756590256507931e-07, 2.179833744510066470e-08, 3.580284771465818339e-10,)
const _DAWSON_B = (1.000000000000000001e+00, -1.091856339098055081e-01, 4.306752089643668630e-02, -1.497994816990499693e-03, 3.339251349042084215e-04, -1.857560947421670758e-06, 6.710614181693725933e-07, 4.960681954907270745e-09, 2.335008164471329286e-10,)
const _DAWSON_A1 = (6.262545365722985257e-01, 1.885254515728064716e-01, 3.762507762042703303e-02, 4.987142830383369467e-03, 6.094443149907482678e-04, 3.731642154042092311e-05, 4.239542422752186968e-06, 5.982569991142279706e-08, 1.174473438334967860e-08,)
const _DAWSON_B1 = (1.000651029373811169e+00, -4.210501168492219455e-02, 3.981798325257397695e-02, 1.025635538449943440e-03, 4.099305572015143676e-04, 1.596933865686699852e-05, 2.244388174606844525e-06, 3.090088866526039173e-08, 5.894053468010250910e-09, -1.103344321928670385e-13,)
const _DAWSON_A2 = (-4.779664691344904517e+00, 3.275776248263644817e+00, -8.957008446633805629e-01, 1.364310927775567374e-01, -1.254364750946865384e-02, 7.203420550034988907e-04, -2.404381247095368252e-05, 3.964910408777002362e-07,)
const _DAWSON_B2 = (-1.813546033357404877e+00, 1.455784804742178438e+00, -4.178728305708361334e-01, 6.532985079094380749e-02, -6.100423692880221490e-03, 3.543095925141692234e-04, -1.192279451351772837e-05, 1.982456062768223261e-07, -2.989826117053882457e-16,)
const _DAWSON_A3 = (-7.277595909795730480e+01, 2.030133303927516449e+03, -2.623905050877405481e+04, 1.572452946972041743e+05, -3.959782020608236881e+05, 3.242235333091555331e+05, -3.695969139799184010e+04,)
const _DAWSON_B3 = (4.999999999999999987e-01, -3.613797954897864312e+01, 9.972476621892572722e+02, -1.263834541306106681e+04, 7.275923849721134328e+04, -1.668382015827219241e+05, 1.031530728542497468e+05,)
const _ERF_A = (4.891837297874833514e-01, 1.110175596090841322e-01, 1.526977188817169289e-02, 1.388143322498740953e-03, 8.446154421878829637e-05, 3.239915935631695227e-06, 6.200069065781009292e-08,)
const _ERF_B = (1.128379167095512575e+00, 1.758583405424390318e-01, 5.411290829613803886e-02, 3.805760982281134630e-03, 4.314620532020313106e-04, 1.423201530735308681e-05, 6.188698151904667110e-07, 2.293915472550064153e-09,)
const _ERF_A1 = (3.040844316651787557e+01, 3.358927547920980388e+02, 1.703170048554673596e+03, 4.133195606956137156e+03, 4.555611776312694034e+03, 1.935778559646575488e+03, 2.051089518116697193e+02,)
const _ERF_B1 = (5.641895835477562749e-01, 1.687403209467950089e+01, 1.813522721872712655e+02, 8.779664433851135986e+02, 1.965115658619443782e+03, 1.865558781108286245e+03, 5.828865035607128761e+02, 2.558559157228883880e+01,)


@inline function _horner(c::NTuple{N,Float64}, y) where {N}
    r = c[N] * one(y)
    for i in (N - 1):-1:1
        r = r * y + c[i]
    end
    return r
end

"""`Dawson_Int` (routines.cpp:100-): exp(-x^2) int_0^x exp(y^2) dy (rational approximations, accuracy ~1e-15)."""
function dawson_jc(x)
    xa = abs(x)
    if xa < 2.5
        y = x * x
        return x * _horner(_DAWSON_B, y) / (_horner(_DAWSON_A, y) * y + 1)
    elseif xa < 4.0
        y = x * x
        return x * _horner(_DAWSON_B1, y) / (_horner(_DAWSON_A1, y) * y + 1)
    elseif xa < 5.5
        y = x * x
        return x * _horner(_DAWSON_B2, y) / (_horner(_DAWSON_A2, y) * y + 1)
    else
        y = 1 / (x * x)
        return _horner(_DAWSON_B3, y) / ((_horner(_DAWSON_A3, y) * y + 1) * x)
    end
end

"""`erf_JC` (routines.cpp:167-)."""
function erf_jc(x)
    if abs(x) < 2.0
        y = x * x
        return x * _horner(_ERF_B, y) / (_horner(_ERF_A, y) * y + 1)
    else
        y = 1 / (x * x)
        f = 1 - exp(-x * x) * _horner(_ERF_B1, y) / ((_horner(_ERF_A1, y) * y + 1) * abs(x))
        return x < 0 ? -f : f
    end
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Voigtprofile_Dawson (Voigtprofiles.cpp:100-300) for a line with atomic data (nu21, Gamma, AM)
# ---------------------------------------------------------------------------------------------------------------------------------------------------
"""Atomic data of a native `Voigtprofile_Dawson` object: line frequency `nu21` [Hz], damping `Gamma` [1/s], mass number `AM`."""
struct VoigtLine
    nu21::Float64
    Gamma::Float64
    AM::Float64
end

const DPESC_X_WING = 30.0
const DPESC_SQRT_PI = sqrt(3.1415926535897931)
# Voigtprofile_Base: k_mHc2 = const_kB/const_mH_gr/const_cl/const_cl, with const_kB = 1.3806504e-23*1e7, const_mH_gr = 1.0078250321*1.660538782e-24
const DPESC_K_MHC2 = (1.3806504e-23 * 1.0e+7) / (1.0078250321 * 1.660538782e-24) / 2.99792458e+10 / 2.99792458e+10

voigt_DnuT(l::VoigtLine, Tm) = l.nu21 * sqrt(2.0 * DPESC_K_MHC2 * Tm / l.AM)
voigt_a(l::VoigtLine, Tm) = l.Gamma / (4.0 * 3.1415926535897931) / voigt_DnuT(l, Tm)
voigt_x2nu(l::VoigtLine, xD, Tm) = xD * voigt_DnuT(l, Tm) + l.nu21

function voigt_phi_wings(x, a)
    x2 = x * x
    d = a * a / x2
    K0 = 1.0 + (1.5 + (15.0 / 4.0 + 105.0 / 8.0 / x2) / x2) / x2
    K1 = 1.0 + (5.0 + 105.0 / 4.0 / x2) / x2
    return a / (x2 * 3.1415926535897931) * (K0 - d * K1)
end

function _voigt_Hn(x)
    x2 = x * x
    F = dawson_jc(x)
    s = DPESC_SQRT_PI
    H0 = exp(-x2)
    H1 = 2.0 / s * (2.0 * x * F - 1.0)
    H2 = H0 * (1.0 - 2.0 * x2)
    H3 = 4.0 / 3.0 / s * (x2 - 1.0 + F * x * (3.0 - 2.0 * x2))
    H4 = H0 / 6.0 * (3.0 + 4.0 * x2 * (x2 - 3.0))
    H5 = 2.0 / 15.0 / s * (-4.0 + (9.0 - 2.0 * x2) * x2 + F * x * (15.0 - 4.0 * x2 * (5.0 - x2)))
    H6 = H0 / 90.0 * (15.0 - 2.0 * (45.0 - 2.0 * (15.0 - 2.0 * x2) * x2) * x2)
    H7 = -2.0 * (24.0 - (87.0 - (40.0 - 4.0 * x2) * x2) * x2 - x * (105.0 - (210.0 - (84.0 - 8.0 * x2) * x2) * x2) * F) / 315.0 / s
    H8 = H0 * (105.0 - (840.0 - (840.0 - (224.0 - 16.0 * x2) * x2) * x2) * x2) / 2520.0
    return (H0, H1, H2, H3, H4, H5, H6, H7, H8)
end

"""`Voigtprofile_Dawson::phi(x, a)`: Voigt profile from the Mihalas expansion in `a`, wing expansion for |x| >= 30."""
function voigt_phi(x, a)
    abs(x) >= DPESC_X_WING && return voigt_phi_wings(x, a)
    H = _voigt_Hn(x)
    r = H[1] + (H[2] + (H[3] + (H[4] + (H[5] + (H[6] + (H[7] + (H[8] + H[9] * a) * a) * a) * a) * a) * a) * a) * a
    return r / DPESC_SQRT_PI
end

function voigt_xi_wings(x, a)
    x2 = x * x
    a2 = a * a
    return -a / 3.1415926535897931 / x * (1.0 + ((0.5 - a2 / 3.0) + ((0.75 - a2) + (1.875 - 15.0 * a2 / 4.0) / x2) / x2) / x2)
end

function _voigt_IHn(t1, t2)
    t12 = t1 * t1; t22 = t2 * t2
    F1 = dawson_jc(t1); F2 = dawson_jc(t2)
    Exp1 = exp(-t12); Exp2 = exp(-t22)
    E1 = erf_jc(t1); E2 = erf_jc(t2)
    s = DPESC_SQRT_PI
    IH0 = (E2 - E1) * s * 0.5
    IH1 = 2.0 / s * (F1 - F2)
    IH2 = t2 * Exp2 - t1 * Exp1
    IH3 = -2.0 * (t2 - t1 + (1.0 - 2.0 * t22) * F2 - (1.0 - 2.0 * t12) * F1) / 3.0 / s
    IH4 = ((3.0 - 2.0 * t22) * t2 * Exp2 - (3.0 - 2.0 * t12) * t1 * Exp1) / 6.0
    IH5 = ((5.0 - 2.0 * t12) * t1 - (5.0 - 2.0 * t22) * t2 + (3.0 + 4.0 * t12 * (-3.0 + t12)) * F1 - (3.0 + 4.0 * t22 * (-3.0 + t22)) * F2) / 15.0 / s
    IH6 = (t2 * (15.0 + 4 * t22 * (-5.0 + t22)) * Exp2 - t1 * (15.0 + 4 * t12 * (-5.0 + t12)) * Exp1) / 90.0
    IH7 = ((33.0 - (28.0 - 4.0 * t12) * t12) * t1 - (33.0 - (28.0 - 4.0 * t22) * t22) * t2
           + (15.0 - (90.0 - (60.0 - 8.0 * t12) * t12) * t12) * F1 - (15.0 - (90.0 - (60.0 - 8.0 * t22) * t22) * t22) * F2) / 315.0 / s
    IH8 = (t2 * (105.0 - (210.0 - (84.0 - 8.0 * t22) * t22) * t22) * Exp2 - t1 * (105.0 - (210.0 - (84.0 - 8.0 * t12) * t12) * t12) * Exp1) / 2520.0
    return (IH0, IH1, IH2, IH3, IH4, IH5, IH6, IH7, IH8)
end

function voigt_xi_dawson(xmin, xmax, a)
    IH = _voigt_IHn(xmin, xmax)
    r = IH[1] + (IH[2] + (IH[3] + (IH[4] + (IH[5] + (IH[6] + (IH[7] + (IH[8] + IH[9] * a) * a) * a) * a) * a) * a) * a) * a
    return r / DPESC_SQRT_PI
end

"""`Voigtprofile_Dawson::xi_Int(xmin, xmax, a)`: int_xmin^xmax phi(y) dy (wing/Dawson split at |x| = 30, mirrored case analysis of the native code)."""
function voigt_xi_int(xmin, xmax, a)
    xmin == xmax && return zero(a)
    fac = 1.0
    if xmin > xmax
        xmin, xmax = xmax, xmin
        fac = -1.0
    end
    r = zero(a)
    xw = DPESC_X_WING
    if xmin <= 0 && xmax <= 0
        ab = min(xmin, -xw); bb = min(xmax, -xw)
        ab < bb && (r += voigt_xi_wings(bb, a) - voigt_xi_wings(ab, a))
        ab = max(xmin, -xw); bb = max(xmax, -xw)
        ab < bb && (r += voigt_xi_dawson(ab, bb, a))
    elseif xmin > 0 && xmax > 0
        ab = min(-xmax, -xw); bb = min(-xmin, -xw)
        ab < bb && (r += voigt_xi_wings(bb, a) - voigt_xi_wings(ab, a))
        ab = max(-xmax, -xw); bb = max(-xmin, -xw)
        ab < bb && (r += voigt_xi_dawson(ab, bb, a))
    else
        r += 1.0 - voigt_xi_int(xmax, xmax * 1.0e+6, a)
        r += -voigt_xi_int(-1.0e+6, xmin, a)
    end
    return r * fac
end

voigt_xi_int(x, a) = voigt_xi_int(-1.0e+4, x, a)

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# H 1s photoionization cross-section shape (HILyc = Photoionization_cross_section_SH(n=1,l=0,Z=1)), only the ratio sigma(nu)/sigma_nuc enters
# ---------------------------------------------------------------------------------------------------------------------------------------------------
"""H 1s ionization edge `nu_ion` [Hz]. `sig_phot_ion(nu)/sig_phot_ion_nuc()` is evaluated from the closed-form 1s Gaunt factor (Photoionization_Lyc::g_phot_ion)."""
struct HLycCross
    nu_ion::Float64
    thres::Float64
    xi_max::Float64
end
const NATIVE_HLYC = HLycCross(3288051231691769.0, 1.0e-8, 1.0e+8)

function _lyc_shape(c::HLycCross, nu)       # Kramers(nu) * g(nu), up to the nu-independent constant
    nu <= c.nu_ion * (1.0 + c.thres) && (nu = c.nu_ion * (1.0 + c.thres))
    E = nu / c.nu_ion - 1.0
    y = sqrt(E)
    return (c.nu_ion / nu)^3 * 8.0 * sqrt(3.0) * 3.1415926535897931 / (1.0 + E) * exp(-4.0 * atan(y) / y) / (1.0 - exp(-min(2.0 * 3.1415926535897931 / y, 500.0)))
end

"""`sig_phot_ion(nu) / sig_phot_ion_nuc()`; 0 below the edge or above `xi_max nu_ion` (native branches, decided on the primal value)."""
function lyc_ratio(c::HLycCross, nu)
    nu < c.nu_ion && return zero(nu)
    nu > c.xi_max * c.nu_ion && return zero(nu)
    return _lyc_shape(c, nu) / _lyc_shape(c, c.nu_ion * (1.0 + c.thres))
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Patterson driver (Patterson.cpp:304-345, 352-420, 519-552) and DC_sumprod (routines.cpp:1111)
# ---------------------------------------------------------------------------------------------------------------------------------------------------
function dc_sumprod(y::AbstractVector, w::AbstractVector, lo::Int, M::Int)
    if M <= 4
        r = y[lo] * w[lo]
        for i in 1:(M - 1)
            r += y[lo + i] * w[lo + i]
        end
        return r
    end
    N = M ÷ 2
    return dc_sumprod(y, w, lo, N) + dc_sumprod(y, w, lo + N, M - N)
end

"""
    patterson_integrate(f, a, b, epsrel, epsabs; level = 0) -> (r, level)

`Integrate_using_Patterson`: 1-point value, then the 3, 7, ... 255-point nested rules until `|r_prev - r| <= max(epsabs, |r| epsrel)`; the native code returns
the last rule (255 points) when never converged (its refinement branch is dead, `compute_integral_function_Patterson` always returns 0). With `level > 0`
the rule `level` (1-based into `PATTERSON_NPOINTS`) is applied without tests, which is the fixed-order differentiable evaluation.
Convergence decisions use primal values only.
"""
function patterson_integrate(f, a::Real, b::Real, epsrel, epsabs; level::Int = 0)
    a >= b && return (zero(f(0.5 * (a + b))), 1)
    xc = (a + b) / 2.0
    Dx = (b - a) / 2.0
    f0 = f(xc)
    fvals = [f0]
    r = 2.0 * f0 * Dx
    level == 1 && return (r, 1)
    for it in 2:8
        xnew = PATTERSON_X[it]
        w = PATTERSON_W[it]
        nx = length(xnew)
        loc = Vector{typeof(f0)}(undef, length(fvals) + nx)
        for k in 1:nx
            loc[2k - 1] = fvals[k]
            Del = Dx * xnew[k]
            loc[2k] = f(xc + Del) + f(xc - Del)
        end
        fvals = loc
        r1 = dc_sumprod(fvals, w, 1, length(fvals)) * Dx
        Dr = r - r1
        r = r1
        level == it && return (r, it)
        level == 0 && abs(_primal(Dr)) <= max(epsabs, abs(_primal(r)) * epsrel) && return (r, it)
    end
    return (r, 8)
end

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Inner (analytic) and outer (Patterson) integral of Pesc.HI.cpp
# ---------------------------------------------------------------------------------------------------------------------------------------------------
struct DPescIntegrand{A,B,C,D,E}
    Te::A
    eta_S::B
    eta_c::C
    a::D
    xDm::E
    line::VoigtLine
    hlyc::HLycCross
end

function inner_int_appr(p::DPescIntegrand, xD)
    phi = voigt_phi(xD, p.a)
    DnuT = voigt_DnuT(p.line, p.Te)
    nu = voigt_x2nu(p.line, xD, p.Te)
    eta_cont = p.eta_c * lyc_ratio(p.hlyc, nu) * DnuT / nu / phi
    xi = p.eta_S / (p.eta_S + eta_cont)
    ym = p.xDm - voigt_xi_int(xD, p.a)
    dtau_c = (p.eta_S + eta_cont) * ym
    dtau_S = p.eta_S * ym
    return phi * (1.0 - exp(-dtau_S) - xi * (1.0 - exp(-dtau_c)))
end

(p::DPescIntegrand)(xD) = inner_int_appr(p, xD) + inner_int_appr(p, -xD)

const DPESC_EDGES = (0.0, 1.0, 2.0, 4.0, 10.0, 30.0, 100.0, 1.0e+4)
const DPESC_EPSREL = 1.0e-6

"""
    dpesc_appr_I_sym(Te, eta_S, eta_c, line, hlyc, epsabs; orders = nothing) -> (Dpij, orders)

`DPesc_appr_I_sym`: int_0^1e4 [Inner(x) + Inner(-x)] dx over the seven native sub-intervals, `epsabs` for the first one and `|r| epsrel` for the following.
`orders` is the 7-tuple of Patterson levels (1-based, rule sizes `PATTERSON_NPOINTS`); if given, those fixed rules are used (differentiable evaluation).
"""
function dpesc_appr_I_sym(Te, eta_S, eta_c, line::VoigtLine, hlyc::HLycCross, epsabs; orders = nothing)
    a = voigt_a(line, Te)
    xDm = voigt_xi_int(1.0e+4, a)
    f = DPescIntegrand(Te, eta_S, eta_c, a, xDm, line, hlyc)
    r = zero(f(0.5))
    used = Int[]
    for k in 1:7
        eabs = k == 1 ? epsabs : abs(_primal(r)) * DPESC_EPSREL
        r1, lv = patterson_integrate(f, DPESC_EDGES[k], DPESC_EDGES[k + 1], DPESC_EPSREL, eabs; level = orders === nothing ? 0 : orders[k])
        r += r1
        push!(used, lv)
    end
    return r, Tuple(used)
end

"""Native `p_ij` of the singlet escape (shared with `sobolev_p`)."""
@inline _pij(tau) = sobolev_p(tau)

"""
    DPescModel(singlet::VoigtLine, triplet::VoigtLine, hlyc::HLycCross)

Atomic/profile data of the explicit DPesc integral (`HeI_Atoms.nP_S_profile(2)`, `nP_T_profile(2)`, `HeI_Atoms.HILyc`).
"""
struct DPescModel
    singlet::VoigtLine
    triplet::VoigtLine
    hlyc::HLycCross
end
const NATIVE_DPESC = DPescModel(VoigtLine(5130495483823985.0, 1800874601.4423001, 4.0), VoigtLine(5069096363992708.0, 10216177.609656001, 4.0), NATIVE_HLYC)

"""
    dpesc_coh(model, triplet, eta_S, eta_c, Tg, Te, pd; orders = nothing) -> (Pesc - PS, orders)

Native `DPesc_coh` (ODEdef_Xi.HeI.cpp:305-323) with the profile of `call_DP_Singlet` (`triplet = false`) or `call_DP_Triplet`: returns `Pesc - p_ij(eta_S)` with
`P = p_ij(pd eta_S) + Dpij`, `Pesc = pd P / (1 - (1 - pd) P)`. Note native `call_DP_*(T, eta, tauS, pd)` maps to `eta_S = tauS`, `eta_c = eta`, `Tg = Te = T`.
"""
function dpesc_coh(model::DPescModel, triplet::Bool, eta_S, eta_c, Tg, Te, pd; orders = nothing)
    line = triplet ? model.triplet : model.singlet
    PS = _pij(eta_S)
    P_d = _pij(pd * eta_S)
    Dpij, used = dpesc_appr_I_sym(Te, pd * eta_S, eta_c, line, model.hlyc, _primal(P_d) * 1.0e-5; orders = orders)
    P = P_d + Dpij
    Pesc = pd * P / (1.0 - (1.0 - pd) * P)
    return Pesc - PS, used
end

"""Callable fallback `(Tg, eta, tauS, pd) -> Pesc - PS` of the table interpolation (native `call_DP_Singlet`/`call_DP_Triplet`); `orders` pins the Patterson levels."""
function dpesc_fallback(model::DPescModel, triplet::Bool; orders = nothing)
    return (Tg, eta, tauS, pd) -> dpesc_coh(model, triplet, tauS, eta, Tg, Tg, pd; orders = orders)[1]
end
