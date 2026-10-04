# Pure-Julia port of the DEFAULT approximate H-I absorption of HeI 2^1P-1^1S and 2^3P1-1^1S photons of ORIGINAL CosmoRec v3.0b
# (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3): Modules/ODEdef_CosmoRec.cpp:281-321 (branch + output slots/signs),
#   Development/Recombination/ODEdef_Xi.HeI.cpp:377-475 (evaluate_HI_abs_HeI, evaluate_HI_abs_HeI_Intercombination),
#   Modules/DPesc_HI_abs_tabulation.cpp:34-50 (pd), 74-257 (DP_interpol_S*), 262-... (DP_interpol_T*), 686-745 (read_DP_Data),
#   Development/Recombination/ODEdef_Xi.HeI.cpp:216-292 (fcorr_Loaded), Development/Recombination/Sobolev.cpp (p_ij, tau_S_gw, exp_nu),
#   Rec_database/Effective_Rates.HeI/get_effective_rates.HeI.cpp:440-486 (get_Bitot_HeI), Development/Simple_routines/routines.cpp (locate_JC, GSL spline).
# Production default: CAMB/CosmoRec runmode 0 sets HI_absorption = 2, which setup_functions.cpp:411 remaps to 1 with _HI_abs_appr_flag = 1 (fcorr applied to the singlet line).
# Native fallback to the explicit DPesc_coh integral (src/DPescCoh.jl, Chunk 3e) is OPT-IN through the `fallback` keyword; without it DPTableDomainError is thrown,
# NOT ported:
#   HeI diffusion/feedback, the ODE, initialization. CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748),
#   Chluba & Sunyaev 2006 (A&A 446, 39) and Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK3C_RESULTS.md NOTICE). Table data are NOT bundled.

"""Query left the pre-tabulated DP table (native would run the expensive explicit integral, which is not ported) or the Bitot table."""
struct DPTableDomainError <: Exception
    reason::Symbol
    value::Float64
end
Base.showerror(io::IO, e::DPTableDomainError) = print(io, "DPTableDomainError(", e.reason, "): ", e.value)

"""Constants of the absorption branch: `cl`, `h_kb`, `pi`, threshold H 1s cross section `sig_c` (`HILyc.sig_phot_ion_nuc()`), `zcrit_HI`, `zcrit_HI_Int` (HeI_routines.cpp:60-61)."""
struct HIAbsConstants
    cl::Float64
    h_kb::Float64
    pi::Float64
    sig_c::Float64
    zcrit_HI::Float64
    zcrit_HI_Int::Float64
end
const NATIVE_HIABS_CONSTANTS = HIAbsConstants(29979245800.0, 4.7992373449498863e-11, 3.1415926535897931, 6.3111866125271474e-18, 3400.0, 3400.0)

"""A line `upper -> 1^1S0` with `gw = 3`, `gwp = 1` (hard-coded natively) and the 1-based X indices of lower and upper states."""
struct HIAbsLine
    A21::Float64
    lambda21::Float64
    Dnu::Float64
    gw::Float64
    gwp::Float64
    lower::Int
    upper::Int
end
const NATIVE_HIABS_SINGLET = HIAbsLine(1798900000.0, 5.8433431808919833e-06, 5130495483823985.0, 3.0, 1.0, 1, 3)
const NATIVE_HIABS_TRIPLET = HIAbsLine(177.58000000000001, 5.9141203179626764e-06, 5069096363992708.0, 3.0, 1.0, 1, 6)
"""`nP_S_profile(2).Get_A21()` used in pd (equal to the singlet `A21`)."""
const NATIVE_HIABS_PROFILE_A21 = 1798900000.0

"""One ln(T) sheet of the DP table: ln-axes `eta`, `tauS` (singlet) and `tauT` (triplet) and `DP_S[ieta, itau]`, `DP_T[ieta, itau]`."""
struct DPSheet
    eta::Vector{Float64}
    tauS::Vector{Float64}
    tauT::Vector{Float64}
    DP_S::Matrix{Float64}
    DP_T::Matrix{Float64}
end

"""Explicit DP table (native `DP_collection`): `lnT` descending, one `DPSheet` per entry. Nothing is cached or loaded by default."""
struct DPTable
    lnT::Vector{Float64}
    sheets::Vector{DPSheet}
end

"""Read a native `DP_Coll_Data.31.fac_50.neff_*.dat` (read_DP_Data, cpp:686-745)."""
function read_native_dp_table(path::AbstractString)
    tok = split(read(path, String))
    pos = Ref(1)
    function block(marker)
        v = Float64[]
        while tok[pos[]] != marker
            push!(v, parse(Float64, tok[pos[]])); pos[] += 1
        end
        pos[] += 1
        return v
    end
    function mat(ne, nt)
        M = Matrix{Float64}(undef, ne, nt)
        for e in 1:ne, t in 1:nt
            M[e, t] = parse(Float64, tok[pos[]]); pos[] += 1
        end
        return M
    end
    lnT = block("TB")
    sheets = DPSheet[]
    for _ in lnT
        eta = block("tS"); tauS = block("tT"); tauT = block("DP_S")
        DP_S = mat(length(eta), length(tauS))
        tok[pos[]] == "DP_T" || throw(ArgumentError("DP file: expected DP_T"))
        pos[] += 1
        DP_T = mat(length(eta), length(tauT))
        pos[] <= length(tok) && tok[pos[]] == "TB" && (pos[] += 1)
        push!(sheets, DPSheet(eta, tauS, tauT, DP_S, DP_T))
    end
    return DPTable(lnT, sheets)
end

"""Natural cubic spline (GSL `gsl_interp_cspline`) of the `fcorr` knots `(z, f)` in ascending z; evaluation via [`fcorr`](@ref)."""
struct FcorrSpline
    z::Vector{Float64}
    f::Vector{Float64}
    M::Vector{Float64}
end

function FcorrSpline(z::Vector{Float64}, f::Vector{Float64})
    n = length(z)
    n >= 3 || throw(ArgumentError("need >= 3 knots"))
    issorted(z) || throw(ArgumentError("z knots must be ascending"))
    h = diff(z)
    M = zeros(n)
    sub = zeros(n); dia = ones(n); sup = zeros(n); rhs = zeros(n)
    for i in 2:(n - 1)
        sub[i] = h[i - 1]; dia[i] = 2 * (h[i - 1] + h[i]); sup[i] = h[i]
        rhs[i] = 6 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1])
    end
    for i in 2:n
        w = sub[i] / dia[i - 1]
        dia[i] -= w * sup[i - 1]; rhs[i] -= w * rhs[i - 1]
    end
    M[n] = rhs[n] / dia[n]
    for i in (n - 1):-1:1
        M[i] = (rhs[i] - sup[i] * M[i + 1]) / dia[i]
    end
    return FcorrSpline(z, f, M)
end

"""Read the native `f.corr.dat` (two columns z, f; ascending or descending)."""
function read_native_fcorr(path::AbstractString)
    v = parse.(Float64, split(read(path, String)))
    z = v[1:2:end]; f = v[2:2:end]
    if z[1] > z[2]
        z = reverse(z); f = reverse(f)
    end
    return FcorrSpline(z, f)
end

# calc_spline_JC (routines.cpp:392-412): clamp-free except the 1e-14 edge nudges; branches on primal values.
function spline_eval(s::FcorrSpline, x)
    xmin = s.z[1]; xmax = s.z[end]
    xp = _primal(x)
    if xp <= xmin * (1.0 + 1.0e-14) && xp >= xmin * (1.0 - 1.0e-14)
        x = xmin * (1.0 + 1.0e-14)
    end
    xp = _primal(x)
    if xp >= xmax * (1.0 - 1.0e-14) && xp <= xmax * (1.0 + 1.0e-14)
        x = xmax * (1.0 - 1.0e-14)
    end
    xp = _primal(x)
    i = clamp(searchsortedlast(s.z, xp), 1, length(s.z) - 1)
    h = s.z[i + 1] - s.z[i]
    a = (s.z[i + 1] - x) / h
    b = (x - s.z[i]) / h
    return a * s.f[i] + b * s.f[i + 1] + ((a^3 - a) * s.M[i] + (b^3 - b) * s.M[i + 1]) * h^2 / 6
end

"""`fcorr_Loaded(Tg)`: spline in `z = Tg/2.725 - 1`, clamped to the knot range (ODEdef_Xi.HeI.cpp:281-290)."""
function fcorr(s::FcorrSpline, Tg)
    z = Tg / 2.725 - 1.0
    zp = _primal(z)
    zp > s.z[end] && (z = s.z[end])
    zp < s.z[1] && (z = s.z[1])
    return spline_eval(s, z)
end

"""`ln Bitot` series of the resolved 2^1P state (native `BitotVec`) on the HeI effective-rate grid `lgTg`, for the death probability pd."""
struct BitotSeries
    lgTg::Vector{Float64}
    lnBitot::Vector{Float64}
end

"""Read the Bitot column of a native `HeI_Rates_n2_l1_S0_J1.nS_*.dat` file."""
function read_native_helium_bitot(path::AbstractString)
    toks = split(read(path, String))
    N = parse(Int, toks[1]); neq = parse(Int, toks[2])
    lgTg = parse.(Float64, toks[3:(2 + N)])
    lnB = zeros(N); p = 2 + N
    for i in 1:N
        lnB[i] = parse(Float64, toks[p + 2]); p += 3 + neq
    end
    return BitotSeries(lgTg, lnB)
end

"""Native `get_Bitot_HeI(Tg, Bitot, ires)`: 4-point Lagrange of `ln Bitot` in ln Tg (no shift). Errors where native exits / reads out of bounds."""
function helium_Bitot(s::BitotSeries, Tg)
    logTg = log(Tg)
    N = length(s.lgTg)
    lT = _primal(logTg)
    !(s.lgTg[1] <= lT <= s.lgTg[N]) && throw(DPTableDomainError(:Bitot_out_of_table, exp(lT)))
    j = locate(s.lgTg, logTg)
    j + 3 > N && throw(DPTableDomainError(:Bitot_stencil_exceeds_table, exp(lT)))
    a = lagrange_weights(s.lgTg, j, logTg)
    fx = a[1] * s.lnBitot[j] + a[2] * s.lnBitot[j + 1] + a[3] * s.lnBitot[j + 2] + a[4] * s.lnBitot[j + 3]
    return exp(fx)
end

"""Death probability of 2^1P: `pd = 1/(1 + A21/Bitot(Tg f_t)/f_b)` (pd_2P_Singlet_HeI); the triplet one is fixed to 1 natively (pd_2P_Triplet_HeI)."""
pd_singlet(b::BitotSeries, Tg; A21 = NATIVE_HIABS_PROFILE_A21, f_t = 1.0, f_b = 1.0) = 1.0 / (1.0 + A21 / helium_Bitot(b, Tg * f_t) / f_b)

"""Sobolev escape probability with the native branches (Sobolev.cpp:73-78); this is the global `p_ij`, not `ODE_HI_effective_funcs::p_ij`."""
function sobolev_p(tau)
    t = _primal(tau)
    t <= 1.0e-10 && return 1.0 - 0.5 * tau + tau * tau / 6.0
    t >= 5.0e+2 && return 1.0 / tau
    return (1.0 - exp(-tau)) / tau
end

"""`tau_S_gw` (Sobolev.cpp:52-56)."""
tau_S_gw(A21, lambda21, gwi, Ni, gwj, Nj, Hz, pi_ = 3.1415926535897931) = A21 * lambda21^3 / (8.0 * pi_ * Hz) * abs(Ni * (Nj / Ni * gwi / gwj - 1.0))

"""`exp_nu` (Sobolev.cpp:22-23), argument capped at 700."""
exp_nu(Dnu, Tg, h_kb = 4.7992373449498863e-11) = exp(min(700.0, h_kb * Dnu / Tg))

# cubic weights of DP_interpol_*_tau/_eta (uniform-spacing form, DD3 from the first two axis entries)
@inline function uniform_weights(ax, i, D, x)
    DD3 = D^3
    D1 = x - ax[i]; D2 = x - ax[i + 1]; D3 = x - ax[i + 2]; D4 = x - ax[i + 3]
    return (D2 * D3 * D4 / (-6.0 * DD3), D1 * D3 * D4 / (2.0 * DD3), D1 * D2 * D4 / (-2.0 * DD3), D1 * D2 * D3 / (6.0 * DD3))
end

@inline function clamp_start(i, n)   # `if(i>1) i--; if(i>n-4) i=n-4` (0-based native) -> 1-based
    i0 = i - 1
    i0 > 1 && (i0 -= 1)
    i0 > n - 4 && (i0 = n - 4)
    return i0 + 1
end

# DP_interpol_*_eta + _tau on one sheet. Where native calls the explicit integral (eta outside the sheet: the whole eta-level returns call_DP_*;
# tau outside: every eta row returns call_DP_*, then the four eta weights are applied) `fb` is that precomputed value, else DPTableDomainError.
function sheet_out_of_range(sh::DPSheet, triplet::Bool, lgeta, lgtau)
    tax = triplet ? sh.tauT : sh.tauS
    neta = length(sh.eta); ntau = length(tax)
    le = _primal(lgeta); lt = _primal(lgtau)
    eta_out = !(min(sh.eta[1], sh.eta[neta]) <= le <= max(sh.eta[1], sh.eta[neta]))
    tau_out = !(min(tax[1], tax[ntau]) <= lt <= max(tax[1], tax[ntau]))
    return eta_out, tau_out
end

function sheet_dp(sh::DPSheet, triplet::Bool, lgeta, lgtau, fb = nothing)
    tax = triplet ? sh.tauT : sh.tauS
    DP = triplet ? sh.DP_T : sh.DP_S
    neta = length(sh.eta); ntau = length(tax)
    le = _primal(lgeta); lt = _primal(lgtau)
    eta_out, tau_out = sheet_out_of_range(sh, triplet, lgeta, lgtau)
    if eta_out
        fb === nothing && throw(DPTableDomainError(:eta_out_of_table, exp(le)))
        return fb
    end
    if tau_out
        fb === nothing && throw(DPTableDomainError(:tau_out_of_table, exp(lt)))
        Dlge = sh.eta[2] - sh.eta[1]
        ie = clamp_start(Int(floor((le - sh.eta[1]) / Dlge)) + 1, neta)
        we = uniform_weights(sh.eta, ie, Dlge, lgeta)
        return we[1] * fb + we[2] * fb + we[3] * fb + we[4] * fb
    end
    Dlge = sh.eta[2] - sh.eta[1]
    Dlgt = tax[2] - tax[1]
    ie = clamp_start(Int(floor((le - sh.eta[1]) / Dlge)) + 1, neta)
    it = clamp_start(Int(floor((lt - tax[1]) / Dlgt)) + 1, ntau)
    wt = uniform_weights(tax, it, Dlgt, lgtau)
    we = uniform_weights(sh.eta, ie, Dlge, lgeta)
    out = zero(wt[1] * we[1] * DP[1, 1])
    for a in 1:4
        row = wt[1] * DP[ie + a - 1, it] + wt[2] * DP[ie + a - 1, it + 1] + wt[3] * DP[ie + a - 1, it + 2] + wt[4] * DP[ie + a - 1, it + 3]
        out += we[a] * row
    end
    return out
end

"""1-based first sheet of the 4-sheet T stencil (native `locate_JC` on the DESCENDING `lnT`, then `iT--` and the upper clamp); errors where native uses the explicit integral."""
function dp_sheet_start(t::DPTable, lnT)
    nT = length(t.lnT)
    lp = _primal(lnT)
    !(min(t.lnT[1], t.lnT[nT]) <= lp <= max(t.lnT[1], t.lnT[nT])) && throw(DPTableDomainError(:T_out_of_table, exp(lp)))
    j = if lp >= t.lnT[1]
        1
    elseif lp <= t.lnT[nT]
        nT
    else
        jl = 0; ju = nT
        while ju - jl > 1
            jm = (ju + jl) >> 1
            if lp >= t.lnT[jm + 1]
                ju = jm
            else
                jl = jm
            end
        end
        jl + 1
    end
    return clamp_start(j, nT)
end

"""4 x 4 x 4 interpolated DP at table coordinates `lgT = log(T)`, `lgeta`, `lgtau = log(tauS pd)` (native DP_interpol_*_tau/_eta and the T Lagrange of DP_interpol_S/T).
`fallback` is a zero-argument thunk giving native `call_DP_*` at the sheet level (evaluated at most once); without it an off-table sheet query throws."""
function dp_lookup(t::DPTable, triplet::Bool, lgT, lgeta, lgtau; fallback = nothing)
    j = dp_sheet_start(t, lgT)
    w = lagrange4(lgT, t.lnT[j], t.lnT[j + 1], t.lnT[j + 2], t.lnT[j + 3])
    fb = nothing
    if fallback !== nothing
        need = false
        for k in 0:3
            eo, to = sheet_out_of_range(t.sheets[j + k], triplet, lgeta, lgtau)
            need |= eo | to
        end
        need && (fb = fallback())
    end
    DP = w[1] * sheet_dp(t.sheets[j], triplet, lgeta, lgtau, fb)
    DP += w[2] * sheet_dp(t.sheets[j + 1], triplet, lgeta, lgtau, fb)
    DP += w[3] * sheet_dp(t.sheets[j + 2], triplet, lgeta, lgtau, fb)
    DP += w[4] * sheet_dp(t.sheets[j + 3], triplet, lgeta, lgtau, fb)
    return DP
end

"""`true` if `lnT` lies inside the T axis of the table (native `DP_interpol_*` otherwise returns `call_DP_*(Tg, |eta|, |tauS|, pd)` directly)."""
dp_T_in_table(t::DPTable, lnT) = (lp = _primal(lnT); min(t.lnT[1], t.lnT[end]) <= lp <= max(t.lnT[1], t.lnT[end]))

"""
    dp_correction(table, triplet, T, tauS, eta; pd, fc)

Native `DP_interpol_S` (`triplet = false`) / `DP_interpol_T`: returns `Pesc - p_ij(tauS)` with `P = p_ij(pd tauS) + DP`, `Pesc = fc pd P / (1 - (1 - pd) P)`.
`T` is the already f_t-scaled temperature used as the table coordinate; `fc` is `fcorr(T)` (singlet, only if the approximation flag is on) or 1.
`fallback(Tg, |eta|, |tauS|, pd)` (native `call_DP_Singlet/Triplet`, see [`dpesc_fallback`](@ref)) replaces the `DPTableDomainError`s: T outside the table returns it directly
(with the unscaled `Tg`, no fcorr), eta/tau outside a sheet inserts it as that sheet's DP, exactly like the native nested interpolation.
"""
function dp_correction(t::DPTable, triplet::Bool, T, tauS, eta; pd, fc, fallback::F = nothing, Tg = T) where {F}
    if fallback !== nothing && !dp_T_in_table(t, log(T))
        return fallback(Tg, abs(eta), abs(tauS), pd)
    end
    sheet_fb = fallback === nothing ? nothing : () -> fallback(T, abs(eta), abs(tauS), pd)
    DP = dp_lookup(t, triplet, log(T), log(abs(eta)), log(abs(tauS) * pd); fallback = sheet_fb)
    PS = sobolev_p(tauS)
    P = sobolev_p(pd * tauS) + DP
    Pesc = fc * pd * P / (1.0 - (1.0 - pd) * P)
    return Pesc - PS
end

"""Continuum optical-depth parameter `eta_c = cl NH X1s_HI sig_c / H`."""
hi_abs_eta(NH, XH1s, Hz, c::HIAbsConstants = NATIVE_HIABS_CONSTANTS) = c.cl * NH * XH1s * c.sig_c / Hz

# shared body of evaluate_HI_abs_HeI / _Intercombination: returns dXe_HeI_HI_abs_dt
function line_channel(line::HIAbsLine, Tg, Xi, Xj, NH, Hz, eta, Dpij_fun, c::HIAbsConstants)
    efac = exp_nu(line.Dnu, Tg, c.h_kb)
    dum = Xi * efac - Xj * line.gw / line.gwp
    tauS = tau_S_gw(line.A21, line.lambda21, line.gw, Xi * NH, line.gwp, Xj * NH, Hz, c.pi)
    Dpij = Dpij_fun(tauS, eta)
    return Dpij * line.A21 / (efac - 1.0) * dum
end

"""Ly-alpha-like singlet channel (`evaluate_HI_abs_HeI` with `DP_interpol_S`): `dXe = DP A21/(efac - 1) (X2P efac - 3 X1s)`. `fcorr_on` is `_HI_abs_appr_flag == 1`."""
function hi_abs_singlet(dp::DPTable, bitot::BitotSeries, fc::FcorrSpline, Tg, X1s, X2P, NH, Hz, XH1s;
                        fcorr_on::Bool = true, f_t = 1.0, f_b = 1.0, fallback::F = nothing, line::HIAbsLine = NATIVE_HIABS_SINGLET, c::HIAbsConstants = NATIVE_HIABS_CONSTANTS) where {F}
    eta = hi_abs_eta(NH, XH1s, Hz, c)
    T = Tg * f_t
    pd = pd_singlet(bitot, Tg; f_t = f_t, f_b = f_b)
    fcv = fcorr_on ? fcorr(fc, T) : one(T)
    return line_channel(line, Tg, X2P, X1s, NH, Hz, eta, (tauS, e) -> dp_correction(dp, false, T, tauS, e; pd = pd, fc = fcv, fallback = fallback === nothing ? nothing : fallback(false), Tg = Tg), c)
end

"""Intercombination channel (`evaluate_HI_abs_HeI_Intercombination` with `DP_interpol_T`; pd = 1, no fcorr)."""
function hi_abs_triplet(dp::DPTable, Tg, X1s, X2T, NH, Hz, XH1s;
                        f_t = 1.0, fallback::F = nothing, line::HIAbsLine = NATIVE_HIABS_TRIPLET, c::HIAbsConstants = NATIVE_HIABS_CONSTANTS) where {F}
    eta = hi_abs_eta(NH, XH1s, Hz, c)
    T = Tg * f_t
    return line_channel(line, Tg, X2T, X1s, NH, Hz, eta, (tauS, e) -> dp_correction(dp, true, T, tauS, e; pd = 1.0, fc = 1.0, fallback = fallback === nothing ? nothing : fallback(true), Tg = Tg), c)
end

"""
    hi_absorption_rhs!(g, iXe, iHI1s, iHeI, dp, bitot, fc, z, Tg, NH, Hz, XH1s, X; spin_forbidden, fcorr_on, diffusion_correction, f_t, f_b)

Applies ODEdef_CosmoRec.cpp:284-321 to the derivative vector `g` (ACCUMULATES, like native). `iXe`, `iHI1s`, `iHeI` are the 1-based positions of the electron equation,
the H 1s population and the first HeI level (native `g[0]`, `g[iHI]`, `g[iHeI]`); `X` is the 7-element HeI population vector (native level order).
Active only for `z <= zcrit_HI f_t` and not while the HeI diffusion correction is on; the intercombination channel additionally needs `z <= zcrit_HI_Int f_t` and `spin_forbidden`.
Returns `(dXe_S, dXe_T)` (0 for inactive channels). Signs: HeI 1s `+d`, 2^1P `-d` (S) or 2^3P1 `-d` (T), electron `+d`, H 1s `-d`.
"""
function hi_absorption_rhs!(g, iXe::Int, iHI1s::Int, iHeI::Int, dp::DPTable, bitot::BitotSeries, fc::FcorrSpline, z, Tg, NH, Hz, XH1s, X;
                            spin_forbidden::Bool = true, fcorr_on::Bool = true, diffusion_correction::Bool = false, f_t = 1.0, f_b = 1.0, fallback::F = nothing,
                            S::HIAbsLine = NATIVE_HIABS_SINGLET, T::HIAbsLine = NATIVE_HIABS_TRIPLET, c::HIAbsConstants = NATIVE_HIABS_CONSTANTS) where {F}
    zero_ = zero(Tg * NH * Hz * XH1s * X[1])
    if !(_primal(z) <= c.zcrit_HI * f_t) || diffusion_correction
        return (zero_, zero_)
    end
    dS = hi_abs_singlet(dp, bitot, fc, Tg, X[S.lower], X[S.upper], NH, Hz, XH1s; fcorr_on = fcorr_on, f_t = f_t, f_b = f_b, fallback = fallback, line = S, c = c)
    g[iHeI + S.lower - 1] += dS
    g[iHeI + S.upper - 1] += -dS
    g[iXe] += dS
    g[iHI1s] += -dS
    dT = zero_
    if _primal(z) <= c.zcrit_HI_Int * f_t && spin_forbidden
        dT = hi_abs_triplet(dp, Tg, X[T.lower], X[T.upper], NH, Hz, XH1s; f_t = f_t, fallback = fallback, line = T, c = c)
        g[iHeI + T.lower - 1] += dT
        g[iHeI + T.upper - 1] += -dT
        g[iXe] += dT
        g[iHI1s] += -dT
    end
    return (dS, dT)
end
