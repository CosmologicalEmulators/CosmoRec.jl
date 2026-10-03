# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 7a: atomic data, two-photon / Raman profile ratios and the frequency grid of the HI radiation PDE.
# Sources (ORIGINAL CosmoRec v3.0b, read-only): Development/Hydrogenic/Oscillator_strength.cpp (A_SH), Development/Simple_routines/routines.cpp
# (log10factorial, init_xarr, polint_JC/locate_JC), Development/Line_profiles/{HI_Transition_Data,HI_matrix_elements,nsnd_2gamma_profiles,
# Raman_profiles}.cpp, PDE_Problem/Solve_PDEs_grid.cpp (init_PDE_xarr_cores_HI), PDE_Problem/Solve_PDEs.cpp:207-437 (arm_PDE_solver, init_Solve_PDE_Data).
# Production configuration only (nShells = 3, nS_2gamma = 3, nS_Raman = 2): matrix elements and profiles beyond n = 3 are NOT ported and throw.
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# physical_consts.h (NIST 2008) values used by the profile/atomic code
const PDE_ALPHA = 1.0 / 137.035999679
const PDE_CL = 2.99792458e10
const PDE_A0 = 5.2917720859e-9
const PDE_ME_MP = 5.4461702177e-4
const PDE_EH_INF_HZ = 3.289841960361e15
const PDE_RY_INF_ICM = 1.0973731568527e5
const PDE_PI = 3.1415926535897932384626433832795
const PDE_FOURPI = 4.0 * PDE_PI

"""Native `dlog10factorial` / `log10factorial` (routines.cpp:263-301): recursive `log10(n) + log10(n-1) + ...` below 20, Stirling series above."""
function pde_log10factorial(n::Int)
    if n >= 20
        m = n + 1.0
        return 0.5 * log10(2.0 * PDE_PI * m) + n * log10(m) - m * 0.43429448190325176 +
               log10(1.0 + 8.333333333333e-2 / m + 3.472222222222e-3 / m / m - 2.681327160494e-3 / m / m / m - 2.294720936214e-4 / m / m / m / m +
                     7.840392217201e-4 / m / m / m / m / m + 6.972813758366e-5 / m / m / m / m / m / m - 5.921664373537e-4 / m / m / m / m / m / m / m)
    end
    r = 0.0
    for d in 2:n
        r = log10(Float64(d)) + r
    end
    return r
end

_lC(n, l) = sqrt(1.0 * n * n - l * l) / n

function _R0_root(n::Int, np::Int)
    log10r = (np + 2.0) * log10(4.0 * np * n) - log10(4.0) +
             0.5 * (pde_log10factorial(n + np) - pde_log10factorial(n - np - 1) - pde_log10factorial(2 * np - 1)) +
             log10(1.0 * n - np) * (n - np - 2.0) - log10(1.0 * n + np) * (n + np + 2.0)
    return 10.0^(log10r / np)
end

"""Native `A_SH(Z = 1, Np = 1.0, n, l, np, lp)` (Oscillator_strength.cpp:46-165): hydrogen dipole rate n,l -> np,lp [1/s] from the Storey & Hummer (1991) recursions."""
function hydrogen_A_SH(n::Int, l::Int, np::Int, lp::Int)
    (abs(l - lp) != 1 || n <= np || l >= n || lp >= np || l < 0 || lp < 0 || n < 2 || np < 1) && return 0.0
    Rp = zeros(np + 1); Rm = zeros(np + 1)          # 0-based index lt -> lt + 1
    lt = np - 1
    f = _R0_root(n, np); f1 = f; f2 = f1^2
    Rp[lt + 1] = 1.0
    Rm[lt + 1] = f2 * _lC(n, lt + 1) * Rp[lt + 1] / ((lt + 1.0) * 2.0 * _lC(n, lt))
    for lt in (np - 2):-1:lp
        Rp[lt + 1] = ((2.0 * lt + 3.0) * _lC(n, lt + 2) * Rp[lt + 2] * f1 + _lC(np, lt + 2) * Rm[lt + 3]) / ((lt + 2.0) * 2.0 * _lC(np, lt + 1))
        Rm[lt + 1] = (_lC(n, lt + 1) * Rp[lt + 1] * f2 + (2.0 * lt + 1.0) * _lC(np, lt + 1) * Rm[lt + 2] * f1) / ((lt + 1.0) * 2.0 * _lC(n, lt))
    end
    mu_red = 1.0 / (1.0 + PDE_ME_MP)
    fac = PDE_ALPHA^4 * PDE_CL / 6.0 / PDE_A0 * mu_red * (1.0 / np / np - 1.0 / n / n)^3
    A = lp == l + 1 ? fac * (l + 1.0) * (Rm[lp + 1] * f^(lp - 1.0))^2 / (2.0 * l + 1.0) :
        lp == l - 1 ? fac * l * (Rp[lp + 1] * f^(lp + 1.0))^2 / (2.0 * l + 1.0) : 0.0
    return A * (1.0 + PDE_ME_MP) / (1.0 + PDE_ME_MP / 1.0)       # Z^4 (1 + me/mp)/(1 + me/(Np mp)) with Z = 1, Np = 1
end

"""Native `HI_Transition_Data` tables (vacuum, nmax_setup = 10): `Gamma_np[n+1] = A_tot(np)/4pi`, `A_npks[k+1][n+1]`, `A_npkd[k+1][n+1]` (0-based native indices shifted by one)."""
struct HITransitionData
    Gamma_np::Vector{Float64}
    A_npks::Matrix{Float64}
    A_npkd::Matrix{Float64}
end

function hi_transition_data(; nmax::Int = 10)
    G = zeros(nmax + 1)
    G[3] = hydrogen_A_SH(2, 1, 1, 0) / PDE_FOURPI
    for d in 3:nmax
        Atot = hydrogen_A_SH(d, 1, 1, 0)
        Atot += hydrogen_A_SH(d, 1, 2, 0)
        for k in 3:(d - 1)
            Atot += hydrogen_A_SH(d, 1, k, 0) + hydrogen_A_SH(d, 1, k, 2)
        end
        G[d + 1] = Atot / PDE_FOURPI
    end
    S = zeros(nmax + 1, nmax + 1)               # rows k = 0..nmax, the native `dum` row is reused (stale entries kept, as in the C++ loop)
    dum = zeros(nmax + 1)
    for nv in 2:nmax
        dum[nv + 1] = hydrogen_A_SH(nv, 1, 1, 0)
    end
    S[2, :] .= dum
    for kv in 2:nmax
        dum[kv + 1] = 0.0
        for nv in 1:(kv - 1)
            dum[nv + 1] = hydrogen_A_SH(kv, 0, nv, 1)
        end
        for nv in (kv + 1):nmax
            dum[nv + 1] = 3.0 * hydrogen_A_SH(nv, 1, kv, 0)
        end
        S[kv + 1, :] .= dum
    end
    D = zeros(nmax + 1, nmax + 1); dum .= 0.0
    for kv in 3:nmax
        dum[kv + 1] = 0.0
        for nv in 1:(kv - 1)
            dum[nv + 1] = hydrogen_A_SH(kv, 2, nv, 1)
        end
        for nv in (kv + 1):nmax
            dum[nv + 1] = 3.0 / 5.0 * hydrogen_A_SH(nv, 1, kv, 2)
        end
        D[kv + 1, :] .= dum
    end
    return HITransitionData(G, S, D)
end

gamma_np(t::HITransitionData, n::Int) = t.Gamma_np[n + 1]
A_npks(t::HITransitionData, k::Int, n::Int) = t.A_npks[k + 1, n + 1]
A_npkd(t::HITransitionData, k::Int, n::Int) = t.A_npkd[k + 1, n + 1]

# HI_matrix_elements.cpp (only the elements needed by the production configuration)
hi_R1snp(n::Int) = (dn = Float64(n); 16.0 * ((dn - 1.0) / (dn + 1.0))^n * sqrt(dn^7 / (dn * n - 1.0)^5))
hi_R2snp(n::Int) = (dn = Float64(n); 256.0 * sqrt(2.0 * dn^7 * (dn * n - 1.0)) * ((dn - 2.0) / (n + 2.0))^n / (dn * n - 4.0)^3)
hi_R3snp(n::Int) = (dn = Float64(n); 432.0 * sqrt(3.0 * dn^7 * (dn * n - 1.0)) * (7.0 * n * n - 27.0) * ((dn - 3.0) / (n + 3.0))^n / (dn * n - 9.0)^4)
hi_R3dnp(n::Int) = (dn = Float64(n); 864.0 * sqrt(1.2 * dn^11 * (dn * n - 1.0)) * ((dn - 3.0) / (dn + 3.0))^n / (dn * n - 9.0)^4)
"""Native `Rksnp(k, n)` for k <= 3 (k == n: `-1.5 n sqrt(n^2 - 1)`); larger k is not ported (production nS_2gamma = 3)."""
function hi_Rksnp(k::Int, n::Int)
    k == n && return -1.5 * n * sqrt(n * n - 1.0)
    k == 1 && return hi_R1snp(n)
    k == 2 && return hi_R2snp(n)
    k == 3 && return hi_R3snp(n)
    throw(ArgumentError("Rksnp(k = $k) is not ported (production nS_2gamma = 3, nS_Raman = 2)"))
end
"""Native `Rkdnp(k, n)` for k == 3 (k == n: `-1.5 n sqrt(n^2 - 4)`)."""
function hi_Rkdnp(k::Int, n::Int)
    k == n && return -1.5 * n * sqrt(n * n - 4.0)
    k == 3 && return hi_R3dnp(n)
    throw(ArgumentError("Rkdnp(k = $k) is not ported (production nS_2gamma = 3, nS_Raman = 2)"))
end

# ------------------------------------------------------------------------------------------------------------------
# non-resonant matrix-element tables (two-photon-data/*.Mnr_*.5000.dat; read by the native init routines, never recomputed here)
# ------------------------------------------------------------------------------------------------------------------
const PDE_C_SD = 9.0 / 1024.0 * PDE_ALPHA^6 * PDE_CL * PDE_RY_INF_ICM / (1.0 + PDE_ME_MP)
const PDE_NU_1SC = PDE_EH_INF_HZ / (1.0 + PDE_ME_MP)
_pde_nuij(n, np) = PDE_NU_1SC * ((1.0 * np)^-2.0 - (1.0 * n)^-2.0)
_pde_lorentzian(a, b) = a / (a * a + b * b)
_pde_G_ns(n) = PDE_C_SD * (4.0 / 3.0 * (1.0 - 1.0 / n / n))^5
_pde_G_nd(n) = PDE_C_SD / 2.5 * (4.0 / 3.0 * (1.0 - 1.0 / n / n))^5

function _read_mnr(path::AbstractString, nmin::Int)
    isfile(path) || throw(ArgumentError("two-photon table not found: $path"))
    cols = Dict(n => Float64[] for n in nmin:8)
    nrow = 0
    for line in eachline(path)
        f = split(line); isempty(f) && continue
        length(f) == 1 + 9 - nmin || throw(ArgumentError("malformed row in $path"))
        for n in nmin:8
            push!(cols[n], parse(Float64, f[2 + n - nmin]))
        end
        nrow += 1
    end
    nrow == 500 || throw(ArgumentError("$path: expected 500 rows (native nsnd_2gamma/Raman_profiles_xpts), found $nrow"))
    return cols
end

"""Two-photon (`nsnd_2gamma_profiles`) and Raman (`Raman_profiles`) profile data for the production levels: non-resonant splines (GSL cspline of the
stored `x(1-x)Mnr` / `x(1/(n^2-1)-x)Mnr` columns), resonance strengths `kappa_res` and positions `yr`, plus the vacuum transition data."""
struct HIProfileData
    td::HITransitionData
    A2s1s::Float64
    tg_ns::Dict{Int,NaturalCubicSpline{Float64}}
    tg_nd::Dict{Int,NaturalCubicSpline{Float64}}
    tg_kappa_ns::Dict{Int,Vector{Float64}}
    tg_kappa_nd::Dict{Int,Vector{Float64}}
    tg_yr::Dict{Int,Vector{Float64}}
    R_ns::Dict{Int,NaturalCubicSpline{Float64}}
    R_kappa_ns::Dict{Int,Vector{Float64}}
    R_yr::Dict{Int,Vector{Float64}}
end

"""
    load_hi_profile_data(dir; A2s1s = 8.2206)

Native `init_nsnd_2gamma_profiles` + `init_Raman_profiles` restricted to the production levels (2gamma: 2s, 3s, 3d; Raman: 2s) from the ORIGINAL tables in
`dir` (= `Development/Line_profiles/two-photon-data`; never redistributed). `A2s1s` = `const_HI_A2s_1s` (2s-1s normalisation, 8.2206 s^-1 in production).
"""
function load_hi_profile_data(dir::AbstractString; A2s1s::Float64 = 8.2206)
    td = hi_transition_data()
    tg_x = init_xarr_linear(5.0e-5, 0.5, 500)
    tns = _read_mnr(joinpath(dir, "DGamma.Mnr_ns.5000.dat"), 2); tnd = _read_mnr(joinpath(dir, "DGamma.Mnr_nd.5000.dat"), 3)
    tg_ns = Dict(n => natural_cubic_spline(tg_x, tns[n]) for n in 2:3)
    tg_nd = Dict(3 => natural_cubic_spline(tg_x, tnd[3]))
    tg_kappa_ns = Dict(3 => [n < 2 ? 0.0 : hi_R1snp(n) * hi_Rksnp(3, n) for n in 0:2])
    tg_kappa_nd = Dict(3 => [n < 2 ? 0.0 : hi_R1snp(n) * hi_Rkdnp(3, n) for n in 0:2])
    tg_yr = Dict(3 => [n < 2 ? 0.0 : -(1.0 / 3 / 3 - 1.0 / n / n) / (1.0 - 1.0 / 3 / 3) for n in 0:2])
    R_x = init_xarr_linear(5.0e-5, 1.0, 500)
    rns = _read_mnr(joinpath(dir, "Raman.Mnr_ns.5000.dat"), 2)
    _read_mnr(joinpath(dir, "Raman.Mnr_nd.5000.dat"), 3)      # native read_Mnr requires both files to parse
    ni = 2; res_up = 10
    xs = _pde_nuij(res_up, 1) / _pde_nuij(ni, 1) - 1.0 - 5.0e-5
    R_ns = Dict(2 => natural_cubic_spline(R_x .* xs, rns[2]))
    R_kappa_ns = Dict(2 => [n <= ni ? 0.0 : hi_R1snp(n) * hi_Rksnp(ni, n) for n in 0:res_up])
    R_yr = Dict(2 => [n <= ni ? 0.0 : (1.0 / ni / ni - 1.0 / n / n) / (1.0 - 1.0 / ni / ni) for n in 0:res_up])
    return HIProfileData(td, A2s1s, tg_ns, tg_nd, tg_kappa_ns, tg_kappa_nd, tg_yr, R_ns, R_kappa_ns, R_yr)
end

function _tg_Mnr(s::NaturalCubicSpline, y)
    y > 0.5 && return _tg_Mnr(s, 1.0 - y)
    return spline_eval_native(s, max(5.0e-5, y)) / y / (1.0 - y)
end

function _tg_sigma_nsd(y, ni::Int, Mnr, kappa, tgamma)
    r = Mnr * Mnr
    fre = zeros(typeof(y), ni); fim = zeros(typeof(y), ni)
    for n in 2:(ni - 1)
        dy1 = (Float64(ni)^-2.0 - Float64(n)^-2.0) / (1.0 - Float64(ni)^-2.0) + y
        dy2 = (1.0 - Float64(n)^-2.0) / (1.0 - Float64(ni)^-2.0) - y
        g = tgamma(n) / _pde_nuij(ni, 1)
        fre[n] = _pde_lorentzian(dy1, g) + _pde_lorentzian(dy2, g)
        fim[n] = _pde_lorentzian(g, dy1) - _pde_lorentzian(g, dy2)
    end
    for n in 2:(ni - 1)
        r += kappa[n + 1]^2 * (fre[n]^2 + fim[n]^2)
        for m in 2:(n - 1)
            r += 2.0 * kappa[n + 1] * kappa[m + 1] * (fre[n] * fre[m] + fim[n] * fim[m])
        end
    end
    r_int = zero(r)
    for n in 2:(ni - 1)
        r_int += kappa[n + 1] * fre[n]
    end
    r += 2.0 * Mnr * r_int
    return r * (y * (1.0 - y))^3
end

"""Native `sigma_2s_1s_2gamma(y)` = `G_ns(2) A2s1s Mnr^2 (y(1-y))^3` for 0 <= y <= 1, else 0."""
function sigma_2s1s_2gamma(p::HIProfileData, y)
    (y < 0.0 || y > 1.0) && return zero(y)
    Mnr = _tg_Mnr(p.tg_ns[2], y)
    return _pde_G_ns(2) * p.A2s1s * Mnr * Mnr * (y * (1.0 - y))^3
end

"""Native `sigma_ns_1s_2gamma(n, y)` (n = 2: the 2s-1s profile; n = 3 only beyond that)."""
function sigma_ns1s_2gamma(p::HIProfileData, n::Int, y)
    n == 2 && return sigma_2s1s_2gamma(p, y)
    (y < 0.0 || y > 1.0 || n < 2 || n > 8) && return zero(y)
    n == 3 || throw(ArgumentError("2gamma ns profile n = $n is not ported"))
    return _pde_G_ns(n) * _tg_sigma_nsd(y, n, _tg_Mnr(p.tg_ns[n], y), p.tg_kappa_ns[n], k -> gamma_np(p.td, k))
end

function sigma_nd1s_2gamma(p::HIProfileData, n::Int, y)
    (y < 0.0 || y > 1.0 || n < 3 || n > 8) && return zero(y)
    n == 3 || throw(ArgumentError("2gamma nd profile n = $n is not ported"))
    return _pde_G_nd(n) * _tg_sigma_nsd(y, n, _tg_Mnr(p.tg_nd[n], y), p.tg_kappa_nd[n], k -> gamma_np(p.td, k))
end

function _tg_lorentz_sum(p::HIProfileData, ni::Int, y, Ak)
    r = zero(y)
    for n in 2:(ni - 1)
        g = gamma_np(p.td, n)
        yr = p.tg_yr[ni][n + 1]
        r += Ak(ni, n) * A_npks(p.td, 1, n) / (PDE_FOURPI * g) *
             (_pde_lorentzian(g / _pde_nuij(ni, 1), y - yr) + _pde_lorentzian(g / _pde_nuij(ni, 1), y - (1.0 - yr)))
    end
    return r / PDE_PI
end

sigma_ns1s_2gamma_ratio(p::HIProfileData, n::Int, y) = n == 2 ? one(y) : sigma_ns1s_2gamma(p, n, y) / _tg_lorentz_sum(p, n, y, (a, b) -> A_npks(p.td, a, b))
sigma_nd1s_2gamma_ratio(p::HIProfileData, n::Int, y) = sigma_nd1s_2gamma(p, n, y) / _tg_lorentz_sum(p, n, y, (a, b) -> A_npkd(p.td, a, b))

"""Native `sigma_ns_1s_Raman_ratio(2, y)` = `sigma_2s_1s_Raman_ratio(y)`: full 2s-1s Raman cross section over the sum of Lorentzians of the np (n = 3..10) resonances."""
function sigma_2s1s_Raman_ratio(p::HIProfileData, y)
    ni = 2; nmax = 10
    sig = if y < 0.0
        zero(y)
    else
        Mnr = spline_eval_native(p.R_ns[ni], max(5.0e-5, y)) / y / (1.0 / (ni * ni - 1.0) - y)
        kappa = p.R_kappa_ns[ni]
        fre = zeros(typeof(y), nmax + 1); fim = zeros(typeof(y), nmax + 1)
        for n in (ni + 1):nmax
            dy1 = (Float64(ni)^-2.0 - Float64(n)^-2.0) / (1.0 - Float64(ni)^-2.0) - y
            dy2 = (1.0 - Float64(n)^-2.0) / (1.0 - Float64(ni)^-2.0) + y
            g = gamma_np(p.td, n) / _pde_nuij(ni, 1)
            fre[n + 1] = _pde_lorentzian(dy1, g) + _pde_lorentzian(dy2, g)
            fim[n + 1] = _pde_lorentzian(g, dy1) - _pde_lorentzian(g, dy2)
        end
        r = Mnr * Mnr
        for n in (ni + 1):nmax
            r += kappa[n + 1]^2 * (fre[n + 1]^2 + fim[n + 1]^2)
            for m in (ni + 1):(n - 1)
                r += 2.0 * kappa[n + 1] * kappa[m + 1] * (fre[n + 1] * fre[m + 1] + fim[n + 1] * fim[m + 1])
            end
        end
        r_int = zero(r)
        for n in (ni + 1):nmax
            r_int += kappa[n + 1] * fre[n + 1]
        end
        r += 2.0 * Mnr * r_int
        _pde_G_ns(ni) * r * (y * (1.0 + y))^3
    end
    L = zero(y)
    for n in (ni + 1):nmax
        g = gamma_np(p.td, n)
        L += A_npks(p.td, ni, n) * A_npks(p.td, 1, n) / (PDE_FOURPI * g) * _pde_lorentzian(g / _pde_nuij(ni, 1), y - p.R_yr[ni][n + 1])
    end
    return sig / (L / PDE_PI)
end

# ------------------------------------------------------------------------------------------------------------------
# frequency grid (Solve_PDEs_grid.cpp) and the arm_PDE_solver profile ratios
# ------------------------------------------------------------------------------------------------------------------
# native init_xarr writing npts points at 0-based offset `o` of `a` (accumulated like the C++ loop)
function _init_xarr!(a::Vector{Float64}, o::Int, x0::Float64, xm::Float64, npts::Int, method::Int)
    if method == 0
        dx = (xm - x0) / (npts - 1); a[o + 1] = x0
        for i in 1:(npts - 1)
            a[o + i + 1] = a[o + i] + dx
        end
    else
        dx = (xm / x0)^(1.0 / (npts - 1)); a[o + 1] = x0
        for i in 1:(npts - 1)
            a[o + i + 1] = a[o + i] * dx
        end
    end
    return a
end

# native init_PDE_xarr_HI (Solve_PDEs_grid.cpp:4-80) writing at 0-based offset `o`
function _init_pde_xarr_hi!(a::Vector{Float64}, o::Int, xmin, xmax, xres, xD; core::Int, dec::Int)
    xD *= xres
    core_width = 10.0; enhance_fac_up = 2.4
    Doppler_core = (min(xmax, xres + core_width * xD) - max(xmin, xres - core_width * xD)) / xD
    ncore = trunc(Int, core * Doppler_core)
    distance_treshold = 40.0
    Doppler_up = (xmax - xres) / xD - core_width
    Doppler_low = (xres - xmin) / xD - core_width
    diluted_core = core ÷ 2
    init_up = 1; init_low = 1; min_np_up = 10.0; min_np_low = 10.0
    if Doppler_up < distance_treshold
        init_up = 0; min_np_up = Doppler_up * diluted_core
    end
    if Doppler_low < distance_treshold
        init_low = 0; min_np_low = Doppler_low * diluted_core
    end
    dec_up = log10(xmax - xres) - log10(core_width * xD)
    dec_low = log10(xres - xmin) - log10(core_width * xD)
    nup = trunc(Int, max(min_np_up, enhance_fac_up * dec * dec_up))
    nlow = trunc(Int, max(min_np_low, dec * dec_low))
    Doppler_up <= 0.0 && (nup = 0)
    Doppler_low <= 0.0 && (nlow = 0)
    nused = nup + nlow + ncore
    if nup > 0
        _init_xarr!(a, o + nlow + ncore, core_width * xD, xmax - xres, nup, init_up)
        for k in 0:(nup - 1)
            a[o + nlow + ncore + k + 1] += xres
        end
    end
    _init_xarr!(a, o + nlow, max(xmin, xres - core_width * xD), min(xmax, xres + core_width * xD), ncore + 1, 0)
    if nlow > 0
        Dx = zeros(max(nused, nlow + 1))
        _init_xarr!(Dx, 0, core_width * xD, xres - xmin, nlow + 1, init_low)
        for k in 0:nlow
            a[o + k + 1] = xres - Dx[nlow - k + 1]
        end
    end
    for k in 1:(nused - 1)
        a[o + k] > a[o + k + 1] && error("init_PDE_xarr_HI: grid non-monotonic at $(k - 1)")
    end
    return nused
end

"""
    hi_pde_grid(nresmax; core = 25, dec = 50, npts_2s1s = 400, xmin = 1e-4)

Native `init_Solve_PDE_Data(nresmax)` grid (`init_PDE_xarr_cores_HI`, Solve_PDEs_grid.cpp:85-192): `x = nu/nu_Lya` from `xmin` to
`(1 - 1/(nresmax+1)^2)/0.75`, resonances `(1 - 1/n^2)/0.75` (n = 2..nresmax). Returns `(x, resonances)`.
"""
function hi_pde_grid(nresmax::Int; core::Int = 25, dec::Int = 50, npts_2s1s::Int = 400, xmin::Float64 = 1.0e-4)
    res = [(1.0 - 1.0 / n / n) / 0.75 for n in 2:nresmax]
    xmax = (1.0 - 1.0 / (nresmax + 1) / (nresmax + 1)) / 0.75
    a = zeros(100000)
    nres = length(res); np_2s = npts_2s1s; xcrit = 0.90; xD_f_xarr = 2.0e-5
    xmax_all = xmax
    xmax = max(xmax * 0.98, 0.75 * xmax + 0.25 * res[nres])
    dx = zeros(nres + 1)
    dx[1] = xmax - xmin
    nres > 1 && (dx[1] = (res[1] + res[2]) * 0.5 - xmin)
    for i in 1:(nres - 2)
        dx[i + 1] = (res[i + 2] - res[i]) * 0.5
    end
    nres > 1 && (dx[nres] = xmax - (res[nres - 1] + res[nres]) * 0.5)
    xl = xmin; istart = 0
    if xmin <= xcrit
        xcrit_2 = 0.05; np_2s_log = 0
        if xmin <= xcrit_2
            dens_log = 30
            dlog = log10(xcrit_2 / xmin)
            np_2s_log = trunc(Int, dlog * dens_log)
            _init_xarr!(a, istart, xl, xcrit_2, np_2s_log, 1)
            istart += np_2s_log - 1
            xl = xcrit_2
        end
        _init_xarr!(a, istart, xl, xcrit, np_2s - np_2s_log, 0)
        istart += np_2s - np_2s_log - 1
        xl = xcrit
        dx[1] += xmin - xl
    end
    if nres > 1
        nused = _init_pde_xarr_hi!(a, istart, xl, xl + dx[1], res[1], xD_f_xarr; core = core, dec = dec)
        istart += nused - 1
        for i in 1:(nres - 1)
            xl += dx[i]
            nused = _init_pde_xarr_hi!(a, istart, xl, xl + dx[i + 1], res[i + 1], xD_f_xarr; core = core, dec = dec)
            istart += nused - 1
        end
    else
        nused = _init_pde_xarr_hi!(a, istart, xl, xl + dx[1], res[1], xD_f_xarr; core = core, dec = dec)
        istart += nused - 1
    end
    nused = nres == 0 ? npts_2s1s : 100
    _init_xarr!(a, istart, xmax, xmax_all, nused, 0)
    istart += nused - 1
    for k in 1:(istart - 1)
        a[k] > a[k + 1] && error("init_PDE_xarr_cores_HI: grid non-monotonic at $k")
    end
    return a[1:(istart + 1)], res
end

"""Static data of the armed native PDE solver (`arm_PDE_solver`, Solve_PDEs.cpp:207-366) for the production run: grid, resonances, `index_2` /
`index_emission` (0-based native indices), the profile ratios (`[m]` = native index m-1, i.e. n = m+1) and the vacuum A matrices."""
struct HIPDESetup
    x::Vector{Float64}
    resonances::Vector{Float64}
    nresmax::Int
    n2g::Int
    nR::Int
    index_2::Int
    index_emission::Vector{Int}
    ratio_R_ns::Vector{Vector{Float64}}
    ratio_R_nd::Vector{Vector{Float64}}
    ratio_2g_ns::Vector{Vector{Float64}}
    ratio_2g_nd::Vector{Vector{Float64}}
    A_npns::Vector{Vector{Float64}}
    A_npnd::Vector{Vector{Float64}}
    prof::HIProfileData
end

"""
    hi_pde_setup(prof; nShells = 3, nS_2gamma = 3, nS_Raman = 2)

Production PDE setup of `compute_DPesc_with_diffusion_equation_effective` (Solve_PDEs.cpp:658-737): `nresmax`, two-photon/Raman maxima, grid and
`arm_PDE_solver` ratios. Only the production configuration (nresmax = 3, two-gamma up to n = 3, Raman 2s) is ported.
"""
function hi_pde_setup(prof::HIProfileData; nShells::Int = 3, nS_2gamma::Int = 3, nS_Raman::Int = 2)
    nR_in = min(nS_Raman, nS_2gamma)
    n2g = nS_2gamma < 2 ? 0 : min(8, min(nS_2gamma, nShells))
    nR = nR_in < 2 ? 0 : min(7, min(nR_in, nShells - 1))
    nresmax = min(8, nShells)
    nS_2gamma >= 3 && (nresmax = min(nS_2gamma, nresmax))
    (nresmax == 3 && n2g == 3 && nR == 2) || throw(ArgumentError("only the production PDE configuration (nresmax = 3, nS_2gamma = 3, nS_Raman = 2) is ported"))
    x, res = hi_pde_grid(nresmax)
    npts = length(x)
    nm = nresmax - 1
    Rns = [zeros(npts) for _ in 1:nm]; Rnd = [zeros(npts) for _ in 1:nm]; Gns = [zeros(npts) for _ in 1:nm]; Gnd = [zeros(npts) for _ in 1:nm]
    for k in 1:npts
        Rns[1][k] = sigma_2s1s_Raman_ratio(prof, x[k] - 1.0)
    end
    for ni in max(2, nR + 1):nresmax
        Rns[ni - 1] .= 1.0; Rnd[ni - 1] .= 1.0
    end
    for k in 1:npts
        Gns[1][k] = sigma_2s1s_2gamma(prof, x[k])
    end
    for ni in 3:n2g
        xres_loc = (1.0 - 1.0 / ni / ni) / 0.75
        for k in 1:npts
            xs = x[k] / xres_loc
            Gns[ni - 1][k] = sigma_ns1s_2gamma_ratio(prof, ni, xs)
            Gnd[ni - 1][k] = sigma_nd1s_2gamma_ratio(prof, ni, xs)
        end
    end
    for ni in max(2, n2g + 1):nresmax
        Gns[ni - 1] .= 1.0; Gnd[ni - 1] .= 1.0
    end
    index_2 = findfirst(>=(0.5), x) - 1
    iem = Int[]; idx = 0
    for k in eachindex(res)
        while idx < npts
            if x[idx + 1] >= res[k]
                push!(iem, idx); break
            end
            idx += 1
        end
    end
    isempty(iem) && push!(iem, npts)
    Ans = [zeros(nresmax) for _ in 1:nm]; And = [zeros(nresmax) for _ in 1:nm]
    for k in 0:(length(res) - 1), m in 0:(length(res) - 1)
        Ans[m + 1][k + 1] = A_npks(prof.td, m + 2, k + 2)
        And[m + 1][k + 1] = A_npkd(prof.td, m + 2, k + 2)
    end
    return HIPDESetup(x, res, nresmax, n2g, nR, index_2, iem, Rns, Rnd, Gns, Gnd, Ans, And, prof)
end
