# Pure-Julia port of the DEFAULT neutral-helium effective-population base RHS terms of ORIGINAL CosmoRec v3.0b
# (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), i.e. the calls of Modules/ODEdef_CosmoRec.cpp:217-250 (fcn_HeI_effective):
#   ODE_HeI_effective::evaluate_effective_Rci_Ric_terms, evaluate_effective_Rij_terms (Development/Recombination/ODE_effective.cpp:~459-535),
#   ODE_effective::evaluate_2s_two_photon_decay, evaluate_Ly_n_channel (He weight w = gw/gwp), rate lookup get_rates_HeI
#   (Rec_database/Effective_Rates.HeI/get_effective_rates.HeI.cpp:236-305).
# EXCLUDED here (later dependencies): HeI diffusion correction (Diffusion_correction_HeI_is_on, ODEdef_CosmoRec.cpp:~252-262), HeI feedback
# (HeISTfeedback = 0 in production), H-I absorption of HeI photons (ODEdef_CosmoRec.cpp:~284-322; Chunk3c), population-to-radiation coupling,
# initialization/Saha, the He on/off switch.
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748) and Chluba & Sunyaev 2006
# (A&A 446, 39) cited (see docs/CHUNK3B_RESULTS.md NOTICE).
#
# State X (length 7, per H nucleus), native HeI level index + 1: 1 = 1^1S0, 2 = 2^1S0, 3 = 2^1P1, 4 = 2^3S1, 5 = 2^3P0, 6 = 2^3P1, 7 = 2^3P2.
# The four RESOLVED states (table order m = 1..4) are X[2], X[3], X[4], X[6]; X[1], X[5], X[7] are only touched through 1^1S0 couplings (X[5], X[7] never).
# All component functions ACCUMULATE into dX like native (`+=`).

"""Physical constants of the helium rate lookup (native physical_consts.h / Definitions.h values)."""
struct HeliumConstants
    twopi::Float64
    me_gr::Float64
    kB::Float64
    h::Float64
    h_kb::Float64
end

const NATIVE_HELIUM_CONSTANTS = HeliumConstants(6.2831853071795864769, 9.1093821500000007e-28, 1.3806504000000002e-16, 6.62606896e-27, 4.7992373449498863e-11)

"""Fixed atomic data of the production helium atom `Gas_of_HeI_Atoms(2,10,10,10,-2)` read from the native objects (non-differentiated).
`res_index` are the 1-based X indices of the resolved states; `ly_*` the `k = 1..nres-1` Lyman/intercombination channels (native `Get_Trans_Data(ik,1,0,0,0)`)."""
struct HeliumAtom
    gw::Vector{Float64}
    nu_ion::Vector{Float64}
    res_index::Vector{Int}
    Dnu_2s::Float64
    A2s1s::Float64
    mu_red::Float64
    ly_index::Vector{Int}
    ly_A21::Vector{Float64}
    ly_lambda21::Vector{Float64}
    ly_nu21::Vector{Float64}
    ly_w::Vector{Float64}
    intercombination_index::Int
end

const NATIVE_HELIUM_ATOM = HeliumAtom(
    [1.0, 1.0, 3.0, 3.0, 1.0, 3.0, 5.0],
    [5945204290713440.0, 960331708156763.88, 814708806889455.5, 1152842422987567.8, 876078309783460.0, 876107926720732.25, 876110217886640.0],
    [2, 3, 4, 6], 4984872582556676.0, 51.299999999999997, 0.99986292543644184,
    [3, 4, 6], [1798900000.0, 0.0001272426, 177.58000000000001],
    [5.8433431808919833e-06, 6.2556306529970167e-06, 5.9141203179626764e-06],
    [5130495483823985.0, 4792361867725872.0, 5069096363992708.0], [3.0, 3.0, 3.0], 6,
)

"""Explicit 1-D (in log Tg) helium effective-rate table: native `res_state_Data_HeI` for the `nres` resolved states, `lgTg` ascending.
`B[i, m]`, `R[i, m, j]` are the natural-log tabulated values at `lgTg[i]` (native `BiVec`, `RijMatrix`); `A` is NOT tabulated
(native enforces detailed balance explicitly) and is built from `gw`, `nu_ion`, `mu_red` of the resolved states."""
struct HeliumRateTable
    lgTg::Vector{Float64}
    B::Matrix{Float64}
    R::Array{Float64,3}
    gw::Vector{Float64}
    nu_ion::Vector{Float64}
    mu_red::Float64
end

nres(t::HeliumRateTable) = size(t.B, 2)

"""Read one native `HeI_Rates_n*_l*_S*_J*.nS_*.dat` file: returns `(lgTg, B, R)` with `B[i]`, `R[i, j]` natural-log values (cpp:105-172)."""
function read_native_helium_rate_file(path::AbstractString)
    toks = split(read(path, String))
    N = parse(Int, toks[1]); neq = parse(Int, toks[2])
    lgTg = parse.(Float64, toks[3:(2 + N)])
    B = zeros(N); R = zeros(N, neq)
    p = 2 + N
    for i in 1:N
        B[i] = parse(Float64, toks[p + 1])
        R[i, :] = parse.(Float64, toks[(p + 3):(p + 2 + neq)])
        p += 3 + neq
    end
    return lgTg, B, R
end

"""Build a [`HeliumRateTable`](@ref) from the four native res_2 files in resolved order (explicit path; nothing is cached or loaded by default)."""
function load_helium_rate_table(dir::AbstractString; nS::Int = 30, atom::HeliumAtom = NATIVE_HELIUM_ATOM)
    names = ["HeI_Rates_n2_l0_S0_J0", "HeI_Rates_n2_l1_S0_J1", "HeI_Rates_n2_l0_S1_J1", "HeI_Rates_n2_l1_S1_J1"]
    data = [read_native_helium_rate_file(joinpath(dir, string(n, ".nS_", nS, ".dat"))) for n in names]
    lgTg = data[1][1]
    N = length(lgTg); neq = size(data[1][3], 2)
    B = zeros(N, 4); R = zeros(N, 4, neq)
    for m in 1:4
        data[m][1] == lgTg || throw(ArgumentError("Tg grids of the resolved states differ"))
        B[:, m] = data[m][2]; R[:, m, :] = data[m][3]
    end
    return HeliumRateTable(lgTg, B, R, atom.gw[atom.res_index], atom.nu_ion[atom.res_index], atom.mu_red)
end

"""Detailed-balance `A[m] = exp(log(qnl / qe))`, `qe = (2 pi me mu kB Tg / h / h)^1.5`, `qnl = gw/4 exp(h_kb nu_ion / Tg)` (cpp:271-279)."""
@inline function helium_A(gw, nu_ion, mu_red, Tg, c::HeliumConstants = NATIVE_HELIUM_CONSTANTS)
    qe = (c.twopi * c.me_gr * mu_red * c.kB * Tg / c.h / c.h)^1.5
    qln = gw / 4.0 * exp(c.h_kb * nu_ion / Tg)
    return exp(log(qln / qe))
end

"""
    helium_stencil_start(table, logTg) -> j

1-based first stencil index (native `lx = row` is `j - 1`, no HI-style shift). Throws `RateTableDomainError` where native exits
(out of table: `logTg < lgTg[1]` or `> lgTg[N]`, or `lx + 3 > N`) and ALSO where native reads one element past the end (`lx + 3 == N`, undefined behaviour natively).
"""
function helium_stencil_start(t::HeliumRateTable, logTg)
    N = length(t.lgTg)
    lT = _primal(logTg)
    Tgv = exp(Float64(_val(logTg)))
    if lT < t.lgTg[1] || lT > t.lgTg[N]
        throw(RateTableDomainError(:out_of_table, Tgv, Tgv))
    end
    j = locate(t.lgTg, logTg)
    if j + 3 > N
        throw(RateTableDomainError(:stencil_exceeds_table, Tgv, Tgv))
    end
    return j
end

"""`(A, B, R)` at `Tg` (= Te): `A` (nres), `B` (nres), `R` (nres x nres, `R[m, i] = 0` for `i <= m`). Native `get_rates_HeI`."""
function get_helium_rates(t::HeliumRateTable, Tg, c::HeliumConstants = NATIVE_HELIUM_CONSTANTS)
    logTg = log(Tg)
    j = helium_stencil_start(t, logTg)
    a = lagrange_weights(t.lgTg, j, logTg)
    n = nres(t)
    T = typeof(a[1] * logTg + Tg)
    A = zeros(T, n); B = zeros(T, n); R = zeros(T, n, n)
    _fill_herates!(A, B, R, t, Tg, c, j, a)
    return (A = A, B = B, R = R)
end
# private: the rate loop into A, B, R of element type T (R[m, i <= m] untouched; shared by get_helium_rates and the private workspace RHS)
function _fill_herates!(A, B, R, t::HeliumRateTable, Tg, c::HeliumConstants, j, a)
    T = eltype(A)
    n = nres(t)
    for m in 1:n
        A[m] = helium_A(t.gw[m], t.nu_ion[m], t.mu_red, Tg, c)
        fx = zero(T)
        for k in 1:4
            fx += a[k] * t.B[j + k - 1, m]
        end
        B[m] = exp(fx)
        for i in (m + 1):n
            fx = zero(T)
            for k in 1:4
                fx += a[k] * t.R[j + k - 1, m, i]
            end
            R[m, i] = exp(fx)
        end
    end
    return nothing
end

"""Effective continuum, `dX[idx_m] += B[m] (A[m] Xe Nc - X[idx_m])`, `Nc = NH XHeII` (native evaluate_effective_Rci_Ric_terms)."""
function helium_continuum!(dX, Xe, Nc, X, A, B, atom::HeliumAtom = NATIVE_HELIUM_ATOM)
    for m in eachindex(A)
        i = atom.res_index[m]
        dX[i] += B[m] * (A[m] * Xe * Nc - X[i])
    end
    return dX
end

"""Interlevel effective rates, resolved `i < j`: `d = R[i,j] (X_i - gw_i/gw_j exp(-h_kb (nu_ion_j - nu_ion_i)/Tg) X_j)`, `dX_i -= d`, `dX_j += d`."""
function helium_interlevel!(dX, Tg, X, R, atom::HeliumAtom = NATIVE_HELIUM_ATOM, c::HeliumConstants = NATIVE_HELIUM_CONSTANTS)
    n = length(atom.res_index)
    for i in 1:n, j in (i + 1):n
        ii, jj = atom.res_index[i], atom.res_index[j]
        xij = c.h_kb * (atom.nu_ion[jj] - atom.nu_ion[ii]) / Tg
        d = R[i, j] * (X[ii] - atom.gw[ii] / atom.gw[jj] * exp(-xij) * X[jj])
        dX[ii] -= d
        dX[jj] += d
    end
    return dX
end

"""2^1S-1^1S two-photon decay (native evaluate_2s_two_photon_decay, Lambda_ind = 1): `d = A2s1s (X2s - X1s exp(-h_kb Dnu/Tg))`, `dX[1] += d`, `dX[2] -= d`."""
function helium_two_photon!(dX, Tg, X, atom::HeliumAtom = NATIVE_HELIUM_ATOM, c::HeliumConstants = NATIVE_HELIUM_CONSTANTS)
    d = atom.A2s1s * (X[2] - X[1] * exp(-c.h_kb * atom.Dnu_2s / Tg))
    dX[1] += d
    dX[2] -= d
    return dX
end

"""Sobolev Lyman / intercombination channels k = 1..3 (native loop, `spin_forbidden` skips 2^3P1); reuses [`lyman_channel!`](@ref) with weight `gw/gwp`."""
function helium_lyman!(dX, Tg, X, NH, Hz, atom::HeliumAtom = NATIVE_HELIUM_ATOM, c::RHSConstants = NATIVE_RHS_CONSTANTS; spin_forbidden::Bool = true)
    for k in eachindex(atom.ly_index)
        idx = atom.ly_index[k]
        (!spin_forbidden && idx == atom.intercombination_index) && continue
        atom.ly_A21[k] == 0 && continue
        lyman_channel!(dX, idx, Tg, X[1], X[idx], NH, Hz, atom.ly_A21[k], atom.ly_lambda21[k], atom.ly_nu21[k], c, atom.ly_w[k])
    end
    return dX
end

"""Combined base helium RHS (native call order of fcn_HeI_effective with diffusion, feedback and H-I absorption OFF) from already-evaluated rates; accumulates into `dX` (length 7)."""
function helium_base_rhs!(dX, Tg, Xe, NH, Hz, X, fHe, A, B, R, atom::HeliumAtom = NATIVE_HELIUM_ATOM, c::HeliumConstants = NATIVE_HELIUM_CONSTANTS, rc::RHSConstants = NATIVE_RHS_CONSTANTS; spin_forbidden::Bool = true)
    helium_continuum!(dX, Xe, NH * (fHe - X[1]), X, A, B, atom)
    helium_interlevel!(dX, Tg, X, R, atom, c)
    helium_two_photon!(dX, Tg, X, atom, c)
    helium_lyman!(dX, Tg, X, NH, Hz, atom, rc; spin_forbidden)
    return dX
end

"""As [`helium_base_rhs!`](@ref), evaluating the rates from the explicit `table` at `Tg`."""
function helium_base_rhs!(dX, table::HeliumRateTable, Tg, Xe, NH, Hz, X, fHe, atom::HeliumAtom = NATIVE_HELIUM_ATOM, c::HeliumConstants = NATIVE_HELIUM_CONSTANTS, rc::RHSConstants = NATIVE_RHS_CONSTANTS; spin_forbidden::Bool = true)
    r = get_helium_rates(table, Tg, c)
    return helium_base_rhs!(dX, Tg, Xe, NH, Hz, X, fHe, r.A, r.B, r.R, atom, c, rc; spin_forbidden)
end
