# Chunk 6a: coefficients of the HI diffusion PDE of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3):
#   PDE_Problem/Load.Populations.HI.cpp (compute_Xi_HI_splines, calc_HI_Xe/rho/Xi), PDE_Problem/HI_pd_Rp_splines_effective.cpp (calc_HI_Rp_Rm_pd_all,
#   set_up_splines_for_HI_pd_Rp_effective with the X_Data overload, calc_HI_pd_i / Dnem_i), Rec_database/Effective_Rates.HI/get_effective_rates.HI.cpp:350-460 (get_rates_all).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39) and
# Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK6A_RESULTS.md NOTICE).
#
# SCOPE: explicit population splines from the rows of a previous pass, the full effective rates (`get_rates_all`: A, B, Bitot, all R), death probabilities `pd` and emission
# coefficients `Rp/Rm`, and the `pd`/`Dnem` coefficient splines used by the PDE. Production setting a_scale = m_scale = 1 (f_a = f_b = f_t = 1). NOT here: the PDE itself (phase 7),
# the correction integrals (phase 8), the feedback into the ODE (phase 9).

"""`ln Bitot` of the resolved HI states on the effective-rate `lgTg` grid (the `Bitot` column of the native `Rates_n*_l*.nS_*.dat`), `lnBitot[i, m]`."""
function load_hi_bitot_table(dir::AbstractString, nS_effective::Integer, states::ResolvedStates)
    parts = [open(read_native_rate_file, joinpath(dir, "Rates_n$(states.n[m])_l$(states.l[m]).nS_$(nS_effective).dat")) for m in eachindex(states.n)]
    return reduce(hcat, [p.Bitot for p in parts])          # the native file stores ln Bitot (like Bi, Rij, Ai)
end

"""
    get_rates_all(table, lnBitot, Tg, Te) -> (A, B, Bitot, R)

Native `get_rates_all` (get_effective_rates.HI.cpp:350-460): same stencil, detailed-balance switch and `eps_A_effective` clamp as `get_rates` (Chunk 2), plus `Bitot` and ALL
`R[m, i]` (`i != m`; `R[m, m] = 0`), with `R[m, i]` the rate from resolved state m to i (native `RijVec[m][i]`).
"""
function get_rates_all(t::AtomicRateTable, lnBitot::AbstractMatrix, Tg::Real, Te::Real)
    lx, ly, a, b, db = _setup(t, Tg, Te)
    n, neq = n_resolved(t), n_eq(t)
    A = [_A(t, m, lx, ly, a, b, Tg, db) for m in 1:n]
    B = [_B(t, m, lx, a) for m in 1:n]
    Bt = [exp(a[1] * lnBitot[lx, m] + a[2] * lnBitot[lx + 1, m] + a[3] * lnBitot[lx + 2, m] + a[4] * lnBitot[lx + 3, m]) for m in 1:n]
    T = typeof(a[1])
    R = zeros(T, n, neq)
    for m in 1:n, i in 1:neq
        i == m && continue
        fx = zero(T)
        for k in 1:4
            fx += a[k] * t.R[k + lx - 1, m, i]
        end
        R[m, i] = exp(fx)
    end
    return A, B, Bt, R
end

"""
    HIPopulationSplines(rows; zs, ze)

Native `compute_Xi_HI_splines(zs, ze, X_Data)`: `rows` are the stored rows of a pass (native `pass_on_the_Solution_CosmoRec`, z descending, columns `z, Xe, X(1s, 2s, 2p, 3s, 3p, 3d), rho, ...`).
As natively, the LAST row is not used (`nz = size - 1`); the range must satisfy `zs <= z_1` and `ze >= z_nz` (native exit otherwise: here an error). Natural cubic splines (4b) on ascending z
of `ln Xe`, `ln X_i` and `rho - 1`.
"""
struct HIPopulationSplines{S}
    lnXe::S
    lnX::Vector{S}
    rho1::S
end

function HIPopulationSplines(rows::AbstractMatrix; zs::Real, ze::Real)
    nz = size(rows, 1) - 1
    (zs > rows[1, 1] || ze < rows[nz, 1]) && throw(ArgumentError("compute_Xi_HI_splines: check redshift range $zs $ze $(rows[1, 1]) $(rows[nz, 1])"))
    za = [Float64(_primal(rows[nz + 1 - k, 1])) for k in 1:nz]      # redshifts are not differentiated
    col(c, f) = [f(rows[nz + 1 - k, c]) for k in 1:nz]
    lnXe = natural_cubic_spline(za, col(2, log))
    lnX = [natural_cubic_spline(za, col(2 + i, log)) for i in 1:6]
    rho1 = natural_cubic_spline(za, col(9, x -> x - 1.0))
    return HIPopulationSplines(lnXe, lnX, rho1)
end

hi_Xe(p::HIPopulationSplines, z) = exp(spline_eval_native(p.lnXe, z))
hi_rho(p::HIPopulationSplines, z) = spline_eval_native(p.rho1, z) + 1.0
"""Native `calc_HI_Xi(z, i)` with the native 0-based HI level index `i` (0 = 1s, 1 = 2s, 2 = 2p, 3 = 3s, 4 = 3p, 5 = 3d)."""
hi_Xi(p::HIPopulationSplines, z, i::Integer) = exp(spline_eval_native(p.lnX[i + 1], z))

"""Atomic data of the resolved HI levels used by the PDE coefficients (native objects; record LVL6): HI index, n, l, `A21(1,0)`, `lambda21(1,0)` [cm], `Dnu_1s` [Hz], and `A2s1s`."""
struct HIPDELevels
    index::Vector{Int}
    n::Vector{Int}
    l::Vector{Int}
    A21::Vector{Float64}
    lambda21::Vector{Float64}
    Dnu_1s::Vector{Float64}
    A2s1s::Float64
end

const NATIVE_HI_PDE_LEVELS = HIPDELevels([1, 2, 3, 4, 5], [2, 2, 3, 3, 3], [0, 1, 0, 1, 2], [0.0, 626490298.80265749, 0.0, 167252720.15153983, 0.0],
    [0.0, 1.2156844561319915e-05, 0.0, 1.025733759861368e-05, 0.0], [2466038423768827.0, 2466038423768827.0, 2922712205948239.0, 2922712205948239.0, 2922712205948239.0],
    8.2205999999999992)

"""
    hi_rp_rm_pd(z, pops, table, lnBitot, levels, Tg, Te, Ne, Xp, NH) -> (Rp_Rm, pd)

Native `calc_HI_Rp_Rm_pd_all` (f_a = f_b = f_t = 1): `Rp = Ne Xp B A + sum_{j != nr} R[j, nr] X_j(z)`, `Rp_Rm = Rp NH / Bitot`, `pd = Bitot / (A_{2s1s} + Bitot)` for 2s and
`Bitot / (A21(1,0) + Bitot)` otherwise (so `pd = 1` for 3s, 3d whose native `A21(1,0)` is 0).
"""
function hi_rp_rm_pd(z, pops::HIPopulationSplines, t::AtomicRateTable, lnBitot::AbstractMatrix, lv::HIPDELevels, Tg, Te, Ne, Xp, NH)
    A, B, Bt, R = get_rates_all(t, lnBitot, Tg, Te)
    n = length(A)
    T = promote_type(eltype(A), eltype(B), typeof(Ne), typeof(NH), typeof(hi_Xi(pops, z, 0)))
    RpRm = Vector{T}(undef, n); pd = Vector{T}(undef, n)
    for nr in 1:n
        Rp = Ne * Xp * B[nr] * A[nr]
        Rm = Bt[nr]
        for j in 1:n
            j == nr && continue
            Rp += R[j, nr] * hi_Xi(pops, z, lv.index[j])
        end
        RpRm[nr] = Rp * NH / Rm
        pd[nr] = lv.index[nr] == 1 ? Rm / (lv.A2s1s + Rm) : Rm / (lv.A21[nr] + Rm)
    end
    return RpRm, pd
end

"""
    _hi_pd_only(z, pops, table, lnBitot, levels, Tg, Te, Ne, NH) -> pd

Private: exactly the `pd` output of [`hi_rp_rm_pd`](@ref) (same `_setup` finiteness/table-domain checks, same `Bitot` interpolation and `Rm / (A + Rm)` expressions, same promoted element type), without the A/B/R rates
and the `Rp` sums whose results [`hi_pde_coefficients`](@ref) discarded. `pd` depends only on `Bitot` and the atomic decay rate. The element type of the rate arrays that `hi_rp_rm_pd` promotes over is obtained from the
types of the stencil weights and one `log_qnl_qe` value (no rate is evaluated). Public `get_rates_all` / `hi_rp_rm_pd` are unchanged.
"""
function _hi_pd_only(z, pops::HIPopulationSplines, t::AtomicRateTable, lnBitot::AbstractMatrix, lv::HIPDELevels, Tg, Te, Ne, NH)
    lx, ly, a, b, db = _setup(t, Tg, Te)
    n = n_resolved(t)
    lq = log_qnl_qe(t, 1, Tg)
    TA = db ? typeof(lq) : typeof(zero(a[1] * b[1]) + b[1] * (a[1] * t.A[1]) + lq)      # eltype of `A` in get_rates_all
    TB = typeof(zero(a[1]) + a[1] * t.B[1])                                              # eltype of `B`
    T = promote_type(TA, TB, typeof(Ne), typeof(NH), typeof(hi_Xi(pops, z, 0)))
    pd = Vector{T}(undef, n)
    for nr in 1:n
        Rm = exp(a[1] * lnBitot[lx, nr] + a[2] * lnBitot[lx + 1, nr] + a[3] * lnBitot[lx + 2, nr] + a[4] * lnBitot[lx + 3, nr])
        pd[nr] = lv.index[nr] == 1 ? Rm / (lv.A2s1s + Rm) : Rm / (lv.A21[nr] + Rm)
    end
    return pd
end

"""The `pd` and `Dnem` coefficient splines of the PDE (native `HI_pd_eff_Data`, `HI_Dnem_eff_Data`): natural cubic splines of `ln pd_i` and `Dnem_i` on the coefficient grid `z`."""
struct HIPDECoefficientSplines{S}
    z::Vector{Float64}
    lnpd::Vector{S}
    Dnem::Vector{S}
end

"""
    hi_pde_coefficients(rows, pops, cosmos, table, lnBitot, levels; zs, ze) -> HIPDECoefficientSplines

Native `set_up_splines_for_HI_pd_Rp_effective(ze, zs, cosmos, HA, X_Data, 1, 1)`: grid = `init_xarr(ze/1.0001, zs*1.0001, nz, linear)` with `nz` = number of stored rows with
`z >= ze/1.0001` (the native collection keeps the zero placeholders of the rows above `zs*1.0001`; the values are then overwritten by the linear grid); at each grid node `Tg = TCMB`, `Te = Tg rho(z)`, `NH`, `H` from the Cosmos accessors (4b), `pd` from
[`hi_rp_rm_pd`](@ref) (computed by its private pd-only kernel `_hi_pd_only`, bitwise identical) and `Dnem`: 2s `N2s/N1s - e^{-x}`; np (n >= 2) `(N_np/3/N1s - e^{-x}) (1 + (1/pd - 1) P_S)`, `P_S = (1 - e^{-tau_S})/tau_S`,
`tau_S = A21 lambda21^3 / (8 pi H) (3 N1s - N_np)`; ns, nd (n > 2) `N_i/N1s/(2l+1) - e^{-x}`; `x = h nu21 / (kB Tg)` with the level's `Dnu_1s`.
"""
function hi_pde_coefficients(rows::AbstractMatrix, pops::HIPopulationSplines, cosmos, t::AtomicRateTable, lnBitot::AbstractMatrix, lv::HIPDELevels;
                             zs::Real, ze::Real, h_kb::Float64 = 4.7992373449498863e-11, pi_::Float64 = 3.1415926535897931)
    lo = ze / 1.0001; hi = zs * 1.0001
    # native: the in-range redshifts are written to positions nmax-1-k and only the LEADING zeros (rows below the range) are erased, so rows ABOVE the range
    # still count: nz = number of stored rows with z >= ze/1.0001
    nz = count(r -> _primal(r) >= lo, view(rows, :, 1))
    zg = init_xarr_linear(lo, hi, nz)
    n = length(lv.index)
    T = promote_type(typeof(hi_Xi(pops, zg[1], 0)), typeof(cosmos_H(cosmos, zg[1])), typeof(cosmos_NH(cosmos, zg[1])))
    lnpd = [Vector{T}(undef, nz) for _ in 1:n]; Dn = [Vector{T}(undef, nz) for _ in 1:n]
    for k in 1:nz
        z = zg[k]
        Tg = cosmos_TCMB(cosmos, z); Te = Tg * hi_rho(pops, z); NH = cosmos_NH(cosmos, z); Hz = cosmos_H(cosmos, z)
        pd = _hi_pd_only(z, pops, t, lnBitot, lv, Tg, Te, NH * hi_Xe(pops, z), NH)
        N1s = NH * hi_Xi(pops, z, 0)
        for m in 1:n
            i = lv.index[m]
            ex = exp(-h_kb * lv.Dnu_1s[m] / Tg)
            Ni = NH * hi_Xi(pops, z, i)
            if i == 1
                Dn[m][k] = Ni / N1s - ex
            elseif lv.l[m] == 1
                nL = Ni / 3.0 / N1s
                tauS = lv.A21[m] * lv.lambda21[m]^3 / (8.0 * pi_ * Hz) * (N1s * 3.0 - Ni)
                PS = (1.0 - exp(-tauS)) / tauS
                Dn[m][k] = (nL - ex) * (1.0 + (1.0 / pd[m] - 1.0) * PS)
            else
                Dn[m][k] = Ni / N1s / (2.0 * lv.l[m] + 1.0) - ex
            end
            lnpd[m][k] = log(pd[m])
        end
    end
    return HIPDECoefficientSplines(zg, [natural_cubic_spline(zg, lnpd[m]) for m in 1:n], [natural_cubic_spline(zg, Dn[m]) for m in 1:n])
end

"""Native `calc_HI_pd_i_splines_effective(z, i)` for resolved position `m` (1-based)."""
hi_pd(c::HIPDECoefficientSplines, z, m::Integer) = exp(spline_eval_native(c.lnpd[m], z))
"""Native `calc_HI_Dnem_i_splines_effective(z, i)` for resolved position `m` (1-based)."""
hi_Dnem(c::HIPDECoefficientSplines, z, m::Integer) = spline_eval_native(c.Dnem[m], z)
