# Chunk 4a: explicit-input Saha per-level initialization and ODE-vector packing of ORIGINAL CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3):
#   Modules/HI_routines.cpp:150-175 (Set_Hydrogen_Levels_to_Saha), Modules/HeI_routines.cpp:140-171 (Set_HeI_Levels_to_Saha),
#   Development/Hydrogenic/Atom.cpp:1346-1359 (Ni_NeNc_LTE, Xi_Saha), Development/Helium.v1.1/HeI_Atom.cpp:960-972, 1332-1344, 1710-1722 (HeI Ni_NeNc_LTE, Xi_Saha),
#   Modules/ODEdef_CosmoRec.cpp:87-107 (copy_LI_to_ysol, flag_He >= 1 branch).
# CosmoRec is (c) J. Chluba et al.; use must be acknowledged and Chluba & Thomas 2010 (MNRAS 412, 748), Chluba & Sunyaev 2006 (A&A 446, 39) and
# Rubino-Martin, Chluba & Sunyaev 2008 cited (see docs/CHUNK4A_RESULTS.md NOTICE).
#
# SCOPE (nothing else is claimed): every scalar the native routines read from the Cosmos object (Xe_Seager, Xp, NH, Te, TCMB, NHeII/NH, NHeI/NH, fHe) is an EXPLICIT argument;
# atomic level data and constants are explicit structs. NOT here: Recfast++ history, Cosmos accessors/splines (z>=3500 branches included), H(z), the sampled-HeI
# state switch (copy_LI_to_ysol's flag_He = 0 branch), any ODE.
# State layout (native Data_Level_I, 1-based): X[1] = Xe, X[2:1+nH] = H levels (1s, 2s, 2p, 3s, 3p, 3d), X[2+nH:1+nH+nHe] = HeI levels, X[end] = rho = Te/Tg.

"""Constants of the Saha equations: `lambdac = const_lambdac` [cm], `kb_mec2 = const_kB/me/c^2` (`const_kb_mec2`), `kB` [erg/K], `pi`."""
struct SahaConstants
    lambdac::Float64
    kb_mec2::Float64
    kB::Float64
    pi::Float64
end

"""Native-printed values (physical_consts.h, harness record CONST)."""
const NATIVE_SAHA_CONSTANTS = SahaConstants(2.4263102175000002e-10, 1.6863720498875446e-10, 1.3806504000000002e-16, 3.1415926535897931)

"""
Per-level data of one species (native objects, harness records LVLH / LVLHE): `g` = statistical-weight numerator (hydrogen: `g_i = 2(2 l_i + 1)`, helium: `gw_i`),
`Eion` ionization energy [erg] of the level, `mu_red`, `ME_scale`.
"""
struct SahaLevels
    g::Vector{Float64}
    Eion::Vector{Float64}
    mu_red::Vector{Float64}
    ME_scale::Vector{Float64}
    function SahaLevels(g, Eion, mu_red, ME_scale)
        n = length(g)
        (n >= 1 && length(Eion) == n && length(mu_red) == n && length(ME_scale) == n) || throw(DimensionMismatch("SahaLevels: need equal non-empty vectors"))
        return new(g, Eion, mu_red, ME_scale)
    end
end

const NATIVE_SAHA_HYDROGEN = SahaLevels(
    [2.0, 2.0, 6.0, 2.0, 6.0, 10.0],   # 2(2l+1) for l = 0, 0, 1, 0, 1, 2
    [2.1786854208347312e-11, 5.446713552086828e-12, 5.446713552086828e-12, 2.420761578705257e-12, 2.420761578705257e-12, 2.420761578705257e-12],
    fill(0.99945567942448077, 6), ones(6))

const NATIVE_SAHA_HELIUM = SahaLevels(
    [1.0, 1.0, 3.0, 3.0, 1.0, 3.0, 5.0],
    [3.9393333611555142e-11, 6.3632241227213117e-12, 5.3983167367688556e-12, 7.6388133947291136e-12, 5.804955294985449e-12, 5.8051515388541983e-12, 5.805166720277502e-12],
    fill(0.99986292543644184, 7), ones(7))

"""Native `Gas_of_Atoms::Ni_NeNc_LTE(i, TM)` for 1-based level `i` (gc = 1): `g/2/gc lambdac^3 (2 pi kb_mec2 TM ME mu)^(-3/2) exp(Eion/kB/TM)`."""
function saha_lte_hydrogen(lv::SahaLevels, i::Integer, T, c::SahaConstants = NATIVE_SAHA_CONSTANTS)
    return lv.g[i] / 2.0 / 1.0 * c.lambdac^3 * (2.0 * c.pi * c.kb_mec2 * T * lv.ME_scale[i] * lv.mu_red[i])^(-1.5) * exp(lv.Eion[i] / c.kB / T)
end

"""Native `Electron_Level_HeI_*::Ni_NeNc_LTE(TM)`: `gw/4 lambdac^3 (2 pi kb_mec2 mu ME TM)^(-3/2) exp(Eion/kB/TM)` (same for singlet, triplet and non-j triplet levels)."""
function saha_lte_helium(lv::SahaLevels, i::Integer, T, c::SahaConstants = NATIVE_SAHA_CONSTANTS)
    return lv.g[i] / 4.0 * c.lambdac^3 * (2.0 * c.pi * c.kb_mec2 * lv.mu_red[i] * lv.ME_scale[i] * T)^(-1.5) * exp(lv.Eion[i] / c.kB / T)
end

"""
Explicit inputs of the native Saha initialization at one redshift (values of the Cosmos accessors; `z` itself is NOT an input and 4a has no z branch):
`fHe`, `Xe_Seager = cosmos.Xe_Seager(z)` (unclipped), `Xp_raw = cosmos.Xp(z)`, `NH` [cm^-3], `Te = cosmos.Te(z)` [K], `Tg = cosmos.TCMB(z)` [K],
`XHeII = cosmos.NHeII(z)/cosmos.NH(z)`, `XHeI1s = cosmos.NHeI(z)/cosmos.NH(z)`.
"""
struct SahaInputs{T}
    fHe::T
    Xe_Seager::T
    Xp_raw::T
    NH::T
    Te::T
    Tg::T
    XHeII::T
    XHeI1s::T
end
SahaInputs(a...) = SahaInputs(promote(a...)...)

"""
    saha_initial_state(inp, H, He, c = NATIVE_SAHA_CONSTANTS) -> X

Native `Set_Hydrogen_Levels_to_Saha(z)` followed by `Set_HeI_Levels_to_Saha(z)` into a fresh `Level_I.X` of length `1 + nH + nHe + 1`:
H: `X[1] = min(Xe_Seager, 1 + fHe)`, `Xp = min(Xp_raw, 1)`, `X_i = Xe Xp NH f_i(Te)`, ground `X_1s = 1 - Xp` (override), `X[end] = rho = 1`;
He (never touches `X[1]`): `X_i = Xe_Seager XHeII NH f^He_i(Tg)` with the UNCLIPPED `Xe_Seager`, ground `X = XHeI1s` (override). No other checks or smoothing (like native).
"""
function saha_initial_state(inp::SahaInputs, H::SahaLevels, He::SahaLevels, c::SahaConstants = NATIVE_SAHA_CONSTANTS)
    nH = length(H.g); nHe = length(He.g)
    T = promote_type(typeof(inp.fHe), typeof(c.lambdac), eltype(H.Eion), eltype(He.Eion))
    X = Vector{T}(undef, 1 + nH + nHe + 1)
    Xe = min(inp.Xe_Seager, 1.0 + inp.fHe)
    X[1] = Xe
    Xp = min(inp.Xp_raw, 1.0)
    for i in 1:nH
        X[1 + i] = Xe * Xp * inp.NH * saha_lte_hydrogen(H, i, inp.Te, c)
    end
    X[2] = 1.0 - Xp
    X[end] = one(T)
    for i in 1:nHe
        X[1 + nH + i] = inp.Xe_Seager * inp.XHeII * inp.NH * saha_lte_helium(He, i, inp.Tg, c)
    end
    X[2 + nH] = inp.XHeI1s
    return X
end

"""Native `get_HI_index(k)` for the production 3-shell atom (0-based level offsets of the 5 resolved states, harness record RESHI)."""
const NATIVE_RESHI = (1, 2, 3, 4, 5)
"""Native `get_HeI_index(k)` (harness record RESHE)."""
const NATIVE_RESHE = (1, 2, 3, 5)

"""
    pack_ysol(X, nH, nHe, resHI = NATIVE_RESHI, resHeI = NATIVE_RESHE) -> y

Native `copy_LI_to_ysol` with `flag_He >= 1` (`helium = true`, default) or the `flag_He = 0` branch of Chunk 4d (`helium = false`: `[rho, H 1s, resolved H]`): `y[1] = rho = X[end]`, `y[2] = X[H 1s]`, then the resolved H levels `X[H 1s + resHI[k]]`,
`He 1s`, and `X[He 1s + resHeI[k]]`. `res*` are the native 0-based level offsets; native `index_HI = 2`, `index_HeI = 2 + nH` in 1-based form.
"""
function pack_ysol(X::AbstractVector, nH::Integer, nHe::Integer, resHI = NATIVE_RESHI, resHeI = NATIVE_RESHE; helium::Bool = true)
    length(X) == 1 + nH + nHe + 1 || throw(DimensionMismatch("pack_ysol: length(X) = $(length(X)) != 1 + nH + nHe + 1 = $(1 + nH + nHe + 1)"))
    all(k -> 0 <= k < nH, resHI) || throw(BoundsError(resHI, 0))
    all(k -> 0 <= k < nHe, resHeI) || throw(BoundsError(resHeI, 0))
    iHI = 2; iHe = 2 + nH
    if !helium            # native `flag_He = 0` branch (after the sampled-HeI switch, Chunk 4d): only rho, H 1s and the resolved H levels
        y0 = Vector{eltype(X)}(undef, 2 + length(resHI))
        y0[1] = X[end]; y0[2] = X[iHI]
        for k in eachindex(resHI)
            y0[2 + k] = X[iHI + resHI[k]]
        end
        return y0
    end
    y = Vector{eltype(X)}(undef, 2 + length(resHI) + 1 + length(resHeI))
    y[1] = X[end]
    y[2] = X[iHI]
    for k in eachindex(resHI)
        y[2 + k] = X[iHI + resHI[k]]
    end
    y[2 + length(resHI) + 1] = X[iHe]
    for k in eachindex(resHeI)
        y[3 + length(resHI) + k] = X[iHe + resHeI[k]]
    end
    return y
end
