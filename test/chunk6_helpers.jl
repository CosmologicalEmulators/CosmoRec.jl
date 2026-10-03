# Phase 6 helpers: native PDE-coefficient fixture, HI Bitot table (from COSMOREC_NATIVE_DATA_DIR), production range.
using CosmoRec
using SHA
@isdefined(RM5) || include("chunk5_helpers.jl")

const FIX6A = joinpath(@__DIR__, "fixtures", "native_hi_pde_coefficients.txt")
const FIX6A_SHA256 = "98aa73169cfc42044f645d5fd75d49db0b73177a3548c39b2b2ce96543e8b855"
const ZS6 = 2500.0      # min(Diff_corr_zmax f_t 1.25, zstart/1.002), Diff_corr_zmax = 2000
const ZE6 = 500.0       # max(Diff_corr_zmin f_t, zend 1.002), Diff_corr_zmin = 500
const HDIR6 = joinpath(DATADIR5, "Effective_Rates.HI", "Effective_Rate_Tables.nS_3")
const LNBITOT6 = load_hi_bitot_table(HDIR6, 500, NATIVE_HYDROGEN_RESOLVED)
const HTAB6 = EFF5.htable

function read_fixture6a(path = FIX6A)
    d = Dict{String,Any}("POP" => Vector{Vector{Float64}}(), "GRA" => Vector{Vector{Float64}}(), "RPD" => Vector{Vector{Float64}}(), "PDN" => Vector{Vector{Float64}}(), "LVL6" => Vector{Vector{Float64}}())
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t')
        nums = [x for x in f[2:end] if tryparse(Float64, x) !== nothing]
        haskey(d, f[1]) && push!(d[f[1]], parse.(Float64, nums))
        f[1] == "CFG6" && (d["CFG6"] = Dict(String(f[k]) => parse(Float64, f[k + 1]) for k in 2:2:length(f)))
    end
    return d
end
const FX6A = read_fixture6a()
"""Native-layout rows `[z, Xe, X1s, 2s, 2p, 3s, 3p, 3d, rho]` from a Julia 5b pass."""
pass_rows6(pass) = reduce(vcat, [[pass.z[m] pass.observables[m]'] for m in eachindex(pass.z)])
# Source-derived residual scale of Dnem (HI_pd_Rp_splines_effective.cpp): Dnem is a DIFFERENCE of the population ratio term and exp(-x), so near equilibrium it cancels;
# errors are measured relative to the sum of the magnitudes of its two terms (times the escape factor for np states), never relative to Dnem itself.
function dnem_scale6(pops, cosmos, z, m, pd; lv = NATIVE_HI_PDE_LEVELS, h_kb = 4.7992373449498863e-11)
    Tg = cosmos_TCMB(cosmos, z); NH = cosmos_NH(cosmos, z); Hz = cosmos_H(cosmos, z)
    N1s = NH * hi_Xi(pops, z, 0); Ni = NH * hi_Xi(pops, z, lv.index[m]); ex = exp(-h_kb * lv.Dnu_1s[m] / Tg)
    if lv.index[m] == 1
        return Ni / N1s + ex
    elseif lv.l[m] == 1
        tauS = lv.A21[m] * lv.lambda21[m]^3 / (8.0 * pi * Hz) * (N1s * 3.0 - Ni); PS = (1.0 - exp(-tauS)) / tauS
        return (Ni / 3.0 / N1s + ex) * abs(1.0 + (1.0 / pd - 1.0) * PS)
    else
        return Ni / N1s / (2.0 * lv.l[m] + 1.0) + ex
    end
end
