# Chunk 4b helpers: parse the native history / Hubble / accessor text fixtures and build the explicit accessor object.
using CosmoRec
using SHA

const FIX4B_HIST = joinpath(@__DIR__, "fixtures", "native_recfast_history.txt")
const FIX4B_HUB = joinpath(@__DIR__, "fixtures", "native_hubble_input.txt")
const FIX4B_ACC = joinpath(@__DIR__, "fixtures", "native_cosmos_accessors.txt")
const FIX4B_SHA = Dict(FIX4B_HIST => "f893af2d8a1f52c0278e422a5490616acb347081a10c09981de0b33020e97716",
                       FIX4B_HUB => "bcfa199b78d756f433202f0997bf4315761b28e6193b363745fdbc25c1f2ed99",
                       FIX4B_ACC => "b34bffd9234a2fa334533ce5defd3af06239e0a56c97587374fce39e44622ee9")

function read_history4b(path = FIX4B_HIST)
    rows = [parse.(Float64, split(l, '\t')) for l in eachline(path) if !isempty(l) && !startswith(l, "#")]
    return reduce(hcat, rows)'          # M x 7: z, Xe_H, Xe_He, Xe, dXe, dXe_H, TM
end

function read_hubble4b(path = FIX4B_HUB)
    rows = [parse.(Float64, split(l)) for l in eachline(path) if !isempty(l) && !startswith(l, "#")]
    return reduce(hcat, rows)'          # N x 2: z (descending), H [1/s]
end

function read_accessors4b(path = FIX4B_ACC)
    cos2 = Dict{String,Float64}(); const_ = Dict{String,Float64}(); acc = Vector{Dict{String,Float64}}(); sbt = Vector{Dict{String,Float64}}()
    kv(f, start = 2) = Dict(String(f[k]) => parse(Float64, f[k + 1]) for k in start:2:length(f))
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t')
        f[1] == "COS2" && merge!(cos2, kv(f))
        f[1] == "CONST" && merge!(const_, kv(f))
        f[1] == "ACC" && push!(acc, merge(Dict("z" => parse(Float64, f[2])), kv(f, 3)))
        f[1] == "SBT" && push!(sbt, merge(Dict("T" => parse(Float64, f[2])), kv(f, 3)))
    end
    return (cos2 = cos2, const_ = const_, acc = acc, sbt = sbt)
end

const HIST4B = read_history4b()
const HUB4B = read_hubble4b()
const ACC4B = read_accessors4b()

function constants4b(cos2 = ACC4B.cos2, k = ACC4B.const_)
    return CosmosConstants(cos2["Y_p"], cos2["T_CMB0"], cos2["fac_mHemH"], cos2["zsRe"], cos2["Nb0"], cos2["H0"], cos2["O_k"], cos2["O_L"], cos2["O_m"], cos2["O_rel"],
                           cos2["rho_g_gr"], cos2["sigT0"], 2.99792458e+10, cos2["hPlanck"], cos2["kBoltz"], cos2["mElect"], k["EH_inf_ergs"], k["me_mp"], k["me_malp"], 3.93933e-18, k["PI"], 3500.0)
end

splines4b(c = constants4b(), H = HIST4B) = recfast_splines(c, H[:, 1], H[:, 2], H[:, 3], H[:, 4], H[:, 5], H[:, 6], H[:, 7])
accessors4b(; hubble = true) = CosmosAccessors(constants4b(), splines4b(), hubble ? hubble_table(HUB4B[:, 1], HUB4B[:, 2]) : nothing)
