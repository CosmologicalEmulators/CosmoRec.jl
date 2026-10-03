# Chunk 4d helpers: parse the plain-text native sampled-HeI switch fixture (test/fixtures/native_helium_switch.txt).
using CosmoRec
using SHA

const FIX4D = joinpath(@__DIR__, "fixtures", "native_helium_switch.txt")
const FIX4D_SHA256 = "afa4b768ac5f6ec347a0629649994a056e60f2d0b183fa6d435dc390225d71ae"

struct HSW4D
    id::Int
    zs::Float64
    delta::Float64
    Xi0::Float64
    fHe::Float64
    diff::Float64
    cond::Bool
    flag_pre::Int
    Xpre::Vector{Float64}
    ypre::Vector{Float64}
    flag_post::Int
    neq_y::Int
    Xpost::Vector{Float64}
    Xi::Vector{Float64}
    ypost::Vector{Float64}
end

function read_fixture4d(path = FIX4D)
    cfg = Dict{String,Float64}(); rows = HSW4D[]
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t')
        if f[1] == "HSWCFG"
            for k in 2:2:length(f); cfg[String(f[k])] = parse(Float64, f[k + 1]); end
        elseif f[1] == "HSW"
            v(i) = parse(Float64, f[i])
            push!(rows, HSW4D(parse(Int, f[2]), v(3), v(4), v(6), v(8), v(10), f[12] == "1", parse(Int, f[14]), parse.(Float64, f[16:30]), parse.(Float64, f[32:43]),
                              parse(Int, f[45]), parse(Int, f[47]), parse.(Float64, f[49:63]), parse.(Float64, f[65:71]), parse.(Float64, f[73:end])))
        end
    end
    return cfg, rows
end

const FX4D_CFG, FX4D = read_fixture4d()
