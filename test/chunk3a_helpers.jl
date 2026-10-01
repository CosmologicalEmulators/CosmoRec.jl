# Test helpers for Chunk 3a: parse the plain-text native RHS-component fixture.
using CosmoRec

const FIX3A = joinpath(@__DIR__, "fixtures", "native_hrhs_components.txt")

struct Case3a
    id::Int
    label::String
    z::Float64
    Tg::Float64
    Te::Float64
    rho::Float64
    Xe::Float64
    fHe::Float64
    Xp::Float64
    NH::Float64
    Hz::Float64
    X::Vector{Float64}
    A::Vector{Float64}
    B::Vector{Float64}
    R::Matrix{Float64}
    before::Dict{String,Vector{Float64}}
    after::Dict{String,Vector{Float64}}
end

function read_fixture3a(path = FIX3A)
    G = Vector{Vector{String}}()
    cases = Dict{Int,Any}()
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, '\t')
        if f[1] == "G"
            push!(G, f[2:end])
        elseif f[1] == "C"
            id = parse(Int, f[2])
            v = parse.(Float64, f[4:end]); cases[id] = Dict{Symbol,Any}(:label => f[3], :v => v, :A => zeros(5), :B => zeros(5), :R => zeros(5, 5), :b => Dict{String,Vector{Float64}}(), :a => Dict{String,Vector{Float64}}())
        elseif f[1] == "RT"
            c = cases[parse(Int, f[2])]; m = parse(Int, f[3]) + 1
            v = parse.(Float64, f[4:end]); c[:A][m] = v[1]; c[:B][m] = v[2]; c[:R][m, :] = v[3:7]
        elseif f[1] == "K"
            c = cases[parse(Int, f[2])]; v = parse.(Float64, f[4:end]); nb = length(v) ÷ 2
            c[:b][f[3]] = v[1:nb]; c[:a][f[3]] = v[(nb + 1):end]
        end
    end
    out = Case3a[]
    for id in sort(collect(keys(cases)))
        c = cases[id]; v = c[:v]
        push!(out, Case3a(id, c[:label], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10:15], c[:A], c[:B], c[:R], c[:b], c[:a]))
    end
    return G, out
end

grow(G, key) = first(r for r in G if r[1] == key)
