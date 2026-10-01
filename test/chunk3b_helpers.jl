# Test helpers for Chunk 3b: parse the plain-text native helium RHS-component fixture.
using CosmoRec

const FIX3B = joinpath(@__DIR__, "fixtures", "native_herhs_components.txt")
const NLHE = 7

struct Case3b
    id::Int
    label::String
    z::Float64
    Tg::Float64
    Xe::Float64
    fHe::Float64
    NH::Float64
    Hz::Float64
    XHeII::Float64
    X::Vector{Float64}
    A::Vector{Float64}
    B::Vector{Float64}
    R::Matrix{Float64}
    lx::Int
    after::Dict{String,Vector{Float64}}
end

"""Returns `(G, lgTg, cases, W)`; `W[(m, row0)] = (lgTg, Blog, Rlog)` with 0-based native table row, `m` 1-based resolved state."""
function read_fixture3b(path = FIX3B)
    G = Vector{Vector{String}}(); lgTg = Float64[]
    cases = Dict{Int,Any}(); W = Dict{Tuple{Int,Int},Tuple{Float64,Float64,Vector{Float64}}}()
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, '\t')
        if f[1] == "G"
            push!(G, f[2:end])
        elseif f[1] == "T"
            lgTg = parse.(Float64, f[6:end])
        elseif f[1] == "C"
            v = parse.(Float64, f[4:end])
            cases[parse(Int, f[2])] = Dict{Symbol,Any}(:label => f[3], :v => v, :A => zeros(4), :B => zeros(4), :R => zeros(4, 4), :a => Dict{String,Vector{Float64}}(), :lx => -1)
        elseif f[1] == "RT"
            c = cases[parse(Int, f[2])]; m = parse(Int, f[3]) + 1
            v = parse.(Float64, f[4:end]); c[:A][m] = v[1]; c[:B][m] = v[2]; c[:R][m, :] = v[3:6]
        elseif f[1] == "L"
            cases[parse(Int, f[2])][:lx] = parse(Int, f[3])
        elseif f[1] == "W"
            v = parse.(Float64, f[5:end]); W[(parse(Int, f[3]) + 1, parse(Int, f[4]))] = (v[1], v[2], v[3:6])
        elseif f[1] == "K"
            cases[parse(Int, f[2])][:a][f[3]] = parse.(Float64, f[4:end])
        end
    end
    out = Case3b[]
    for id in sort(collect(keys(cases)))
        c = cases[id]; v = c[:v]
        push!(out, Case3b(id, c[:label], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8:14], c[:A], c[:B], c[:R], c[:lx], c[:a]))
    end
    return G, lgTg, out, W
end

grow3b(G, key) = first(r for r in G if r[1] == key)

"""Sparse native table: full native `lgTg` grid; B and R filled (native log values) only where the native stencil windows were captured, NaN elsewhere."""
function window_table(lgTg, W)
    N = length(lgTg)
    B = fill(NaN, N, 4); R = fill(NaN, N, 4, 4)
    for ((m, row0), (lg, b, r)) in W
        @assert lgTg[row0 + 1] == lg
        B[row0 + 1, m] = b; R[row0 + 1, m, :] = r
    end
    atom = NATIVE_HELIUM_ATOM
    return HeliumRateTable(lgTg, B, R, atom.gw[atom.res_index], atom.nu_ion[atom.res_index], atom.mu_red)
end

"""Native sentinel 'before' vector of the harness."""
sentinel3b() = [(-1)^(i % 2) * 1.25e-3 / 2^i for i in 0:(NLHE - 1)]
