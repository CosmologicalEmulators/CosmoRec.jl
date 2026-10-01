# Test helpers for Chunk 3c: parse the plain-text native H-I absorption fixture and build per-case sparse native tables.
using CosmoRec

const FIX3C = joinpath(@__DIR__, "fixtures", "native_hiabs_components.txt")

struct Case3c
    id::Int
    label::String
    z::Float64
    Tg::Float64
    NH::Float64
    Hz::Float64
    XH1s::Float64
    Xe::Float64
    X::Vector{Float64}
    eta::Float64
    tauS_S::Float64
    tauS_T::Float64
    efacS::Float64
    efacT::Float64
    pdS::Float64
    pdT::Float64
    Bitot::Float64
    fcorr::Float64
    inS::Bool
    inT::Bool
    Bwin::Union{Nothing,Tuple{Int,Vector{Float64},Vector{Float64}}}      # (lx0, lgTg[4], lnBitot[4])
    win::Dict{String,Vector{Vector{Float64}}}                             # "S"/"T" -> 4 sheet rows
    sheet0::Dict{String,Int}                                              # 0-based first sheet
    V::Dict{Int,Vector{Float64}}                                          # flag -> dpS dxS dpT dxT
    K::Dict{Tuple{Int,Int,Int},Tuple{Bool,Bool,Vector{Float64}}}          # (flag, spin, kind) -> (activeS, activeT, g)
end

"""Returns `(G, lnT, lgTg, F, cases)`."""
function read_fixture3c(path = FIX3C)
    G = Vector{Vector{String}}(); lnT = Float64[]; lgTg = Float64[]; F = Vector{Tuple{Float64,Float64}}()
    raw = Dict{Int,Dict{Symbol,Any}}()
    for line in eachline(path)
        (isempty(line) || startswith(line, "#") || !occursin(r"^[A-Z]\t", line)) && continue
        f = split(line, '\t')
        t = f[1]
        if t == "G"
            push!(G, f[2:end])
        elseif t == "T"
            lnT = parse.(Float64, f[3:end])
        elseif t == "L"
            lgTg = parse.(Float64, f[4:end])
        elseif t == "F"
            push!(F, (parse(Float64, f[2]), parse(Float64, f[3])))
        else
            id = parse(Int, f[2])
            if t == "C"
                raw[id] = Dict{Symbol,Any}(:label => String(f[3]), :c => parse.(Float64, f[4:end]), :win => Dict{String,Vector{Vector{Float64}}}(),
                                           :s0 => Dict{String,Int}(), :V => Dict{Int,Vector{Float64}}(), :K => Dict{Tuple{Int,Int,Int},Any}(), :B => nothing)
            elseif t == "I"
                raw[id][:I] = parse.(Float64, f[3:end])
            elseif t == "H"
                lx = parse(Int, f[3])
                raw[id][:B] = length(f) >= 12 ? (lx, parse.(Float64, f[5:8]), parse.(Float64, f[9:12])) : nothing
            elseif t == "S"
                raw[id][:s0][String(f[3])] = parse(Int, f[4])
            elseif t == "R"
                push!(get!(raw[id][:win], String(f[3]), Vector{Vector{Float64}}()), parse.(Float64, f[5:end]))
            elseif t == "V"
                raw[id][:V][parse(Int, f[3])] = parse.(Float64, f[4:end])
            elseif t == "K"
                v = parse.(Float64, f[3:end])
                raw[id][:K][(Int(v[1]), Int(v[2]), Int(v[3]))] = (v[4] == 1, v[5] == 1, v[6:end])
            end
        end
    end
    cases = Case3c[]
    for id in sort(collect(keys(raw)))
        r = raw[id]; c = r[:c]; I = r[:I]
        push!(cases, Case3c(id, r[:label], c[1], c[2], c[3], c[4], c[5], c[6], c[7:13], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8], I[9],
                            I[10] == 1, I[11] == 1, r[:B], r[:win], r[:s0], r[:V], r[:K]))
    end
    return G, lnT, lgTg, F, cases
end

grow3c(G, key) = first(r for r in G if r[1] == key)

"""Sparse per-case native DP table: full native `lnT`; only the 4 captured sheets (axis ends, spacing and the stencil window) are populated, NaN elsewhere."""
function window_dp_table(lnT, c::Case3c)
    isempty(c.win) && return nothing
    nan = fill(NaN, 31)
    dummy = DPSheet(nan, nan, nan, fill(NaN, 31, 31), fill(NaN, 31, 31))
    sheets = Vector{DPSheet}(fill(dummy, length(lnT)))
    for k in 1:4
        eta = fill(NaN, 31); tS = fill(NaN, 31); tT = fill(NaN, 31); DS = fill(NaN, 31, 31); DT = fill(NaN, 31, 31)
        for (ch, tax, DP) in (("S", tS, DS), ("T", tT, DT))
            haskey(c.win, ch) || continue
            w = c.win[ch][k]
            neta = Int(w[1]); ntau = Int(w[2]); ie = Int(w[3]); it = Int(w[4])
            @assert neta == 31 && ntau == 31
            eta[1] = w[5]; eta[2] = w[6]; eta[end] = w[7]; tax[1] = w[8]; tax[2] = w[9]; tax[end] = w[10]
            eta[ie + 1:ie + 4] = w[11:14]; tax[it + 1:it + 4] = w[15:18]
            DP[ie + 1:ie + 4, it + 1:it + 4] = permutedims(reshape(w[19:34], 4, 4))
        end
        j = (haskey(c.sheet0, "S") ? c.sheet0["S"] : c.sheet0["T"]) + k
        sheets[j] = DPSheet(eta, tS, tT, DS, DT)
    end
    return DPTable(lnT, sheets)
end

"""Sparse Bitot series on the full native HeI `lgTg` grid."""
function window_bitot(lgTg, c::Case3c)
    ln = fill(NaN, length(lgTg))
    if c.Bwin !== nothing
        lx, lg, lb = c.Bwin
        @assert lgTg[lx + 1:lx + 4] == lg
        ln[lx + 1:lx + 4] = lb
    end
    return BitotSeries(lgTg, ln)
end

fcorr_spline3c(F) = FcorrSpline([f[1] for f in F], [f[2] for f in F])
