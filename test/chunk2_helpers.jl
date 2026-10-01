# Test helpers for Chunk 2 (Base + stdlib only): parse the plain-text native fixtures.
using CosmoRec

const FIXDIR2 = joinpath(@__DIR__, "fixtures")

function parse_header2(path)
    comments = String[]
    for line in eachline(path)
        startswith(line, "#") ? push!(comments, line) : break
    end
    return comments
end

function fixture_states2(comments)
    n = Int[]; l = Int[]; gw = Float64[]; nuion = Float64[]; mu = Float64[]
    for c in comments
        m = match(r"^# state \d+ n=(\d+) l=(\d+) HI_index=\d+ gw=(\S+) nuion=(\S+) mu_red=(\S+) NTg=", c)
        m === nothing && continue
        push!(n, parse(Int, m[1])); push!(l, parse(Int, m[2]))
        push!(gw, parse(Float64, m[3])); push!(nuion, parse(Float64, m[4])); push!(mu, parse(Float64, m[5]))
    end
    return ResolvedStates(n, l, gw, nuion, mu)
end

function fixture_constants2(comments)
    c = first(filter(s -> startswith(s, "# constants"), comments))
    d = Dict(String(m[1]) => parse(Float64, m[2]) for m in eachmatch(r"(\w+)=(\S+)", split(c, ":"; limit = 2)[2]))
    return d
end

# table window -> (AtomicRateTable, original row indices)
function read_table_window2(path = joinpath(FIXDIR2, "native_hrates_table_window.txt"))
    comments = parse_header2(path)
    states = fixture_states2(comments)
    nres = length(states.n)
    lgrho = Float64[]; rows = Int[]; lg = Dict{Int,Float64}()
    Bd = Dict{Tuple{Int,Int},Float64}(); Rd = Dict{Tuple{Int,Int},Vector{Float64}}(); Ad = Dict{Tuple{Int,Int},Vector{Float64}}()
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, '\t')
        if f[1] == "lgrho"
            lgrho = parse.(Float64, f[2:end])
        elseif f[1] == "lgTg"
            r = parse(Int, f[2]); push!(rows, r); lg[r] = parse(Float64, f[3])
        elseif f[1] == "B"
            Bd[(parse(Int, f[2]), parse(Int, f[3]) + 1)] = parse(Float64, f[4])
        elseif f[1] == "R"
            Rd[(parse(Int, f[2]), parse(Int, f[3]) + 1)] = parse.(Float64, f[4:end])
        elseif f[1] == "A"
            Ad[(parse(Int, f[2]), parse(Int, f[3]) + 1)] = parse.(Float64, f[4:end])
        end
    end
    N, M, neq = length(rows), length(lgrho), length(Rd[(rows[1], 1)])
    B = [Bd[(rows[i], m)] for i in 1:N, m in 1:nres]
    R = [Rd[(rows[i], m)][k] for i in 1:N, m in 1:nres, k in 1:neq]
    A = [Ad[(rows[i], m)][j] for i in 1:N, j in 1:M, m in 1:nres]
    return AtomicRateTable([lg[r] for r in rows], lgrho, B, R, A, states), rows, comments
end

# query records: label => (Tg, Te, A, B, R[nres, neq])
function read_queries2(path = joinpath(FIXDIR2, "native_hrates_get_rates.txt"))
    recs = Vector{Any}()
    cur = nothing
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, '\t')
        id = parse(Int, f[1]); label = String(f[2])
        Tg = parse(Float64, f[3]); Te = parse(Float64, f[4]); m = parse(Int, f[5]) + 1
        A = parse(Float64, f[8]); B = parse(Float64, f[9]); R = parse.(Float64, f[10:end])
        if cur === nothing || cur.id != id
            cur = (id = id, label = label, Tg = Tg, Te = Te, A = Float64[], B = Float64[], R = Vector{Vector{Float64}}())
            push!(recs, cur)
        end
        push!(cur.A, A); push!(cur.B, B); push!(cur.R, R)
    end
    return recs
end
