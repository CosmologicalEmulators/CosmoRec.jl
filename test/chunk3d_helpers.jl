# Test helpers for Chunk 3d: parse the plain-text native fcn_effective fixture and build sparse (NaN outside the captured windows) native tables.
using CosmoRec
using SHA

isdefined(@__MODULE__, :read_table_window2) || include(joinpath(@__DIR__, "chunk2_helpers.jl"))
isdefined(@__MODULE__, :window_table) || include(joinpath(@__DIR__, "chunk3b_helpers.jl"))

const FIX3D = joinpath(@__DIR__, "fixtures", "native_fcn_effective.txt")
const FIX3D_SHA256 = "9f21e5e451a9eee5c83fc9df460bdf2f311cf363e9a7b5b8fb6f7c347cfedc34"

struct State3d
    id::Int
    label::String
    z::Float64
    Tg::Float64
    NH::Float64
    Hz::Float64
    fHe::Float64
    X::Vector{Float64}
    Xe::Float64
    Xp::Float64
    XHeII::Float64
    g::Dict{String,Vector{Float64}}   # "on", "on_cached", "off"
    dp_in::Bool
    eta::Float64
end

struct Fixture3d
    config::Dict{String,String}
    states::Vector{State3d}
    htable::AtomicRateTable
    hetable::HeliumRateTable
    dp::DPTable
    bitot::BitotSeries
    fc::FcorrSpline
end

function read_fixture3d(path = FIX3D)
    comments = parse_header2(path)
    rs = fixture_states2(comments)
    cfg = Dict{String,String}()
    bk = Dict{Int,Vector{String}}(); inn = Dict{Int,Vector{String}}(); fr = Dict{Int,Vector{Float64}}()
    go = Dict{Tuple{Int,String},Vector{Float64}}(); di = Dict{Int,Vector{Float64}}()
    lgrho = Float64[]; lgTgfull = Float64[]; Hrows = Int[]; Hlg = Dict{Int,Float64}()
    Bd = Dict{Tuple{Int,Int},Float64}(); Rd = Dict{Tuple{Int,Int},Vector{Float64}}(); Ad = Dict{Tuple{Int,Int},Vector{Float64}}()
    TH = Float64[]; W = Dict{Tuple{Int,Int},Tuple{Float64,Float64,Vector{Float64}}}()
    lnT = Float64[]; Bwin = Dict{Int,Tuple{Int,Vector{Float64},Vector{Float64}}}()
    s0 = Dict{Tuple{Int,String},Int}()
    F = Tuple{Float64,Float64}[]
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, '\t')
        t = f[1]
        if t == "CF"; cfg[String(f[2])] = String(f[3])
        elseif t == "BK"; bk[parse(Int, f[2])] = f[3:end]
        elseif t == "IN"; inn[parse(Int, f[2])] = f[3:end]
        elseif t == "FR"; fr[parse(Int, f[2])] = parse.(Float64, f[3:end])
        elseif t == "GO" && f[2] != "HeI_f_t_f_b"; go[(parse(Int, f[2]), String(f[3]))] = parse.(Float64, f[4:end])
        elseif t == "lgrho"; lgrho = parse.(Float64, f[2:end])
        elseif t == "lgTgfull"; lgTgfull = parse.(Float64, f[2:end])
        elseif t == "lgTg"; r = parse(Int, f[2]); push!(Hrows, r); Hlg[r] = parse(Float64, f[3])
        elseif t == "B"; Bd[(parse(Int, f[2]), parse(Int, f[3]) + 1)] = parse(Float64, f[4])
        elseif t == "R"; Rd[(parse(Int, f[2]), parse(Int, f[3]) + 1)] = parse.(Float64, f[4:end])
        elseif t == "A"; Ad[(parse(Int, f[2]), parse(Int, f[3]) + 1)] = parse.(Float64, f[4:end])
        elseif t == "TH"; TH = parse.(Float64, f[3:end])
        elseif t == "W"; v = parse.(Float64, f[5:end]); W[(parse(Int, f[3]) + 1, parse(Int, f[4]))] = (v[1], v[2], v[3:6])
        elseif t == "DT"; lnT = parse.(Float64, f[3:end])
        elseif t == "DI"; di[parse(Int, f[2])] = parse.(Float64, f[3:end])
        elseif t == "DH"; length(f) >= 12 && (Bwin[parse(Int, f[2])] = (parse(Int, f[3]), parse.(Float64, f[5:8]), parse.(Float64, f[9:12])))
        elseif t == "DS"; s0[(parse(Int, f[2]), String(f[3]))] = parse(Int, f[4])
        elseif t == "F"; push!(F, (parse(Float64, f[2]), parse(Float64, f[3])))
        end
    end
    Ng = length(lgTgfull); M = length(lgrho); nres = length(rs.n); neq = length(Rd[(first(Hrows), 1)])
    Bm = fill(NaN, Ng, nres); Rm = fill(NaN, Ng, nres, neq); Am = fill(NaN, Ng, M, nres)
    lg = copy(lgTgfull)
    for r in Hrows
        @assert lgTgfull[r + 1] == Hlg[r]
    end
    for ((r, m), v) in Bd; Bm[r + 1, m] = v; end
    for ((r, m), v) in Rd; Rm[r + 1, m, :] = v; end
    for ((r, m), v) in Ad; Am[r + 1, :, m] = v; end
    htable = AtomicRateTable(lg, lgrho, Bm, Rm, Am, rs)
    hetable = window_table(TH, W)
    # DP table: NaN everywhere except the captured axis-ends/spacing and the 4x4 stencil windows
    nan = fill(NaN, 31)
    sheets = [DPSheet(copy(nan), copy(nan), copy(nan), fill(NaN, 31, 31), fill(NaN, 31, 31)) for _ in lnT]
    for line in eachline(path)
        startswith(line, "DR\t") || continue
        f = split(line, '\t')
        id = parse(Int, f[2]); ch = String(f[3]); k = parse(Int, f[4])
        w = parse.(Float64, f[5:end])
        j = s0[(id, ch)] + k + 1
        sh = sheets[j]
        tax = ch == "S" ? sh.tauS : sh.tauT; DP = ch == "S" ? sh.DP_S : sh.DP_T
        ie = Int(w[3]); it = Int(w[4])
        sh.eta[1] = w[5]; sh.eta[2] = w[6]; sh.eta[end] = w[7]; tax[1] = w[8]; tax[2] = w[9]; tax[end] = w[10]
        sh.eta[ie + 1:ie + 4] = w[11:14]; tax[it + 1:it + 4] = w[15:18]
        DP[ie + 1:ie + 4, it + 1:it + 4] = permutedims(reshape(w[19:34], 4, 4))
    end
    dp = DPTable(lnT, sheets)
    lnB = fill(NaN, length(TH))
    for (id, (lx, lgw, lbw)) in Bwin
        @assert TH[lx + 1:lx + 4] == lgw
        lnB[lx + 1:lx + 4] = lbw
    end
    bitot = BitotSeries(TH, lnB)
    fc = FcorrSpline([x[1] for x in F], [x[2] for x in F])
    states = State3d[]
    for id in sort(collect(keys(bk)))
        b = bk[id]; x = parse.(Float64, inn[id][3:end]); d = di[id]
        push!(states, State3d(id, String(b[1]), parse(Float64, b[2]), parse(Float64, b[3]), parse(Float64, b[4]), parse(Float64, b[5]), parse(Float64, b[6]),
                              x, fr[id][1], fr[id][2], fr[id][3],
                              Dict(m => go[(id, m)] for m in ("on", "on_cached", "off")), d[8] == 1 && d[9] == 1, d[1]))
    end
    return Fixture3d(cfg, states, htable, hetable, dp, bitot, fc)
end

fixture_sha3d(path = FIX3D) = bytes2hex(sha256(read(path)))

model3d(fx::Fixture3d; hi_absorption = true) = EffectiveModel(fx.htable, fx.hetable, fx.dp, fx.bitot, fx.fc; hi_absorption)
background3d(s::State3d) = EffectiveBackground(s.Tg, s.NH, s.Hz, s.fHe)

const FX3D = read_fixture3d()
const MODEL3D_ON = model3d(FX3D; hi_absorption = true)
const MODEL3D_OFF = model3d(FX3D; hi_absorption = false)
const ABS_SLOTS3D = [1, 2, 8, 10, 13]   # Julia 1-based = native 0, 1, 7, 9, 12
