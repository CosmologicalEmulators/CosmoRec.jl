# Test helpers for Chunk 3e: parse the plain-text native DPesc_coh fixture (test/fixtures/native_dpesc_coh.txt).
using CosmoRec
using SHA

const FIX3E = joinpath(@__DIR__, "fixtures", "native_dpesc_coh.txt")
const FIX3E_SHA256 = "bb82e20962a3b5abfaaa52dbcf111f1ecbb7fde3937688e94446707c63439b5b"

struct DPQuery
    id::Int
    label::String
    triplet::Bool
    code::String
    T::Float64
    eta::Float64      # eta_c (continuum parameter)
    tauS::Float64     # eta_S of the integral
    pd::Float64
    PS::Float64
    Pd::Float64
    Dpij::Float64
    Pesc::Float64
    corr::Float64     # call_DP_*(...) = Pesc - PS
    tab::Float64      # DP_interpol_S/T at the same query
end

function read_fixture3e(path = FIX3E)
    qs = DPQuery[]
    vx = NamedTuple[]; sx = Tuple{Float64,Float64}[]; vt = NamedTuple[]
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, '\t')
        if f[1] == "Q"
            v = parse.(Float64, f[[6, 7, 8, 9, 11, 13, 15, 17, 19, 21]])
            push!(qs, DPQuery(parse(Int, f[2]), String(f[3]), f[4] == "T", String(f[5]), v...))
        elseif f[1] == "VX"
            push!(vx, (ch = String(f[2]), Tm = parse(Float64, f[3]), x = parse(Float64, f[4]), phi = parse(Float64, f[6]), xi = parse(Float64, f[8])))
        elseif f[1] == "VT"
            push!(vt, (ch = String(f[2]), Tm = parse(Float64, f[3]), DnuT = parse(Float64, f[5]), a = parse(Float64, f[7]), xi1e4 = parse(Float64, f[9])))
        elseif f[1] == "SX"
            push!(sx, (parse(Float64, f[2]), parse(Float64, f[3])))
        end
    end
    return (queries = qs, vx = vx, vt = vt, sx = sx)
end

const FX3E = read_fixture3e()
