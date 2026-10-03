# Chunk 4a helpers: parse the plain-text native Saha-initialization fixture (test/fixtures/native_saha_init.txt).
using CosmoRec
using SHA

const FIX4A = joinpath(@__DIR__, "fixtures", "native_saha_init.txt")
const FIX4A_SHA256 = "d61d10b689e572e6848bd5361503e16825d9ecb32f79a41164b3e8cd886e7656"

function read_fixture4a(path = FIX4A)
    f64(x) = parse(Float64, x)
    d = Dict{String,Any}("SAHAIN" => Dict{Float64,Dict{String,Float64}}(), "LTEH" => Dict{Float64,Vector{Float64}}(), "LTEHE" => Dict{Float64,Vector{Float64}}(),
                         "XLI" => Dict{Float64,Vector{Float64}}(), "YSOL" => Dict{Float64,Vector{Float64}}(), "SP" => Dict{Float64,Dict{String,Float64}}(),
                         "LVLH" => Vector{Vector{Float64}}(), "LVLHE" => Vector{Vector{Float64}}(), "CONST" => Dict{String,Float64}())
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, '\t')
        t = f[1]
        if t == "SAHAIN"
            d[t][f64(f[2])] = Dict(String(f[k]) => f64(f[k + 1]) for k in 3:2:length(f))
        elseif t == "SP"
            d[t][f64(f[2])] = Dict(String(f[k]) => f64(f[k + 1]) for k in 3:2:length(f))
        elseif t in ("LTEH", "LTEHE", "XLI")
            d[t][f64(f[2])] = f64.(f[3:end])
        elseif t == "YSOL"
            d[t][f64(f[2])] = f64.(f[4:end])
        elseif t == "LVLH"
            push!(d[t], [f64(f[4]), f64(f[6]), f64(f[8]), f64(f[10]), f64(f[12])])      # n l Eion mu ME
        elseif t == "LVLHE"
            push!(d[t], [f64(f[5]), f64(f[7]), f64(f[9]), f64(f[11])])                  # gw Eion mu ME
        elseif t == "CONST"
            for k in 2:2:length(f); d[t][String(f[k])] = f64(f[k + 1]); end
        elseif t == "LI"
            d[t] = parse.(Int, f[2:end])
        elseif t == "RESHI" || t == "RESHE"
            d[t] = parse.(Int, f[2:end])
        elseif t == "COS"
            d["fHe"] = f64(f[3])
        end
    end
    return d
end

const FX4A = read_fixture4a()

# SahaInputs for the native state at redshift z (SAHAIN record)
function inputs4a(z)
    s = FX4A["SAHAIN"][z]
    return SahaInputs(FX4A["fHe"], s["Xe_clipped"], s["Xp_raw"], s["NH"], s["Te"], s["TCMB"], s["XHeII"], s["XHeI_ground"])
end
