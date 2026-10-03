# Phase 7 helpers: native PDE-setup / PDE-coefficient fixtures and the ORIGINAL two-photon tables (located next to COSMOREC_NATIVE_DATA_DIR, never bundled).
using CosmoRec
using SHA
@isdefined(FX6A) || include("chunk6_helpers.jl")

const FIX7A = joinpath(@__DIR__, "fixtures", "native_hi_pde_setup.txt")
const FIX7A_SHA256 = "43e6f1e610da2fd4f05432442f32087212f807d4c9d33550be687180ac4b82aa"
const FIX7B = joinpath(@__DIR__, "fixtures", "native_hi_pde_def.txt")
const FIX7B_SHA256 = "9a7a418d08ab79854d3db3f8de7934442d644bb75d8c56572284d2cdc97c62c3"

"""ORIGINAL `Development/Line_profiles/two-photon-data` of the native tree whose `Rec_database` is COSMOREC_NATIVE_DATA_DIR (fails loudly if missing)."""
function two_photon_dir7(d = DATADIR5)
    p = joinpath(dirname(d), "Development", "Line_profiles", "two-photon-data")
    isdir(p) || error("missing native two-photon-data directory at $p (derived from COSMOREC_NATIVE_DATA_DIR = $d)")
    return p
end
const TPD7 = two_photon_dir7()

"""Records of a Phase 7 fixture: tag => list of records; each record keeps the numeric fields (labels dropped)."""
function read_fixture7(path)
    d = Dict{String,Vector{Vector{Float64}}}()
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, '\t')
        push!(get!(d, String(f[1]), Vector{Vector{Float64}}()), [parse(Float64, x) for x in f[2:end] if tryparse(Float64, x) !== nothing])
    end
    return d
end
"""Labelled fields of a single-line record (`TAG  key value  key value ...`)."""
function read_kv7(path, tag)
    for l in eachline(path)
        f = split(l, '\t')
        f[1] == tag && return Dict(String(f[k]) => f[k + 1] for k in 2:2:(length(f) - 1))
    end
    error("record $tag not found in $path")
end
"""Vector record `TAG  count  v1 ... vcount`."""
vec7(r::Vector{Float64}) = (n = Int(r[1]); @assert length(r) == n + 1; r[2:end])

const FX7A = read_fixture7(FIX7A)
const PROF7 = load_hi_profile_data(TPD7)        # production const_HI_A2s_1s = 8.2206 (checked against CON7)
const SETUP7 = hi_pde_setup(PROF7)

const FIX7C = joinpath(@__DIR__, "fixtures", "native_hi_pde_march.txt")
const FIX7C_SHA256 = "646abf27060b21ec141cd45d4c812dfe8838d796bb1a65344d14dadf73edfe8b"
const FIX8A = joinpath(@__DIR__, "fixtures", "native_hi_pde_integrals.txt")
const FIX8A_SHA256 = "9dd47147e65f202e82452dda361af39ce8d7b213d581561405a59d4736eed528"
