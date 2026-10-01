# Chunk 3b benchmark: helium base RHS primitives (explicit args, evals = 1). BenchmarkTools only; no @time/@elapsed.
#   julia --project=<test_env with CosmoRec + BenchmarkTools> benchmark/chunk3b_benchmarks.jl
using BenchmarkTools, CosmoRec
include(joinpath(@__DIR__, "..", "test", "chunk3b_helpers.jl"))
G, lgTg, cases, W = read_fixture3b()
const TAB = window_table(lgTg, W)
c = first(q for q in cases if q.label == "phys_z2300")
r = get_helium_rates(TAB, c.Tg); dX = zeros(7); helium_base_rhs!(dX, TAB, c.Tg, c.Xe, c.NH, c.Hz, c.X, c.fHe)   # warm-up
println("get_helium_rates (allocating A, B, R)")
display(@benchmark get_helium_rates($TAB, $(c.Tg)) samples = 2000 evals = 1)
println("\nhelium_base_rhs! from fixed rates (in place, no allocation expected)")
display(@benchmark helium_base_rhs!(d, $(c.Tg), $(c.Xe), $(c.NH), $(c.Hz), $(c.X), $(c.fHe), $(r.A), $(r.B), $(r.R)) setup = (d = zeros(7)) samples = 2000 evals = 1)
println("\nhelium_base_rhs! through the table (rate lookup + RHS)")
display(@benchmark helium_base_rhs!(d, $TAB, $(c.Tg), $(c.Xe), $(c.NH), $(c.Hz), $(c.X), $(c.fHe)) setup = (d = zeros(7)) samples = 2000 evals = 1)
println()
