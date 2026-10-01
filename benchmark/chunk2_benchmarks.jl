# Chunk 2 benchmark: scalar get_rates lookup only (explicit-table, 5 resolved H states). BenchmarkTools only, after the
# fixture-loading preamble; no @time/@elapsed.
#   julia --project=<test_env or benchmark env with CosmoRec> benchmark/chunk2_benchmarks.jl
using BenchmarkTools, CosmoRec
include(joinpath(@__DIR__, "..", "test", "chunk2_helpers.jl"))
const T, _, _ = read_table_window2()
get_rates(T, 3000.0, 2700.0); get_rates(T, 7000.0, 6300.0)       # warm-up outside the measurement
for (lab, Tg, Te) in (("interior (table branch)", 3000.0, 2700.0), ("detailed-balance branch", 7000.0, 6300.0))
    A = zeros(5); B = zeros(5); R = zeros(5, 5)
    println(lab, ": get_rates (allocating)")
    display(@benchmark get_rates($T, $Tg, $Te) samples = 2000 evals = 10)
    println("\n", lab, ": get_rates! (in place)")
    display(@benchmark get_rates!($A, $B, $R, $T, $Tg, $Te) samples = 2000 evals = 10)
    println()
end
