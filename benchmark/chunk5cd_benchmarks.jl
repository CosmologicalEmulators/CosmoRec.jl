# Chunk 5c/5d benchmark: Recfast tail, output assembly and the complete single pass (BenchmarkTools only). REQUIRES COSMOREC_NATIVE_DATA_DIR.
using BenchmarkTools, CosmoRec
include(joinpath(@__DIR__, "..", "test", "chunk5_helpers.jl"))
const SOLVE = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-16, aex = 1.0e-14)
const θ = theta5(); const Hf = hfun5()
const P0 = recombination_pass(RM5, SOLVE)
const INP = recfast_tail_inputs(RM5, P0.z[end], P0.states[end])
println("recfast_tail (200 nodes, Rodas5P 1e-12; first sample includes compilation)")
display(@benchmark recfast_tail($θ, $RFCONST5, $Hf, $(P0.z[end]), $(INP.Xe_Hi), $(INP.Xe_Hei), $(INP.Xei), $(INP.TMi), $(INP.dXei), solve_tail5) samples = 20 evals = 1)
const T0 = recfast_tail(θ, RFCONST5, Hf, P0.z[end], INP.Xe_Hi, INP.Xe_Hei, INP.Xei, INP.TMi, INP.dXei, solve_tail5)
const ROWS = recombination_output_rows(P0, T0)
const ZG = [1.0e4 - i * (1.0e4 / 9999) for i in 0:9999]
println("\nreturn_solution_to_grid (3199 rows -> 10000 grid nodes)"); display(@benchmark return_solution_to_grid($RM5, $ROWS, $ZG) samples = 50 evals = 1)
println("\ncomplete single pass (pass + tail + assembly)"); display(@benchmark recombination_history($RM5, $θ, $SOLVE, solve_tail5, $ZG) samples = 5 evals = 1)
println()
