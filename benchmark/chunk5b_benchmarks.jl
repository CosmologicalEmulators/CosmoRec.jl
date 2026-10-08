# Chunk 5b benchmark: one recombination pass (BenchmarkTools only; setup = table loading outside the timed region). REQUIRES COSMOREC_NATIVE_DATA_DIR.
using BenchmarkTools, CosmoRec
include(joinpath(@__DIR__, "..", "test", "chunk5_helpers.jl"))
const SOLVE = (a...) -> solve_phase5(a...; reltol = 1.0e-12, a1 = 1.0e-16, aex = 1.0e-14, alg = Rodas5P(), tstops = false)
println("one full pass, 3000 nodes, Rodas5P reltol 1e-12 (cold first call includes compilation; reported by the first sample)")
display(@benchmark recombination_pass($RM5, $SOLVE) samples = 10 evals = 1)
const P0 = recombination_pass(RM5, SOLVE)
println("\nfrozen-branch pass (k_switch given)"); display(@benchmark recombination_pass($RM5, $SOLVE; k_switch = $(P0.k_switch)) samples = 10 evals = 1)
println("\nlooser pass, reltol 1e-8"); display(@benchmark recombination_pass($RM5, (a...) -> solve_phase5(a...; reltol = 1.0e-8, a1 = 1.0e-12, aex = 1.0e-10, alg = Rodas5P(), tstops = false)) samples = 10 evals = 1)
println("\ngrid init_xarr_linear (3000)"); display(@benchmark init_xarr_linear(3000.0, 50.0, 3000) samples = 1000 evals = 1)
println()
