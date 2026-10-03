# Chunk 6a benchmark: HI PDE coefficients (BenchmarkTools only). REQUIRES COSMOREC_NATIVE_DATA_DIR.
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk6_helpers.jl"))
const ROWS = NODES5[:, 1:9]
println("population splines (3000 rows)"); display(@benchmark HIPopulationSplines($ROWS; zs = ZS6, ze = ZE6) samples = 50 evals = 1)
const POPS = HIPopulationSplines(ROWS; zs = ZS6, ze = ZE6)
println("\nget_rates_all"); display(@benchmark get_rates_all($HTAB6, $LNBITOT6, 3000.0, 2999.0) samples = 2000 evals = 1)
println("\ncoefficient splines (pd, Dnem; grid of stored rows)"); display(@benchmark hi_pde_coefficients($ROWS, $POPS, $ACC5, $HTAB6, $LNBITOT6, $NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6) samples = 10 evals = 1)
const CS = hi_pde_coefficients(ROWS, POPS, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
println("\nhi_Dnem evaluation"); display(@benchmark hi_Dnem($CS, 1300.0, 2) samples = 5000 evals = 1)
X0 = vec(ROWS[:, 2:9])
function f(x)
    rows = hcat(ROWS[:, 1], reshape(x, size(ROWS, 1), 8)); pops = HIPopulationSplines(rows; zs = ZS6, ze = ZE6)
    cs = hi_pde_coefficients(rows, pops, ACC5, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6)
    return hi_Dnem(cs, 1300.0, 2) + log(hi_pd(cs, 1300.0, 2))
end
println("\nForwardDiff directional derivative (24000 inputs, one direction)"); v = copy(X0); display(@benchmark ForwardDiff.derivative(t -> f($X0 .+ t .* $v), 0.0) samples = 5 evals = 1)
obj(x, w) = w[1] * f(x); backend = AutoMooncake(; config = nothing); w = [1.0]
println("\nMooncake prepare_gradient (cold, first sample includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $X0, Constant($w)) samples = 1 evals = 1)
prep = prepare_gradient(obj, backend, X0, Constant(w)); gradient(obj, prep, backend, X0, Constant(w))
println("\nMooncake hot gradient (prepared)"); display(@benchmark gradient($obj, $prep, $backend, $X0, Constant($w)) samples = 5 evals = 1)
println()
