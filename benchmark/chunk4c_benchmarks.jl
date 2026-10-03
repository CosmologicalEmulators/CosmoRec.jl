# Chunk 4c benchmark: Recfast++ history (BenchmarkTools only; evals = 1). Setup, grid, primal solve, ForwardDiff Jacobian, Mooncake cold prepare and hot prepared VJP are reported separately.
#   julia --project=<env with CosmoRec + BenchmarkTools + ForwardDiff + Mooncake + DifferentiationInterface + SciMLSensitivity + OrdinaryDiffEqRosenbrock> benchmark/chunk4c_benchmarks.jl
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk4c_helpers.jl"))
const C = NATIVE_RECFAST_CONSTANTS; const Hf = hfun4c(); const θ0 = theta4c()
println("RHS (one evaluation)"); u = [0.0, 0.5, 3000.0]; du = zeros(3); p = (θ0, C, Hf)
display(@benchmark recfast_rhs3!($du, $u, $p, 1500.0) samples = 5000 evals = 1)
println("\nnode grid (recfast_grid, 6000 nodes, Saha segments)"); display(@benchmark recfast_grid($θ0, $C, $Hf) samples = 50 evals = 1)
const GRID = recfast_grid(θ0, C, Hf)
println("\nfull 6000-node history, Rodas5P reltol 1e-10, dense output (cold first sample included)"); display(@benchmark recfast_history($θ0, $C, $Hf, $GRID, solve_rodas_dense) samples = 5 evals = 1)
println("\nfull 6000-node history, Rodas5P reltol 1e-10, saveat nodes"); display(@benchmark recfast_history($θ0, $C, $Hf, $GRID, solve_rodas) samples = 5 evals = 1)
zs = [GRID.z[argmin(abs.(GRID.z .- zt))] for zt in [2500.0, 1800.0, 1200.0, 1000.0, 800.0, 500.0, 200.0, 50.0, 5.0]]
const GSEL = RecfastGrid(vcat(GRID.z[1:GRID.jode], zs), GRID.seg_end, GRID.jode)
f(θ) = (h = recfast_history(θ, C, Hf, GSEL, solve_rodas); j = GSEL.jode; vcat(h.Xe_He[j + 1:end], h.Xe_H[j + 1:end], h.TM[j + 1:end]))
println("\nselected-node history (9 ODE nodes), primal"); display(@benchmark f($θ0) samples = 20 evals = 1)
println("\nForwardDiff.jacobian (27 x 6) through the solve"); display(@benchmark ForwardDiff.jacobian($f, $θ0) samples = 5 evals = 1)
function frev(θ)
    h = recfast_history(θ, C, Hf, GSEL, (r, u0, zs_, zn, pp) -> solve_rodas(ODEFUN4C, u0, zs_, zn, pp[1]; alg = ALG_REV4C, sensealg = GaussAdjoint(autojacvec = MooncakeVJP()), reltol = 1.0e-10, abstol = [1.0e-14, 1.0e-14, 1.0e-10]))
    j = GSEL.jode; vcat(h.Xe_He[j + 1:end], h.Xe_H[j + 1:end], h.TM[j + 1:end])
end
w = collect(range(0.5, 1.5; length = 27)); obj(θ, w) = dot(w, frev(θ)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (cold, includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $θ0, Constant($w)) samples = 2 evals = 1)
prep = prepare_gradient(obj, backend, θ0, Constant(w)); gradient(obj, prep, backend, θ0, Constant(w))
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $θ0, Constant($w)) samples = 5 evals = 1)
println()
