# Chunk 3c benchmark: H-I absorption of HeI photons (explicit args, evals = 1). BenchmarkTools only; no @time/@elapsed.
#   julia --project=<test_env with CosmoRec + BenchmarkTools + ForwardDiff + Mooncake + DifferentiationInterface> benchmark/chunk3c_benchmarks.jl
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk3c_helpers.jl"))
G, lnT, lgTg, F, cases = read_fixture3c()
const fc = fcorr_spline3c(F)
const c = first(q for q in cases if q.label == "Tmid20_fe0.5_ft0.5")
const dp = window_dp_table(lnT, c)
const bt = window_bitot(lgTg, c)
const ZC = c.z
X = c.X; x0 = vcat(c.Tg, c.NH, c.Hz, c.XH1s, X)
rhs!(g, x) = hi_absorption_rhs!(g, 1, 2, 3, dp, bt, fc, ZC, x[1], x[2], x[3], x[4], @view(x[5:11]))
f(x) = (g = zeros(eltype(x), 9); rhs!(g, x); g)
g = zeros(9); rhs!(g, x0)   # warm-up
println("fcorr spline"); display(@benchmark fcorr($fc, $(c.Tg)) samples = 2000 evals = 1)
println("\nhelium_Bitot (4-point Lagrange)"); display(@benchmark helium_Bitot($bt, $(c.Tg)) samples = 2000 evals = 1)
println("\ndp_lookup singlet (4 sheets x 4x4 cubic)"); display(@benchmark dp_lookup($dp, false, $(log(c.Tg)), $(log(c.eta)), $(log(c.tauS_S * c.pdS))) samples = 2000 evals = 1)
println("\nhi_abs_singlet"); display(@benchmark hi_abs_singlet($dp, $bt, $fc, $(c.Tg), $(X[1]), $(X[3]), $(c.NH), $(c.Hz), $(c.XH1s)) samples = 2000 evals = 1)
println("\nhi_abs_triplet"); display(@benchmark hi_abs_triplet($dp, $(c.Tg), $(X[1]), $(X[6]), $(c.NH), $(c.Hz), $(c.XH1s)) samples = 2000 evals = 1)
println("\nhi_absorption_rhs! (both channels, in place)"); display(@benchmark rhs!(gg, $x0) setup = (gg = zeros(9)) samples = 2000 evals = 1)
println("\nForwardDiff.jacobian (9 outputs x 11 inputs)"); display(@benchmark ForwardDiff.jacobian($f, $x0) samples = 500 evals = 1)
obj(x, w) = dot(w, f(x)); w = collect(range(0.3, 1.7; length = 9)); backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (warmed after first sample; first sample includes compilation)")
display(@benchmark prepare_gradient($obj, $backend, $x0, Constant($w)) samples = 3 evals = 1)
prep = prepare_gradient(obj, backend, x0, Constant(w)); gradient(obj, prep, backend, x0, Constant(w))   # warm
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $x0, Constant($w)) samples = 500 evals = 1)
println()
