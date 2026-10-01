# Chunk 3e benchmark: explicit DPesc_coh fallback (z = 1200 native fixture query), evals = 1, BenchmarkTools only. Setup/cold/hot are separate.
#   julia --project=<env with CosmoRec + BenchmarkTools + ForwardDiff + Mooncake + DifferentiationInterface> benchmark/chunk3e_benchmarks.jl
using BenchmarkTools, CosmoRec, ForwardDiff, Mooncake, LinearAlgebra
using DifferentiationInterface: prepare_gradient, gradient, Constant
using ADTypes: AutoMooncake
include(joinpath(@__DIR__, "..", "test", "chunk3e_helpers.jl"))
include(joinpath(@__DIR__, "..", "test", "chunk3d_helpers.jl"))
const Q = first(q for q in FX3E.queries if !q.triplet && q.label == "z1200" && q.code == "000")
const LINE = NATIVE_DPESC.singlet
const PD = sobolev_p(Q.pd * Q.tauS)
const ORD = dpesc_appr_I_sym(Q.T, Q.pd * Q.tauS, Q.eta, LINE, NATIVE_DPESC.hlyc, PD * 1e-5)[2]
println("Patterson levels per sub-interval (1-based into 1,3,7,...,255): ", ORD)
println("setup: generating the Patterson rules (module load, BigFloat) is paid once at precompilation")
println("\nprimal integral dpesc_appr_I_sym (adaptive native stopping rule)")
display(@benchmark dpesc_appr_I_sym($(Q.T), $(Q.pd * Q.tauS), $(Q.eta), $LINE, $(NATIVE_DPESC.hlyc), $(PD * 1e-5)) samples = 300 evals = 1)
println("\nprimal integral, pinned orders")
display(@benchmark dpesc_appr_I_sym($(Q.T), $(Q.pd * Q.tauS), $(Q.eta), $LINE, $(NATIVE_DPESC.hlyc), 0.0; orders = $ORD) samples = 300 evals = 1)
println("\nfull fallback channel dpesc_coh (singlet)")
display(@benchmark dpesc_coh($NATIVE_DPESC, false, $(Q.tauS), $(Q.eta), $(Q.T), $(Q.T), $(Q.pd)) samples = 300 evals = 1)
x0 = [Q.tauS, Q.eta, Q.T, Q.pd]
f(x) = [dpesc_coh(NATIVE_DPESC, false, x[1], x[2], x[3], x[3], x[4]; orders = ORD)[1]]
println("\nForwardDiff.jacobian (1 x 4)"); display(@benchmark ForwardDiff.jacobian($f, $x0) samples = 100 evals = 1)
obj(x, w) = dot(w, f(x)); w = [0.7]; backend = AutoMooncake(; config = nothing)
println("\nMooncake prepare_gradient (cold, includes compilation)"); display(@benchmark prepare_gradient($obj, $backend, $x0, Constant($w)) samples = 3 evals = 1)
prep = prepare_gradient(obj, backend, x0, Constant(w)); gradient(obj, prep, backend, x0, Constant(w))
println("\nMooncake hot gradient (prepared, reused)"); display(@benchmark gradient($obj, $prep, $backend, $x0, Constant($w)) samples = 100 evals = 1)
println()
