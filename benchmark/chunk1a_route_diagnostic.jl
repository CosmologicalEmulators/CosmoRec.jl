# Fresh-process, process-local route diagnostic. Run:
#   julia --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl/benchmark \
#       /home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl/benchmark/chunk1a_route_diagnostic.jl MODE
# MODE = count : instrument MooncakeVJP config hook with a counter (behaviour preserved)
#        throw : the hook throws a sentinel; reaching it proves the hook runs
#        none  : no instrumentation (baseline)
# The override re-implements the 6-line body of
# SciMLSensitivityMooncakeExt.get_paramjac_config(::MooncakeLoaded, ::MooncakeVJP, ...)
# (resolved source quoted in docs/CHUNK1A_RESULTS.md); it exists only in this process.
using LinearAlgebra: dot
using Random, ForwardDiff, Mooncake
using SciMLBase: ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using SciMLSensitivity
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP, MooncakeLoaded
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1a_helpers.jl"))

const MODE = isempty(ARGS) ? "count" : ARGS[1]
const HOOK_CALLS = Ref(0)
struct RouteSentinel <: Exception end

if MODE != "none"
    @eval function SciMLSensitivity.get_paramjac_config(::MooncakeLoaded, ::MooncakeVJP, pf, p, f, y, _t)
        HOOK_CALLS[] += 1
        $(MODE == "throw") && throw(RouteSentinel())
        dy_mem = zero(y)
        λ_mem = zero(y)
        cache = Mooncake.prepare_pullback_cache(pf, dy_mem, y, p, _t)
        p_grad_buf = p isa AbstractArray && !(p isa Array) ? similar(p) : nothing
        return cache, pf, λ_mem, dy_mem, p_grad_buf
    end
end

alg = Rodas5P(autodiff = AutoFiniteDiff())
backend = AutoMooncake(; config = nothing)
w = randn(MersenneTwister(20260930), 6)
println("MODE=$MODE  Julia $VERSION")
for (label, sa) in (("GaussAdjoint(MooncakeVJP)", GaussAdjoint(autojacvec = MooncakeVJP())),
                    ("QuadratureAdjoint(MooncakeVJP)", QuadratureAdjoint(autojacvec = MooncakeVJP())),
                    ("sensealg=nothing", nothing))
    HOOK_CALLS[] = 0
    obj = theta -> dot(w, numeric_g(theta, alg; sensealg = sa))
    res = try
        prep = prepare_gradient(obj, backend, THETA_NOMINAL)
        c_prep = HOOK_CALLS[]
        g = gradient(obj, prep, backend, THETA_NOMINAL)
        c_eval1 = HOOK_CALLS[]
        gradient(obj, prep, backend, [1200.0, 0.8, 40.0, 1.3, 0.7])
        ref = ForwardDiff.jacobian(t -> numeric_g(t, alg), THETA_NOMINAL)' * w
        "ok: hook calls after prepare=$c_prep, after 1st gradient=$c_eval1, after 2nd=$(HOOK_CALLS[]); maxrel vs FD-jac=$(maximum(abs.(g .- ref)) / maximum(abs.(ref)))"
    catch e
        frames = unique(string.(getfield.(stacktrace(catch_backtrace()), :func)))
        keep = filter(f -> occursin(r"adjoint|Adjoint|sensitivity|SensitivityFunction|paramjac|solve_up|concrete", f), frames)
        "EXCEPTION $(typeof(e)); hook calls=$(HOOK_CALLS[]); relevant frames: " * join(first(keep, 14), " <- ")
    end
    println("  $label => $res")
end
