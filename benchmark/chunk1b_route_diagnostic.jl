# Fresh-process, process-local route diagnostic for Chunk 1b. Run:
#   julia --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl/benchmark \
#       /home/marcobonici/Desktop/work/CosmologicalEmulators/CosmoRec.jl/benchmark/chunk1b_route_diagnostic.jl MODE
# MODE = count : count MooncakeVJP config hook calls and record the (eltype,size,type) of every matrix for which the
#                SparspakFactorization cache is built (forward ForwardDiff run, Mooncake forward+adjoint run);
#        throw : MooncakeVJP hook throws a sentinel (proves the adjoint path reaches it, stack printed);
# The overrides re-implement the bodies of (resolved, unedited) SciMLSensitivityMooncakeExt.get_paramjac_config and
# LinearSolveSparspakExt.init_cacheval and exist only in this process (never in Pkg.test).
using LinearAlgebra, SparseArrays, Random, ForwardDiff, Mooncake
using SparseArrays: AbstractSparseMatrixCSC
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using LinearSolve, Sparspak
using SciMLSensitivity
using SciMLSensitivity: GaussAdjoint, QuadratureAdjoint, MooncakeVJP, MooncakeLoaded
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
include(joinpath(@__DIR__, "..", "test", "chunk1b_helpers.jl"))

const MODE = isempty(ARGS) ? "count" : ARGS[1]
const HOOK = Ref(0)
const SPK = Tuple{DataType, Tuple{Int, Int}}[]
struct RouteSentinel <: Exception end

@eval function SciMLSensitivity.get_paramjac_config(::MooncakeLoaded, ::MooncakeVJP, pf, p, f, y, _t)
    HOOK[] += 1
    $(MODE == "throw") && throw(RouteSentinel())
    dy_mem = zero(y); λ_mem = zero(y)
    cache = Mooncake.prepare_pullback_cache(pf, dy_mem, y, p, _t)
    p_grad_buf = p isa AbstractArray && !(p isa Array) ? similar(p) : nothing
    return cache, pf, λ_mem, dy_mem, p_grad_buf
end

const EXT = Base.get_extension(LinearSolve, :LinearSolveSparspakExt)
function LinearSolve.init_cacheval(::SparspakFactorization, A::SparseMatrixCSC{Float64, Int}, b, u, Pl, Pr,
        maxiters::Int, abstol, reltol, verbose::Union{LinearSolve.LinearVerbosity, Bool}, assumptions::LinearSolve.OperatorAssumptions)
    push!(SPK, (typeof(A), size(A)))
    return EXT.PREALLOCATED_SPARSEPAK
end
function LinearSolve.init_cacheval(::SparspakFactorization, A::AbstractSparseMatrixCSC{Tv, Ti}, b, u, Pl, Pr,
        maxiters::Int, abstol, reltol, verbose::Union{LinearSolve.LinearVerbosity, Bool}, assumptions::LinearSolve.OperatorAssumptions) where {Tv, Ti}
    push!(SPK, (typeof(A), size(A)))
    return Sparspak.SparseCSCInterface.sparspaklu(SparseMatrixCSC{Tv, Ti}(size(A)..., SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)), factorize = false)
end

ops = MeshOps(mesh_a())
alg = Rodas5P(autodiff = AutoFiniteDiff(), linsolve = SparspakFactorization())
w = randn(MersenneTwister(1), 3 * (ops.N - 1))
backend = AutoMooncake(; config = nothing)
println("MODE=$MODE Julia $VERSION; interior system size $(ops.N - 1)")
empty!(SPK)
ForwardDiff.jacobian(t -> observe(t, ops, alg), THETA1B)
println("ForwardDiff solve: Sparspak caches built for: ", unique(SPK))
for (label, sa) in (("GaussAdjoint(MooncakeVJP)", GaussAdjoint(autojacvec = MooncakeVJP())),
                    ("QuadratureAdjoint(MooncakeVJP)", QuadratureAdjoint(autojacvec = MooncakeVJP())))
    empty!(SPK); HOOK[] = 0
    res = try
        obj = th -> dot(w, observe(th, ops, alg; sensealg = sa))
        prep = prepare_gradient(obj, backend, THETA1B)
        c1 = HOOK[]
        gradient(obj, prep, backend, THETA1B)
        "ok; hook calls after prepare=$c1, after gradient=$(HOOK[]); Sparspak caches: $(unique(SPK)) (count $(length(SPK)))"
    catch e
        frames = unique(string.(getfield.(stacktrace(catch_backtrace()), :func)))
        keep = filter(f -> occursin(r"adjoint|Adjoint|paramjac|concrete", f), frames)
        "EXCEPTION $(typeof(e)); hook calls=$(HOOK[]); frames: " * join(first(keep, 12), " <- ")
    end
    println("  $label => $res")
end
