# Reduced reproducer (retained blocker): Rodas5P(AutoFiniteDiff) + ODEFunction(rhs!; jac=..., tgrad=...) on the
# nonlinear 2-state problem makes `prepare_gradient` (Mooncake, GaussAdjoint(MooncakeVJP)) crash the Julia 1.12.6
# compiler with SIGSEGV (decay_derived in cgutils.cpp via Mooncake _ssa_to_ids). Usage:
#   julia --project=<repo>/benchmark benchmark/chunk1c_hook_blocker.jl {j|tg|jtg}
# j = jac only (works), tg = tgrad only (works), jtg = both (segfault).
using Mooncake, LinearAlgebra
using SciMLBase: ODEFunction, ODEProblem, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using SciMLSensitivity: GaussAdjoint, MooncakeVJP
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoFiniteDiff
const SV = [0.0005, 0.02, 0.5]
rhs!(du, u, p, s) = (du[1] = -p[1] * u[1] * u[1]; du[2] = -p[2] * u[2]; nothing)
jj!(J, u, p, s) = (J[1, 1] = -2 * p[1] * u[1]; J[1, 2] = 0.0; J[2, 1] = 0.0; J[2, 2] = -p[2]; nothing)
tg!(dT, u, p, s) = (dT[1] = 0.0; dT[2] = 0.0; nothing)
const MODE = ARGS[1]
const be = AutoMooncake(; config = nothing)
const alg = Rodas5P(autodiff = AutoFiniteDiff())
const sa = GaussAdjoint(autojacvec = MooncakeVJP())
function obs(th::AbstractVector{T}) where {T}
    f = MODE == "jtg" ? ODEFunction(rhs!; jac = jj!, tgrad = tg!) : MODE == "j" ? ODEFunction(rhs!; jac = jj!) : ODEFunction(rhs!; tgrad = tg!)
    prob = ODEProblem(f, T[th[3], th[4]], (0.0, 0.5), T[th[1], th[2]])
    sol = solve(prob, alg; saveat = SV, abstol = 1.0e-12, reltol = 1.0e-11, sensealg = sa)
    return vcat(sol.u[1], sol.u[2], sol.u[3])
end
th0 = [1000.0, 1.0, 1.0, 1.0]
obj = th -> sum(obs(th))
prep = prepare_gradient(obj, be, th0)
println(MODE, " ", gradient(obj, prep, be, th0))
