# Chunk 1c toy: recombination-shaped quadratic loss + slow mode. theta = [a, c, x0, y0]
#   x' = -a x^2,  y' = -c y   (a >> c > 0);  x = x0/(1+a x0 s),  y = y0 exp(-c s).
# Callers load SciMLBase (ODEProblem, ODEFunction, solve).
const S1C_SAVE = [0.0005, 0.02, 0.5]
const THETA1C = [1000.0, 1.0, 1.0, 1.0]

function nl_rhs!(du, u, p, s)
    du[1] = -p[1] * u[1]^2
    du[2] = -p[2] * u[2]
    return nothing
end

# analytic state Jacobian diag(-2 a x, -c); autonomous, so d f/ds = 0 for the ORIGINAL rhs only
function nl_jac!(J, u, p, s)
    J[1, 1] = -2 * p[1] * u[1]; J[1, 2] = zero(eltype(J))
    J[2, 1] = zero(eltype(J)); J[2, 2] = -p[2]
    return nothing
end
function nl_tgrad!(dT, u, p, s)
    dT[1] = zero(eltype(dT)); dT[2] = zero(eltype(dT))
    return nothing
end

function nl_exact(theta::AbstractVector{T}) where {T}
    a, c, x0, y0 = theta
    return T[v for s in S1C_SAVE for v in (x0 / (1 + a * x0 * s), y0 * exp(-c * s))]
end

nl_odefunction(::Val{:none}) = ODEFunction(nl_rhs!)
nl_odefunction(::Val{:jac}) = ODEFunction(nl_rhs!; jac = nl_jac!)
nl_odefunction(::Val{:jac_tgrad}) = ODEFunction(nl_rhs!; jac = nl_jac!, tgrad = nl_tgrad!)

# u0 and p built from theta inside the objective; hooks = Val(:none) | Val(:jac) | Val(:jac_tgrad) (hooks on the ORIGINAL ode only)
function nl_observe(theta::AbstractVector{T}, alg; sensealg = nothing, abstol = 1.0e-12, reltol = 1.0e-11, hooks::Val = Val(:none)) where {T}
    a, c, x0, y0 = theta
    f = nl_odefunction(hooks)
    prob = ODEProblem(f, T[x0, y0], (0.0, S1C_SAVE[end]), T[a, c])
    kw = (; saveat = S1C_SAVE, abstol = abstol, reltol = reltol)
    sol = sensealg === nothing ? solve(prob, alg; kw...) : solve(prob, alg; kw..., sensealg = sensealg)
    @assert length(sol.u) == 3
    return vcat(sol.u[1], sol.u[2], sol.u[3])
end
