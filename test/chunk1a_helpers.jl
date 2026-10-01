# Shared Chunk 1a toy problem definition (no @testset). Included by
# test/chunk1a_stiff_ad_probe.jl, benchmark/chunk1a_benchmarks.jl and
# benchmark/chunk1a_route_diagnostic.jl. Callers load the packages.
#
# theta = [a, c, b, y10, y20]; y1' = -a*y1 + b*y2, y2' = -c*y2 (a >> c > 0).
# Time is the forward coordinate s (s = z_start - z convention of the physics port).

const S_EARLY = 0.0005
const S_INTERP = 0.02
const S_LATE = 0.5
const S_SAVE = [S_EARLY, S_INTERP, S_LATE]
const THETA_NOMINAL = [1000.0, 1.0, 50.0, 1.0, 1.0]

function stiff_rhs!(du, u, p, s)
    a, c, b = p[1], p[2], p[3]
    du[1] = -a * u[1] + b * u[2]
    du[2] = -c * u[2]
    return nothing
end

# Closed form (a != c), flat [y1,y2] at S_EARLY, S_INTERP, S_LATE. Solver-free oracle.
function analytic_g(theta::AbstractVector{T}) where {T}
    a, c, b, y10, y20 = theta
    y1(s) = y10 * exp(-a * s) + b * y20 * (exp(-c * s) - exp(-a * s)) / (a - c)
    y2(s) = y20 * exp(-c * s)
    return T[y1(S_EARLY), y2(S_EARLY), y1(S_INTERP), y2(S_INTERP), y1(S_LATE), y2(S_LATE)]
end

# u0 and p are built from theta inside the objective. saveat gives the interpolated
# interior observation at S_INTERP.
function numeric_g(theta::AbstractVector{T}, alg; sensealg = nothing,
        abstol = 1.0e-12, reltol = 1.0e-11) where {T}
    a, c, b, y10, y20 = theta
    prob = ODEProblem(stiff_rhs!, T[y10, y20], (0.0, S_LATE), T[a, c, b])
    sol = if sensealg === nothing
        solve(prob, alg; saveat = S_SAVE, abstol = abstol, reltol = reltol)
    else
        solve(prob, alg; saveat = S_SAVE, abstol = abstol, reltol = reltol, sensealg = sensealg)
    end
    @assert length(sol.u) == 3
    return vcat(sol.u[1], sol.u[2], sol.u[3])
end

# Central differences with per-parameter step h_j = hrel * max(|theta_j|, 1).
function fd_jacobian(g::Function, theta::Vector{Float64}, hrel::Float64)
    n = length(theta)
    J = zeros(length(g(theta)), n)
    for j in 1:n
        h = hrel * max(abs(theta[j]), 1.0)
        tp = copy(theta); tp[j] += h
        tm = copy(theta); tm[j] -= h
        J[:, j] = (g(tp) .- g(tm)) ./ (2h)
    end
    return J
end
