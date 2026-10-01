# Chunk 1b shared definitions: nonuniform-mesh method-of-lines radiation-SHAPED toy.
# NOT the CosmoRec PDE. Requires SparseArrays, LinearAlgebra, SciMLBase (ODEFunction/ODEProblem/solve) loaded by the includer.
#
# PDE (x in [0,1], nonuniform mesh x_0..x_N, s in [0, S1B_END]):
#   u_s = kappa(s) u_xx + v(s) u_x - gamma(s) u + src(s,x)
#   kappa = k0 (1 + 0.1 sin s),  v = v0 (1 + 0.5 s),  gamma = g0 (1 + 0.2 cos s)
#   theta = [k0, v0, g0, A, B, c]
# Manufactured solution (degree <= 3 in x, so 5-point stencils are exact):
#   u = A exp(-r s) (1 + x + x^2) + B exp(-s) x^3,   r = c + 0.3 k0 + 0.2 v0 + 0.1 g0
#   Dirichlet data g_L(s)=u(s,0)=A e^{-rs},  g_R(s)=u(s,1)=3A e^{-rs}+B e^{-s}   (theta-dependent)
#   u(0,x) = A (1+x+x^2) + B x^3                                                 (theta-dependent)
#   src = u_s - kappa u_xx - v u_x + gamma u, from CONTINUUM derivatives written out by hand
#   below (never from the assembled matrices).

const S1B_END = 1.0
const S1B_SAVE = [0.05, 0.3, 1.0]            # early, solver-interpolated interior, late
const THETA1B = [2.0, 3.0, 1.0, 1.0, 0.5, 0.4]

# ---------------- fixed mesh data ----------------
mesh_a() = [0.0, 0.07, 0.17, 0.3, 0.45, 0.62, 0.8, 0.9, 1.0]                                 # 9 points
mesh_b() = collect(range(0, 1; length = 12)) .^ 1.4 |> x -> (x .- x[1]) ./ (x[end] - x[1])   # 12 points, power-law clustered

# Fornberg (1988) finite-difference weights on nodes `xs` at point z, derivatives 0..m: W[node, deriv+1].
function fornberg(z, xs, m)
    n = length(xs)
    c = zeros(n, m + 1)
    c1 = 1.0
    c4 = xs[1] - z
    c[1, 1] = 1.0
    for i in 2:n
        mn = min(i - 1, m)
        c2 = 1.0
        c5 = c4
        c4 = xs[i] - z
        for j in 1:(i - 1)
            c3 = xs[i] - xs[j]
            c2 *= c3
            if j == i - 1
                for k in mn:-1:1
                    c[i, k + 1] = c1 * (k * c[i - 1, k] - c5 * c[i - 1, k + 1]) / c2
                end
                c[i, 1] = -c1 * c5 * c[i - 1, 1] / c2
            end
            for k in mn:-1:1
                c[j, k + 1] = (c4 * c[j, k + 1] - k * c[j, k]) / c3
            end
            c[j, 1] = c4 * c[j, 1] / c3
        end
        c1 = c2
    end
    return c
end

# Five-point stencil for node i (0-based index i in 0..N): centered when possible, one-sided at the ends.
function stencil_nodes(i, N)
    lo = clamp(i - 2, 0, N - 4)
    return lo:(lo + 4)
end

struct MeshOps
    x::Vector{Float64}
    N::Int
    D1::SparseMatrixCSC{Float64, Int}   # full (N+1)x(N+1) first-derivative rows (boundary rows zero: not needed)
    D2::SparseMatrixCSC{Float64, Int}
    A1::SparseMatrixCSC{Float64, Int}   # interior-interior blocks, (N-1)x(N-1)
    A2::SparseMatrixCSC{Float64, Int}
    L1L::Vector{Float64}; L1R::Vector{Float64}   # boundary lifting columns (interior rows)
    L2L::Vector{Float64}; L2R::Vector{Float64}
    pattern::SparseMatrixCSC{Float64, Int}       # union pattern of A1, A2, I (Jacobian prototype)
    nz1::Vector{Float64}; nz2::Vector{Float64}; nzI::Vector{Float64}  # A1, A2, I values aligned to pattern.nzval
end

function MeshOps(x::Vector{Float64})
    N = length(x) - 1
    @assert N >= 6 && issorted(x) && x[1] == 0.0 && x[end] == 1.0
    I1, J1, V1, V2 = Int[], Int[], Float64[], Float64[]
    for i in 1:(N - 1)
        nodes = collect(stencil_nodes(i, N))
        W = fornberg(x[i + 1], x[nodes .+ 1], 2)
        for (k, nd) in enumerate(nodes)
            push!(I1, i + 1); push!(J1, nd + 1); push!(V1, W[k, 2]); push!(V2, W[k, 3])
        end
    end
    D1 = sparse(I1, J1, V1, N + 1, N + 1)
    D2 = sparse(I1, J1, V2, N + 1, N + 1)
    int = 2:N
    A1 = D1[int, int]; A2 = D2[int, int]
    pattern = sparse(abs.(A1) .+ abs.(A2) .+ sparse(1.0I, N - 1, N - 1)); fill!(pattern.nzval, 1.0)
    function aligned(M)
        out = zeros(nnz(pattern))
        rv, cp = rowvals(pattern), pattern.colptr
        for j in 1:(N - 1), q in cp[j]:(cp[j + 1] - 1)
            out[q] = M[rv[q], j]
        end
        out
    end
    return MeshOps(x, N, D1, D2, A1, A2, Vector(D1[int, 1]), Vector(D1[int, N + 1]), Vector(D2[int, 1]),
        Vector(D2[int, N + 1]), pattern, aligned(A1), aligned(A2), aligned(sparse(1.0I, N - 1, N - 1)))
end

# ---------------- runtime (differentiable) coefficient / boundary / source maps ----------------
rate_r(th) = th[6] + 0.3 * th[1] + 0.2 * th[2] + 0.1 * th[3]
kappa(s, th) = th[1] * (1 + 0.1 * sin(s))
vel(s, th) = th[2] * (1 + 0.5 * s)
gam(s, th) = th[3] * (1 + 0.2 * cos(s))
bc_left(s, th) = th[4] * exp(-rate_r(th) * s)
bc_right(s, th) = 3 * th[4] * exp(-rate_r(th) * s) + th[5] * exp(-s)
init_u(x, th) = th[4] * (1 + x + x^2) + th[5] * x^3

u_exact(s, x, th) = th[4] * exp(-rate_r(th) * s) * (1 + x + x^2) + th[5] * exp(-s) * x^3
ux_exact(s, x, th) = th[4] * exp(-rate_r(th) * s) * (1 + 2x) + th[5] * exp(-s) * 3x^2
uxx_exact(s, x, th) = th[4] * exp(-rate_r(th) * s) * 2 + th[5] * exp(-s) * 6x
us_exact(s, x, th) = -rate_r(th) * th[4] * exp(-rate_r(th) * s) * (1 + x + x^2) - th[5] * exp(-s) * x^3
src(s, x, th) = us_exact(s, x, th) - kappa(s, th) * uxx_exact(s, x, th) - vel(s, th) * ux_exact(s, x, th) +
                gam(s, th) * u_exact(s, x, th)

# ---------------- operator application (explicit arguments, no hidden state) ----------------
# operator part: k*(A2 u + L2L gL + L2R gR) + v*(A1 u + L1L gL + L1R gR) - g*u   (explicit arguments only)
function apply_ops(ops::MeshOps, u::AbstractVector, gL, gR, k, v, g)
    return k .* (ops.A2 * u .+ ops.L2L .* gL .+ ops.L2R .* gR) .+
           v .* (ops.A1 * u .+ ops.L1L .* gL .+ ops.L1R .* gR) .- g .* u
end

function rhs_interior(ops::MeshOps, u::AbstractVector, s, th)
    xi = @view ops.x[2:(end - 1)]
    return apply_ops(ops, u, bc_left(s, th), bc_right(s, th), kappa(s, th), vel(s, th), gam(s, th)) .+
           src.(s, xi, Ref(th))
end

struct RadRHS
    ops::MeshOps
end
(f::RadRHS)(du, u, th, s) = (du .= rhs_interior(f.ops, u, s, th); nothing)

# Analytic Jacobian J = kappa A2 + v A1 - gamma I on the fixed union pattern (Dual-compatible).
struct RadJac
    ops::MeshOps
end
function (j::RadJac)(J, u, th, s)
    o = j.ops
    k = kappa(s, th); v = vel(s, th); g = gam(s, th)
    J.nzval .= k .* o.nz2 .+ v .* o.nz1 .- g .* o.nzI
    return nothing
end

u0_interior(ops::MeshOps, th) = init_u.(ops.x[2:(end - 1)], Ref(th))

# Build problem + solve from theta INSIDE the differentiated objective. `kind` selects the Jacobian/linear-solve route.
function solve_rad(th::AbstractVector{T}, ops::MeshOps, alg; sensealg = nothing, abstol = 1e-10, reltol = 1e-9,
        jac = :analytic, saveat = S1B_SAVE) where {T}
    f = if jac === :analytic
        ODEFunction(RadRHS(ops); jac = RadJac(ops), jac_prototype = convert(SparseMatrixCSC{T, Int}, ops.pattern))
    elseif jac === :proto
        ODEFunction(RadRHS(ops); jac_prototype = convert(SparseMatrixCSC{T, Int}, ops.pattern))
    else
        ODEFunction(RadRHS(ops))
    end
    prob = ODEProblem(f, u0_interior(ops, th), (0.0, S1B_END), copy(th))
    kw = (; saveat = saveat, abstol = abstol, reltol = reltol)
    return sensealg === nothing ? solve(prob, alg; kw...) : solve(prob, alg; kw..., sensealg = sensealg)
end

function observe(th, ops, alg; kw...)
    sol = solve_rad(th, ops, alg; kw...)
    @assert length(sol.u) == length(S1B_SAVE)
    return vcat(sol.u...)
end

function exact_observe(th::AbstractVector{T}, ops::MeshOps) where {T}
    xi = ops.x[2:(end - 1)]
    return vcat([T[u_exact(s, x, th) for x in xi] for s in S1B_SAVE]...)
end

replace_at(v, j, val) = (w = copy(v); w[j] = val; w)
