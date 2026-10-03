# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 7c: the HI radiation-PDE march (theta scheme, 2nd order 5-point Lagrange stencils, banded elimination) from zs to ze.
# Sources (ORIGINAL CosmoRec v3.0b, read-only): Development/ODE_PDE_Solver/PDE_solver.cpp:49-622 (init_PDE_Stepper_Data, dli_dx_xj, d2li_d2x_xj,
# setup_Lagrange_interpolation_coefficients_O2, compute_lambda_i_02, Step_PDE_O2t), Development/Simple_routines/routines.cpp (locate_JC, polint_routines,
# polint_JC), PDE_Problem/Solve_PDEs.cpp:786-975 (stepping loop: theta = 0.999 for zin > zs-100 else 0.55, dz_out = 10, lower boundary from the redshifted
# spectrum, upper boundary 0 for boundary_up = 0).
# ----------------------------------------------------------------------------------------------------------------------------------------------------

function _dli_dx_xj(i::Int, j::Int, xa, o::Int, np::Int)      # 0-based i, j on xa[o+1 .. o+np]
    x(k) = xa[o + k + 1]
    dl = 0.0
    if i == j
        for k in 0:(i - 1); dl += 1.0 / (x(i) - x(k)); end
        for k in (i + 1):(np - 1); dl += 1.0 / (x(i) - x(k)); end
    elseif i < j
        dl = 1.0 / (x(i) - x(j))
        for k in 0:(i - 1); dl *= (x(j) - x(k)) / (x(i) - x(k)); end
        for k in (i + 1):(j - 1); dl *= (x(j) - x(k)) / (x(i) - x(k)); end
        for k in (j + 1):(np - 1); dl *= (x(j) - x(k)) / (x(i) - x(k)); end
    else
        dl = 1.0 / (x(i) - x(j))
        for k in 0:(j - 1); dl *= (x(j) - x(k)) / (x(i) - x(k)); end
        for k in (j + 1):(i - 1); dl *= (x(j) - x(k)) / (x(i) - x(k)); end
        for k in (i + 1):(np - 1); dl *= (x(j) - x(k)) / (x(i) - x(k)); end
    end
    return dl
end

function _d2li_d2x_xj(i::Int, j::Int, xa, o::Int, np::Int)
    x(k) = xa[o + k + 1]
    d2 = 0.0
    if i == j
        for k in 0:(i - 1)
            dum = 0.0
            for m in 0:(k - 1); dum += 1.0 / (x(i) - x(m)); end
            for m in (k + 1):(i - 1); dum += 1.0 / (x(i) - x(m)); end
            for m in (i + 1):(np - 1); dum += 1.0 / (x(i) - x(m)); end
            d2 += 1.0 / (x(i) - x(k)) * dum
        end
        for k in (i + 1):(np - 1)
            dum = 0.0
            for m in 0:(i - 1); dum += 1.0 / (x(i) - x(m)); end
            for m in (i + 1):(k - 1); dum += 1.0 / (x(i) - x(m)); end
            for m in (k + 1):(np - 1); dum += 1.0 / (x(i) - x(m)); end
            d2 += 1.0 / (x(i) - x(k)) * dum
        end
    else
        for m in 0:(np - 1)
            if m != i && m != j
                dum = 1.0
                for k in 0:(np - 1)
                    (k != i && k != j && k != m) && (dum *= (x(j) - x(k)) / (x(i) - x(k)))
                end
                d2 += dum / (x(i) - x(m))
            end
        end
        d2 *= 2.0 / (x(i) - x(j))
    end
    return d2
end

"""Native `setup_Lagrange_interpolation_coefficients_O2`: 5-point first/second-derivative weights `dli[i, j]`, `d2li[i, j]` (row i = native point i-1;
first interior point uses the stencil x[0..4] (centre 1), the last x[n-5..n-1] (centre 3), all others are centred)."""
struct LagrangeO2
    dli::Matrix{Float64}
    d2li::Matrix{Float64}
end
function lagrange_o2(x::AbstractVector{Float64})
    np = length(x)
    dli = zeros(np, 5); d2li = zeros(np, 5)
    for j in 0:4
        dli[2, j + 1] = _dli_dx_xj(j, 1, x, 0, 5); d2li[2, j + 1] = _d2li_d2x_xj(j, 1, x, 0, 5)
        i = np - 2
        dli[i + 1, j + 1] = _dli_dx_xj(j, 3, x, i - 3, 5); d2li[i + 1, j + 1] = _d2li_d2x_xj(j, 3, x, i - 3, 5)
    end
    for i in 2:(np - 3), j in 0:4
        dli[i + 1, j + 1] = _dli_dx_xj(j, 2, x, i - 2, 5); d2li[i + 1, j + 1] = _d2li_d2x_xj(j, 2, x, i - 2, 5)
    end
    return LagrangeO2(dli, d2li)
end

"""Native `PDE_Stepper_Data` for the HI PDE: Lagrange weights, the coefficient states at the new (`cur`) and previous (`prev`) redshift, the banded
work arrays, and the `PDE_solver_initial_call_of_solver` flag. After a step `prev` holds the coefficients and side products at the new redshift (the
native pointer swap; these feed the correction integrals)."""
mutable struct HIPDEStepper{T}
    LG::LagrangeO2
    cur::HIPDEState{T}
    prev::HIPDEState{T}
    lambda::Matrix{T}
    lambdap::Matrix{T}
    Ui::Vector{T}; Vi::Vector{T}; bi::Vector{T}
    initialized::Bool
end
function HIPDEStepper{T}(x::AbstractVector{Float64}) where {T}
    np = length(x)
    return HIPDEStepper{T}(lagrange_o2(x), HIPDEState{T}(np), HIPDEState{T}(np), zeros(T, np, 5), zeros(T, np, 5), zeros(T, np), zeros(T, np), zeros(T, np), false)
end

@inline function _lambda_02!(L, row, dz, i, center, st::HIPDEState, LG::LagrangeO2)
    @inbounds for j in 1:5
        L[row, j] = LG.d2li[i, j] * st.A[i] + LG.dli[i, j] * st.B[i]
    end
    L[row, center + 1] += st.C[i]
    L[row, center + 1] -= 1.0 / (dz + 1.0e-100)
    return L
end

"""
    hi_pde_step!(y, S, m, theta, zs, ze, y_low, y_up)

Native `Step_PDE_O2t(theta, zs, ze, x, y, y_low, y_up, PDE_D, def_PDE_Lyn_and_2s1s)` (PDE_solver.cpp:386-622): one theta-scheme step of
`du/dz = A u'' + B u' + C u + D` with Dirichlet boundaries, solved by the native banded elimination (operation order kept).
"""
function hi_pde_step!(y::AbstractVector, S::HIPDEStepper, m::HIPDEModel, theta, zs, ze, y_low, y_up)
    kap = theta - 1.0; eta = kap / theta
    dz = ze - zs
    n = length(y)
    if !S.initialized
        hi_pde_rhs_coefficients!(S.prev, m, zs)
        S.initialized = true
    end
    hi_pde_rhs_coefficients!(S.cur, m, ze)
    L = S.lambda; Lp = S.lambdap; b = S.bi; U = S.Ui; V = S.Vi; LG = S.LG
    # native point i (0-based) -> Julia row i+1
    _lambda_02!(L, 2, dz * theta, 2, 1, S.cur, LG); _lambda_02!(Lp, 2, dz * kap, 2, 1, S.prev, LG)
    for i in 2:(n - 3)
        _lambda_02!(L, i + 1, dz * theta, i + 1, 2, S.cur, LG); _lambda_02!(Lp, i + 1, dz * kap, i + 1, 2, S.prev, LG)
    end
    _lambda_02!(L, n - 1, dz * theta, n - 1, 3, S.cur, LG); _lambda_02!(Lp, n - 1, dz * kap, n - 1, 3, S.prev, LG)
    l(i, j) = L[i + 1, j + 1]; lp(i, j) = Lp[i + 1, j + 1]
    @inbounds for i in 1:(n - 2)
        b[i + 1] = eta * S.prev.D[i + 1] - S.cur.D[i + 1]
    end
    i = 1
    for mm in 0:4; b[i + 1] += lp(i, mm) * eta * y[i + mm]; end
    b[i + 1] -= l(i, 0) * y_low
    i = 2
    b[i + 1] -= l(i, 0) * y_low
    @inbounds for i in 2:(n - 3)
        for mm in 0:4; b[i + 1] += lp(i, mm) * eta * y[i + mm - 1]; end
    end
    i = n - 3
    b[i + 1] -= l(i, 4) * y_up
    i = n - 2
    for mm in 0:4; b[i + 1] += lp(i, mm) * eta * y[i + mm - 2]; end
    b[i + 1] -= l(i, 4) * y_up
    # forward elimination (Ui*, Vi*, bi*)
    i = 1
    b[i + 1] /= l(i, 1)
    U[i + 1] = l(i, 2) / l(i, 1)
    V[i + 1] = l(i, 3) / l(i, 1)
    Zlow = l(i, 4) / l(i, 1)
    i = 2
    dum2 = l(i, 2) - l(i, 1) * U[i]
    b[i + 1] -= l(i, 1) * b[i]
    b[i + 1] /= dum2
    U[i + 1] = (l(i, 3) - l(i, 1) * V[i]) / dum2
    V[i + 1] = (l(i, 4) - l(i, 1) * Zlow) / dum2
    i = 3
    dum1 = l(i, 1) - l(i, 0) * U[i - 1]
    dum2 = l(i, 2) - l(i, 0) * V[i - 1]
    b[i + 1] -= l(i, 0) * b[i - 1]
    U[i + 1] = l(i, 3) - l(i, 0) * Zlow
    dum2 -= dum1 * U[i]
    b[i + 1] -= dum1 * b[i]
    b[i + 1] /= dum2
    U[i + 1] -= dum1 * V[i]
    U[i + 1] /= dum2
    V[i + 1] = l(i, 4) / dum2
    @inbounds for i in 4:(n - 4)
        dum1 = l(i, 1) - l(i, 0) * U[i - 1]
        dum2 = l(i, 2) - l(i, 0) * V[i - 1]
        b[i + 1] -= l(i, 0) * b[i - 1]
        U[i + 1] = l(i, 3)
        dum2 -= dum1 * U[i]
        b[i + 1] -= dum1 * b[i]
        b[i + 1] /= dum2
        U[i + 1] -= dum1 * V[i]
        U[i + 1] /= dum2
        V[i + 1] = l(i, 4) / dum2
    end
    i = n - 3
    dum1 = l(i, 1) - l(i, 0) * U[i - 1]
    dum2 = l(i, 2) - l(i, 0) * V[i - 1]
    b[i + 1] -= l(i, 0) * b[i - 1]
    U[i + 1] = l(i, 3)
    dum2 -= dum1 * U[i]
    b[i + 1] -= dum1 * b[i]
    b[i + 1] /= dum2
    U[i + 1] -= dum1 * V[i]
    U[i + 1] /= dum2
    V[i + 1] = 0.0
    i = n - 2
    dum1 = l(i, 1) - l(i, 0) * U[i - 2]
    dum2 = l(i, 2) - l(i, 0) * V[i - 2]
    dum3 = l(i, 3)
    b[i + 1] -= l(i, 0) * b[i - 2]
    dum2 -= dum1 * U[i - 1]
    dum3 -= dum1 * V[i - 1]
    b[i + 1] -= dum1 * b[i - 1]
    dum3 -= dum2 * U[i]
    b[i + 1] -= dum2 * b[i]
    b[i + 1] /= dum3
    U[i + 1] = 0.0; V[i + 1] = 0.0
    # boundaries and back substitution
    y[1] = y_low
    y[n] = y_up
    y[n - 1] = b[n - 1]
    @inbounds for i in (n - 3):-1:1
        y[i + 1] = b[i + 1] - U[i + 1] * y[i + 2] - V[i + 1] * y[i + 3]
    end
    y[2] -= Zlow * y[5]
    S.cur, S.prev = S.prev, S.cur
    return y
end

"""Native `locate_JC(xx, n, x, &j)` (bisection; returns 0 / n-1 at or beyond the ends); 0-based result on `xx[1:n]`."""
function _locate_jc(xx, n::Int, x)
    asc = xx[n] > xx[1]
    if asc
        x <= xx[1] && return 0
        x >= xx[n] && return n - 1
    else
        x >= xx[1] && return 0
        x <= xx[n] && return n - 1
    end
    jl = 0; ju = n
    while ju - jl > 1
        jm = (ju + jl) >> 1
        if (x >= xx[jm + 1]) == asc
            jl = jm
        else
            ju = jm
        end
    end
    return jl
end

"""Native `polint_routines` (Neville; returns `(y, dy)`) on `xa[o+1 .. o+n]`."""
function _polint_routines(xa, ya, o::Int, n::Int, x)
    ns = 0
    dif = abs(x - xa[o + 1])
    T = promote_type(eltype(ya), typeof(x))
    c = Vector{T}(undef, n); d = Vector{T}(undef, n)
    for i in 0:(n - 1)
        dift = abs(x - xa[o + i + 1])
        if dift < dif
            ns = i; dif = dift
        end
        c[i + 1] = ya[o + i + 1]; d[i + 1] = ya[o + i + 1]
    end
    y = ya[o + ns + 1]; ns -= 1
    dy = zero(T)
    for m in 1:(n - 1)
        for i in 0:(n - m - 1)
            ho = xa[o + i + 1] - x
            hp = xa[o + i + m + 1] - x
            w = c[i + 2] - d[i + 1]
            den = ho - hp
            den == 0.0 && error("polint error")
            den = w / den
            d[i + 1] = hp * den
            c[i + 1] = ho * den
        end
        if 2 * (ns + 1) < (n - m)
            dy = c[ns + 2]
        else
            dy = d[ns + 1]; ns -= 1
        end
        y += dy
    end
    return y, dy
end

"""Native `polint_JC(xa, ya, na, x, istart, npol)`: Neville interpolation of degree npol-1 around the bracketing node (0 above the last node)."""
function polint_jc(xa, ya, x, npol::Int)
    na = length(xa)
    j = _locate_jc(xa, na - 1, x)
    j == na - 1 && return (zero(eltype(ya)), zero(eltype(ya)))
    npl = min(na, npol)
    nlow = (npl - 1) ÷ 2
    ks = min(max(j - nlow, 0), na - 1 - npl)
    return _polint_routines(xa, ya, ks, npl, x)
end

"""
    hi_pde_march(m; zs = 2500.0, ze = 500.0, dz_out = 10.0, theta = 0.55, on_step = nothing)

The production HI PDE march of `compute_DPesc_with_diffusion_equation_effective` (Solve_PDEs.cpp:786-975) from `y = 0` at `zs` to `ze` with the
native step control. `on_step(step, zin, zout, y, S)` is called after every step (0-based step; `S.prev` = coefficients/side products at `zout`).
Returns `(y, zout_list)`.
"""
function hi_pde_march(m::HIPDEModel; zs = 2500.0, ze = 500.0, dz_out = 10.0, theta = 0.55, on_step = nothing, T::Type = Float64)
    x = m.setup.x
    y = zeros(T, length(x))
    S = HIPDEStepper{T}(x)
    zin = zs; zout = zs; step = 0
    zl = Float64[]
    while true
        th = zin > zs - 100.0 ? 0.999 : theta
        zout = max(ze, zin - dz_out)
        xeval = x[1] * (1.0 + zin) / (1.0 + zout)
        y_low, _ = polint_jc(x, y, xeval, 6)
        hi_pde_step!(y, S, m, th, zin, zout, y_low, zero(T))
        push!(zl, zout)
        on_step === nothing || on_step(step, zin, zout, y, S)
        zin = zout; step += 1
        zout > ze || break
    end
    return y, zl
end
