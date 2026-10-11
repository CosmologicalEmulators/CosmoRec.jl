# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Faithful port of the ORIGINAL Recfast++ ODE solver (variable-order Gear/BDF with variable step size, Newton iteration with a 4-point numerical
# Jacobian and Gaussian elimination) of CosmoRec v3.0b (SHA 086769055f61ae0c244a53dd381ee65b624d0ac3), (c) J. Chluba et al.:
#   Development/Cosmology/Recfast++/src/ODE_solver.Recfast.cpp (cited below as RF:<line>), driven as in Recfast++.cpp:284-325 (Xe_frac ODE part).
# Native names and data flow are kept: ODE_Solver_info -> RFGearSolver, ODE_solver_Solution -> RFGearSolution, Snptr[0..5] -> S.Sn[1..6] (references
# rotated exactly like the native pointers). Float64 only (no AD/precision changes). Arithmetic is written in the native operation order; `pow` uses the
# C library (_cpow, StandaloneInit.jl) as the native object code does.
# The solver integrates the native 4-state system y = (Xe_He, Xe_H, TM, t_cos) with the Recfast++ RHS (`recfast_rhs`, RecfastPP.jl); the 4th state is
# the cosmic time (RF fcn -> evaluate_Recfast_System, f[3] = -1/(H (1+z))), initialised with t_cos_rad(zs) = 1/2/H(zs) (cosmology.Recfast.cpp:92-93).
# ----------------------------------------------------------------------------------------------------------------------------------------------------

"""Native `ODE_solver_Solution` (z, y, dy)."""
mutable struct RFGearSolution
    z::Float64
    y::Vector{Float64}
    dy::Vector{Float64}
end
RFGearSolution(neq::Int) = RFGearSolution(0.0, zeros(neq), zeros(neq))

"""Native `ODE_Solver_data` (RF:44-62) plus the RHS closure `f!(f, z, y)`."""
mutable struct RFGearSolver{F}
    f!::F
    Stemp::RFGearSolution
    Snew::RFGearSolution
    Sguess::RFGearSolution
    Sn::Vector{RFGearSolution}       # Snptr[0..5]
    Jac::Vector{Float64}             # row-major Jac[J + neq*row]
    F::Vector{Float64}
    dY::Vector{Float64}
    abs_vector::Vector{Float64}
    Dz_z_last::Float64
    order::Int
    count::Int
    n_up::Int
    n_down::Int
end

_rf_f!(S::RFGearSolver, z, y, f) = (S.f!(f, z, y); f)                                      # ODE_solver_f (RF:253-260)

# ---- numerical Jacobian (RF:276-335) ----
_rf_df_dx_4point(fp2, fp1, fm1, fm2, h) = (8.0 * (fp1 - fm1) + fm2 - fp2) / (12.0 * h)
function _rf_jac!(S::RFGearSolver, neq, col, z, Dz, y, r)
    y0 = y[col]; Dyj = y[col] * 1.0e-12
    fp2 = zeros(neq); fp1 = zeros(neq); fm1 = zeros(neq); fm2 = zeros(neq)
    y[col] = y0 + 2.0 * Dyj; _rf_f!(S, z, y, fp2)
    y[col] = y0 + Dyj;       _rf_f!(S, z, y, fp1)
    y[col] = y0 - Dyj;       _rf_f!(S, z, y, fm1)
    y[col] = y0 - 2.0 * Dyj; _rf_f!(S, z, y, fm2)
    y[col] = y0
    for k in 1:neq
        r[k] = -Dz * _rf_df_dx_4point(fp2[k], fp1[k], fm1[k], fm2[k], Dyj)
    end
    r[col] += 1.0
    return r
end
function _rf_compute_Jacobian_Matrix!(S::RFGearSolver, z, Dz, neq, y, Jac)
    Jij = zeros(neq)
    for J in 1:neq
        _rf_jac!(S, neq, J, z, Dz, y, Jij)
        for row in 1:neq
            Jac[J + neq * (row - 1)] = Jij[row]
        end
    end
    return Jac
end

"""Native `ODE_Solver_set_up_solution_and_memory` (RF:353-416)."""
function RFGearSolver(f!, z::Float64, Sz::RFGearSolution)
    neq = length(Sz.y)
    neq == 4 || throw(ArgumentError("please set initial condition for system (native requires neq == 4)"))
    Sn = [RFGearSolution(neq) for _ in 1:6]
    S = RFGearSolver(f!, RFGearSolution(neq), RFGearSolution(neq), RFGearSolution(neq), Sn, zeros(neq * neq), zeros(neq), zeros(neq),
                     zeros(neq), 0.0, 1, 1, 0, 0)
    S.Sn[1].z = Sz.z; S.Sn[1].y = copy(Sz.y); S.Sn[1].dy = copy(Sz.dy)
    _rf_f!(S, z, S.Sn[1].y, S.Sn[1].dy)
    S.abs_vector[1] = 1.0e-10; S.abs_vector[2] = 1.0e-8; S.abs_vector[3] = 1.0e-10              # RF:403-405 (abs_vector[3] stays 0)
    S.Dz_z_last = 0.0; S.order = 1; S.count = 1; S.n_up = 0; S.n_down = 0
    return S
end

# ---- Gear (BDF) coefficients for variable step size (RF:419-497) ----
_gk1(rk) = 2.0 + rk
_gk2(r1, rk) = (r = _gk1(r1) * _gk1(rk); r - 1.0)
_gk3(r1, r2, rk) = (r = _gk2(r1, r2) * _gk1(rk); r - _gk1(r1 + r2))
function _gk4(r1, r2, r3, rk)
    r = _gk3(r1, r2, r3) * _gk1(rk)
    r -= _gk1(r1 + r2) * r3
    return r - _gk2(r1, r2)
end
function _gk5(r1, r2, r3, r4, rk)
    r = _gk4(r1, r2, r3, r4) * _gk1(rk)
    r -= (_gk2(r1, r2 + r3) + r2 * r3) * r4
    return r - _gk3(r1, r2, r3)
end
_rf_delta0(r1, r2, r3, r4, r5, a1, a2, a3, a4, a5) = 1.0 + a1 * r1 + a2 * r2 + a3 * r3 + a4 * r4 + a5 * r5
_rf_alp0(a1, a2, a3, a4, a5) = 1.0 - a1 - a2 - a3 - a4 - a5
_rf_alp1(r1, r2, r3, r4, r5, a2, a3, a4, a5) = -((1.0 + a2 * r2 * (2.0 + r2) + a3 * r3 * (2.0 + r3) + a4 * r4 * (2.0 + r4) + a5 * r5 * (2.0 + r5)) / (r1 * (2.0 + r1)))
function _rf_alp2(r1, r2, r3, r4, r5, a3, a4, a5)
    t = (1.0 + r1)
    return -(t * t + a3 * r3 * (r1 - r3) * _gk2(r1, r3) + a4 * r4 * (r1 - r4) * _gk2(r1, r4) + a5 * r5 * (r1 - r5) * _gk2(r1, r5)) /
           (r2 * (r1 - r2) * _gk2(r1, r2))
end
function _rf_alp3(r1, r2, r3, r4, r5, a4, a5)
    t = (1.0 + r1) * (1.0 + r2)
    return -(t * t + a4 * r4 * (r1 - r4) * (r2 - r4) * _gk3(r1, r2, r4) + a5 * r5 * (r1 - r5) * (r2 - r5) * _gk3(r1, r2, r5)) /
           (r3 * (r1 - r3) * (r2 - r3) * _gk3(r1, r2, r3))
end
function _rf_alp4(r1, r2, r3, r4, r5, a5)
    t = (1.0 + r1) * (1.0 + r2) * (1.0 + r3)
    return -(t * t + a5 * r5 * (r1 - r5) * (r2 - r5) * (r3 - r5) * _gk4(r1, r2, r3, r5)) /
           (r4 * (r1 - r4) * (r2 - r4) * (r3 - r4) * _gk4(r1, r2, r3, r4))
end
function _rf_alp5(r1, r2, r3, r4, r5)
    t = (1.0 + r1) * (1.0 + r2) * (1.0 + r3) * (1.0 + r4)
    return -t * t / (r5 * (r1 - r5) * (r2 - r5) * (r3 - r5) * (r4 - r5) * _gk5(r1, r2, r3, r4, r5))
end
function _rf_alpi(i, r, a)          # r, a 1-based copies of the native 0-based arrays
    i == 2 && return _rf_alp1(r[1], r[2], r[3], r[4], r[5], a[3], a[4], a[5], a[6])
    i == 3 && return _rf_alp2(r[1], r[2], r[3], r[4], r[5], a[4], a[5], a[6])
    i == 4 && return _rf_alp3(r[1], r[2], r[3], r[4], r[5], a[5], a[6])
    i == 5 && return _rf_alp4(r[1], r[2], r[3], r[4], r[5], a[6])
    i == 6 && return _rf_alp5(r[1], r[2], r[3], r[4], r[5])
    return 0.0
end
"""Native `ODE_Solver_compute_ynp1` (RF:513-552): Gear corrector; returns `delta`."""
function _rf_compute_ynp1!(order, Snp1_new::RFGearSolution, Snp1_guess::RFGearSolution, Sn)
    (order < 1 || order > 6) && error("check order for ODE_Solver_compute_ynp1")
    Dznp1 = Snp1_guess.z - Sn[1].z
    a = zeros(6); rho = zeros(5)
    for k in 2:order
        Dzn = Sn[1].z - Sn[k].z; rho[k - 1] = Dzn / Dznp1
    end
    for k in order:-1:2
        a[k] = _rf_alpi(k, rho, a)
    end
    a[1] = _rf_alp0(a[2], a[3], a[4], a[5], a[6])
    delta = _rf_delta0(rho[1], rho[2], rho[3], rho[4], rho[5], a[2], a[3], a[4], a[5], a[6])
    for i in eachindex(Snp1_new.y)
        Snp1_new.y[i] = Dznp1 * delta * Snp1_guess.dy[i]
        Sa = 0.0
        for mo in (order - 1):-1:0
            Sa += a[mo + 1] * Sn[mo + 1].y[i]
        end
        Snp1_new.y[i] += Sa
    end
    return delta
end
_rf_beta0(b1, b2, b3, b4, b5) = 1.0 - b1 - b2 - b3 - b4 - b5
_rf_beta1(r1, r2, r3, r4, r5, b2, b3, b4, b5) = -(1.0 + b2 * r2 + b3 * r3 + b4 * r4 + b5 * r5) / r1
_rf_beta2(r1, r2, r3, r4, r5, b3, b4, b5) = -((1.0 + r1) + b3 * r3 * (r1 - r3) + b4 * r4 * (r1 - r4) + b5 * r5 * (r1 - r5)) / ((r1 - r2) * r2)
_rf_beta3(r1, r2, r3, r4, r5, b4, b5) = ((1.0 + r1) * (1.0 + r2) + b4 * r4 * (r1 - r4) * (r2 - r4) + b5 * r5 * (r1 - r5) * (r2 - r5)) / ((r1 - r3) * r3 * (r3 - r2))
_rf_beta4(r1, r2, r3, r4, r5, b5) = -((1.0 + r1) * (1.0 + r2) * (1.0 + r3) + b5 * r5 * (r1 - r5) * (r2 - r5) * (r3 - r5)) / ((r3 - r4) * r4 * (r4 - r1) * (r4 - r2))
_rf_beta5(r1, r2, r3, r4, r5) = (1.0 + r1) * (1.0 + r2) * (1.0 + r3) * (1.0 + r4) / ((r2 - r5) * r5 * (r5 - r1) * (r5 - r3) * (r5 - r4))
function _rf_betai(i, r, b)
    i == 2 && return _rf_beta1(r[1], r[2], r[3], r[4], r[5], b[3], b[4], b[5], b[6])
    i == 3 && return _rf_beta2(r[1], r[2], r[3], r[4], r[5], b[4], b[5], b[6])
    i == 4 && return _rf_beta3(r[1], r[2], r[3], r[4], r[5], b[5], b[6])
    i == 5 && return _rf_beta4(r[1], r[2], r[3], r[4], r[5], b[6])
    i == 6 && return _rf_beta5(r[1], r[2], r[3], r[4], r[5])
    return 0.0
end
"""Native `ODE_Solver_extrapolate_ynp1` (RF:593-631): predictor from the stored solutions."""
function _rf_extrapolate_ynp1!(order, Snp1::RFGearSolution, Sn)
    (order < 1 || order > 6) && error("check order for ODE_Solver_extrapolate_ynp1")
    Dznp1 = Snp1.z - Sn[1].z
    b = zeros(6); rho = zeros(5)
    for k in 2:order
        Dzn = Sn[1].z - Sn[k].z; rho[k - 1] = Dzn / Dznp1
    end
    for k in order:-1:2
        b[k] = _rf_betai(k, rho, b)
    end
    b[1] = _rf_beta0(b[2], b[3], b[4], b[5], b[6])
    for i in eachindex(Sn[1].y)
        Sa = 0.0
        for mo in (order - 1):-1:0
            Sa += b[mo + 1] * Sn[mo + 1].y[i]
        end
        Snp1.y[i] = Sa
    end
    return Snp1
end

# ---- error control and step size (RF:653-770) ----
function _rf_error_check(rtol, atol, y, Dy)
    Dy = abs(Dy)
    err_aim = rtol * abs(y); fac = 1.0e+10
    err_aim <= atol && (err_aim = atol)
    Dy != 0.0 && (fac = err_aim / Dy)
    return fac
end
"""Native `ODE_Solver_estimate_error_and_next_possible_step_size` (RF:676-770); returns `(redo_run, Dz_z)`."""
function _rf_estimate_error_and_next_step!(S::RFGearSolver, order, zin, z, zend, Dz_z, Dz_z_max, tol, Snew::RFGearSolution, Sn)
    order = min(5, order)
    neq = length(Snew.y)
    S.Stemp.z = z
    _rf_extrapolate_ynp1!(order, S.Stemp, Sn)
    for i in 1:neq
        S.Stemp.y[i] -= Snew.y[i]
    end
    fac = 1.0e+10
    power = order <= 1 ? 1.0 : 1.0 / order
    for c in 1:3
        dum = _rf_error_check(tol, S.abs_vector[c], Snew.y[c], S.Stemp.y[c])
        fac > dum && (fac = dum)
    end
    Dz_z_min = 1.0e-12
    Dz_z_new = Dz_z * min(100.0, _cpow(fac, power))
    Dz_z_new = min(max(Dz_z_new, Dz_z_min), Dz_z_max)
    if Dz_z * 0.9 > Dz_z_new
        Dz_z = Dz_z_new
        S.n_down = min(S.n_down + 1, 5)
        if S.n_down >= 3
            if S.n_down >= 5
                Dz_z = Dz_z_min
            else
                Dz_z /= _cpow(2.0, Float64(S.n_down - 1))
            end
        end
        Dz_z * 0.5 > Dz_z_new && return (1, Dz_z)
    elseif Dz_z * 1.1 < Dz_z_new
        Dz_z = min(Dz_z_new, Dz_z_max)
        S.n_down = 0
        return (0, Dz_z)
    end
    S.n_down = 0
    return (0, Dz_z)
end

"""Native `ODE_Solver_Solve_LA_system_Gaussian` (RF:785-806; used for neq <= 10, RF:815). A is row-major, modified in place."""
function _rf_solve_gaussian!(A, b, x)
    neq = length(b)
    for cdel in 0:(neq - 1)
        p = A[cdel + neq * cdel + 1]
        for c in cdel:(neq - 1); A[c + neq * cdel + 1] /= p; end
        b[cdel + 1] /= p
        for r in (1 + cdel):(neq - 1)
            p = A[cdel + neq * r + 1]
            for c in cdel:(neq - 1); A[c + neq * r + 1] -= p * A[c + neq * cdel + 1]; end
            b[r + 1] -= p * b[cdel + 1]
        end
    end
    for r in (neq - 1):-1:0
        x[r + 1] = b[r + 1]
        for c in (r + 1):(neq - 1); x[r + 1] -= x[c + 1] * A[c + r * neq + 1]; end
    end
    return x
end

"""Native `ODE_Solver_compute_new_Solution` (RF:848-933): Newton iteration of the Gear corrector; returns LAerror (0 ok, 1 = more than 25 Jacobian loops)."""
function _rf_compute_new_Solution!(S::RFGearSolver, order, zout, reltol)
    Sn = S.Sn
    zin = Sn[1].z; Dz = zout - zin
    neq = length(Sn[1].y)
    tolSol = reltol; tolJac = reltol / 2.0
    S.Sguess.z = zout; S.Snew.z = zout
    _rf_extrapolate_ynp1!(order, S.Sguess, Sn)
    Jac_loops = 0
    converged = 0
    while converged == 0
        Jac_loops += 1
        _rf_f!(S, zout, S.Sguess.y, S.Sguess.dy)
        delta = _rf_compute_ynp1!(order, S.Snew, S.Sguess, Sn)
        for i in 1:neq; S.F[i] = S.Snew.y[i] - S.Sguess.y[i]; end
        _rf_compute_Jacobian_Matrix!(S, zout, Dz * delta, neq, S.Sguess.y, S.Jac)
        _rf_solve_gaussian!(S.Jac, S.F, S.dY)                                     # tolJac unused by the Gaussian route (RF:815)
        for i in 1:neq; S.Sguess.y[i] += S.dY[i]; end
        converged = 1
        for k in 1:neq
            if abs(S.dY[k]) >= max(tolSol * abs(S.Sguess.y[k]), S.abs_vector[k])
                converged = 0; break
            end
        end
        if converged == 0
            Jac_loops % 5 == 0 && (tolJac /= 2.0)
            Jac_loops > 25 && return 1
        end
    end
    _rf_f!(S, zout, S.Sguess.y, S.Sguess.dy)
    S.Snew, S.Sguess = S.Sguess, S.Snew                                            # native pointer swap (RF:928-932)
    return 0
end

"""Native `ODE_Solver_Solve_history(zs, zend, Sz)` (RF:940-1011): advances the solver from `zs` to `zend`; fills `Sz` (z, y, dy) from `Snptr[0]`."""
function rf_gear_solve_history!(S::RFGearSolver, zs::Float64, zend::Float64, Sz::RFGearSolution)
    zout = zs; zin = zs
    tolSol = 1.0e-8
    Dz_z = max(1.0e-8, S.Dz_z_last); Dz_z_max = 10.0 / zs
    maxorder = 4
    while true
        zout = max(zin * (1.0 - Dz_z), zend)
        LAerror = _rf_compute_new_Solution!(S, S.order, zout, tolSol)
        if LAerror == 1
            Dz_z /= 2.0
            Dz_z = max(Dz_z, 1.0e-10)
        else
            redo_run, Dz_z = _rf_estimate_error_and_next_step!(S, S.order, zin, zout, zend, Dz_z, Dz_z_max, tolSol, S.Snew, S.Sn)
            if redo_run == 0
                Sptr = S.Sn[6]
                for k in 6:-1:2; S.Sn[k] = S.Sn[k - 1]; end
                S.Sn[1] = S.Snew
                S.Snew = Sptr
                S.count += 1
                if S.count >= S.order + 1
                    S.order = min(S.order + 1, maxorder)
                    S.count = 1
                end
                zin = zout
            end
        end
        zout > zend || break
    end
    S.Dz_z_last = Dz_z
    Sz.z = S.Sn[1].z; Sz.y = copy(S.Sn[1].y); Sz.dy = copy(S.Sn[1].dy)
    return Sz
end

"""
    recfast_gear_solve(rhs!, u0, zs, znodes, p) -> (states, derivatives)

The ORIGINAL Recfast++ ODE driver of `Xe_frac` (Recfast++.cpp:287-325) as a `recfast_history` callback: native 4-state system (Xe_He, Xe_H, TM, t_cos)
with `t_cos(zs) = 1/2/H(zs)` and the native-order RHS `evaluate_recfast_system_native`, one `ODE_Solver_Solve_history` call per node, derivatives = the solver's `Sz.dy` (RHS at the accepted state).
`p = (θ, c, Hfun)` as passed by `recfast_history`; `rhs!` is not used (the native RHS includes the 4th, cosmic-time, component).
"""
function recfast_gear_solve(rhs!, u0, zs, znodes, p)
    θ, c, Hfun = p
    θ64 = Float64.(θ)
    f! = (f, z, y) -> begin
        r = evaluate_recfast_system_native(θ64, z, y, Float64(Hfun(z)), c)          # native-order RHS (RF fcn -> evaluate_Recfast_System)
        f[1] = r[1]; f[2] = r[2]; f[3] = r[3]; f[4] = r[4]
        return f
    end
    z0 = Float64(zs)
    Sz = RFGearSolution(z0, [Float64(u0[1]), Float64(u0[2]), Float64(u0[3]), t_cos_rad_native(Float64(Hfun(z0)))], zeros(4))
    S = RFGearSolver(f!, z0, Sz)
    n = length(znodes)
    states = zeros(3, n); ders = zeros(3, n)
    states[:, 1] .= Sz.y[1:3]; ders[:, 1] .= S.Sn[1].dy[1:3]
    zcur = z0
    for i in 2:n
        ze = Float64(znodes[i])
        rf_gear_solve_history!(S, zcur, ze, Sz)
        zcur = ze
        states[:, i] .= Sz.y[1:3]; ders[:, i] .= Sz.dy[1:3]
    end
    return states, ders
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Native-order Recfast++ RHS (Float64 only): `evaluate_Recfast_System` (src/evalode.Recfast.cpp:99-220, cited EV:<line>) with its helpers Boltzmann (EV:52-53),
# SahaBoltz (EV:58-63), alphaH_func (EV:68-72), betaH_func (EV:77-78), alphaHe_func (EV:83-88), betaHe_func (EV:93-94) and the Recfast++ cosmology
# accessors NH, f_He, TCMB (cosmology.Recfast.cpp:57-120). Same equations as the generic `recfast_rhs` (RecfastPP.jl, kept for AD: it uses an
# overflow-safe 1/SahaBoltz form), but in the EXACT native operation order with the C library exp/pow; inactive native branches (DM annihilation/decay,
# magnetic heating, reionization: all off in the CAMB batch route) and the variation-of-constants factors (all exactly 1.0, Variation_constants.Recfast.cpp:50)
# are omitted because they are arithmetic no-ops. GCC folds pow(x, 2.0) to x*x (exact), reproduced as x*x.
# ----------------------------------------------------------------------------------------------------------------------------------------------------
const _RF_LIBM = Sys.islinux() ? "libm.so.6" : ""
_cexp(x::Float64) = isempty(_RF_LIBM) ? exp(x) : ccall((:exp, _RF_LIBM), Float64, (Float64,), x)
_clog(x::Float64) = isempty(_RF_LIBM) ? log(x) : ccall((:log, _RF_LIBM), Float64, (Float64,), x)
# C-library pow for Float64 (used by the generic Recfast++ helpers in RecfastPP.jl); generic (AD) types keep Julia ^.
_npow(x::Float64, y::Float64) = _cpow(x, y); _npow(x, y) = x^y

_nat_Boltzmann(gj, gi, E, T, c::RecfastConstants) = max(1.0e-300, ((gj / gi) * _cexp(-E / (c.kBoltz * T))))
function _nat_SahaBoltz(gi, gc, ge, E_ion, T, c::RecfastConstants)
    c1 = (c.hPlanck / (2.0 * c.pi * c.mElect)) * (c.hPlanck / c.kBoltz)
    return (ge * gi / (2.0 * gc) * _cpow(c1, 1.5) * _cpow(T, -1.5) * _cexp(min(680.0, E_ion / (c.kBoltz * T))))
end
function _nat_alphaH(F::Float64, TM::Float64)
    a1 = 4.309; a2 = -0.6166; a3 = 0.6703; a4 = 0.5300; t = TM / 1.0e+4
    return F * a1 * 1.0e-19 * _cpow(t, a2) / (1.0 + a3 * _cpow(t, a4))
end
_nat_betaH(F, T, c) = _nat_alphaH(F, T) / _nat_SahaBoltz(2.0, 1.0, 1.0, c.EionH2s, T, c)
function _nat_alphaHe(TM::Float64)
    a1 = _cpow(10.0, -16.744); a2 = 0.711; T0 = 3.0; T1 = _cpow(10.0, 5.114)
    sqrt_TMT0 = sqrt(TM / T0); sqrt_TMT1 = sqrt(TM / T1)
    return a1 / (sqrt_TMT0 * _cpow(1.0 + sqrt_TMT0, 1.0 - a2) * _cpow(1.0 + sqrt_TMT1, 1.0 + a2))
end
_nat_betaHe(T, c) = _nat_alphaHe(T) / _nat_SahaBoltz(1.0, 2.0, 1.0, c.EionHe2s, T, c)
_nat_H0(c::RecfastConstants, h100) = h100 * 100.0 * 1.0e+3 / c.Mpc                                          # cosmology.Recfast.cpp:51
_nat_fHe(c::RecfastConstants, Yp) = Yp / (4.0 * c.fac_mHemH * (1.0 - Yp))                                   # :52
function _nat_NH(c::RecfastConstants, Yp, Omega_b, h100, z)                                                    # :110-114
    mu_H = 1.0 / (1.0 - Yp); H0 = _nat_H0(c, h100)
    return 3.0 * (H0 * H0) * Omega_b / (8.0 * c.pi * c.G * c.mHatom * mu_H) * _cpow(1.0 + z, 3.0)
end

"""Native `evaluate_Recfast_System(z, y, fvec, 4)` for the CAMB-batch configuration (eval_ion_Tr = 0, no DM/B/reionization, constants not varied)."""
function evaluate_recfast_system_native(θ::AbstractVector{Float64}, z::Float64, y, Hz::Float64, c::RecfastConstants = NATIVE_RECFAST_CONSTANTS)
    F = θ[1]; A2s1s = θ[2]; Yp = θ[3]; Omega_b = θ[4]; h100 = θ[5]; T0 = θ[6]
    nHTot = _nat_NH(c, Yp, Omega_b, h100, z)
    fHe = _nat_fHe(c, Yp)
    xHep = Float64(y[1]); xp = Float64(y[2]); TM = Float64(y[3])
    xe = xp + xHep
    nH1 = (1.0 - xp) * nHTot
    nHe1 = nHTot * (fHe - xHep)
    TR = T0 * (1.0 + z)
    Comp = 8.0 * c.sigmaT * c.aRad * _cpow(TR, 4.0) / (3.0 * Hz * (1.0 + z) * c.mElect * c.cLight)
    A2s1sH = A2s1s; A2s1sHe = c.Lam2s1sHe
    lambda21H = c.LyalphaH; lambda21He = c.LyalphaHe
    alphaH = _nat_alphaH(F, TM)
    alphaHe = _nat_alphaHe(TM)
    BH = _nat_Boltzmann(1.0, 1.0, c.E2s1sH, TM, c)
    betaH = _nat_betaH(F, TM, c)
    BHe = _nat_Boltzmann(1.0, 1.0, c.E2s1sHe, TM, c)
    BHe2p2s = 1.0 / _nat_Boltzmann(1.0, 1.0, c.E2p2sHe, TM, c)
    betaHe = _nat_betaHe(TM, c)
    KHe = _cpow(lambda21He, 3.0) / (Hz * 8.0 * c.pi)
    KH = _cpow(lambda21H, 3.0) / (Hz * 8.0 * c.pi)
    CHe = (1.0 + KHe * A2s1sHe * nHe1 * BHe2p2s) / (1.0 + KHe * (A2s1sHe + betaHe) * nHe1 * BHe2p2s)
    CH = (1.0 + KH * A2s1sH * nH1) / (1.0 + KH * (A2s1sH + betaH) * nH1)
    f1 = (alphaHe * xe * xHep * nHTot - betaHe * (fHe - xHep) * BHe) * CHe / ((1.0 + z) * Hz)
    f2 = (alphaH * xe * xp * nHTot - betaH * (1.0 - xp) * BH) * CH / ((1.0 + z) * Hz)
    f3 = Comp * xe / (1.0 + xe + fHe) * (TM - TR) + 2.0 * TM / (1.0 + z)
    f4 = -1.0 / (Hz * (1.0 + z))
    return (f1, f2, f3, f4)
end

"""Native `t_cos_rad(z) = 1.0/2.0/H_z(z)` (cosmology.Recfast.cpp:92-93)."""
t_cos_rad_native(Hz::Float64) = 1.0 / 2.0 / Hz
