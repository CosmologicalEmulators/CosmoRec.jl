# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Chunk 7b: coefficients of the HI radiation PDE  du/dz = A u'' + B u' + C u + D  (u = nu21 Dn(x), x = nu/nu_Lya).
# Source (ORIGINAL CosmoRec v3.0b, read-only): PDE_Problem/define_PDE.cpp:275-1297 (def_PDE_Lyn_and_2s1s with the production flags: line scattering,
# electron scattering, line emission/absorption and 2s-1s on; recoil and diffusion on; use_full_width = 1; e_sc_efficiency = 1; no nD-1s quadrupole lines),
# Development/Line_profiles/Voigtprofiles.cpp (Voigtprofile_Dawson::dphi_dx).
# ----------------------------------------------------------------------------------------------------------------------------------------------------

"""Native constants of `physical_consts.h` used by `define_PDE.cpp` (`const_sigT` is the CosmoRec.cpp run-time value)."""
struct HIPDEConstants
    h_kb::Float64
    kB::Float64
    me_gr::Float64
    cl::Float64
    mH_gr::Float64
    sigT::Float64
end
const NATIVE_HI_PDE_CONSTANTS = HIPDEConstants(4.7992373449498863e-11, 1.3806504000000002e-16, 9.1093821500000007e-28, 29979245800.0, 1.673532551272445e-24, 6.6524585583988891e-25)

"""HI atom data used by the production PDE (Atom.cpp `Gas_of_Atoms`, nShells = 3): Ly-n Voigt lines (n = 2, 3: `nu21`, `Gamma` = A_tot(np), A21(np-1s)),
`nu21` = Dnu_1s(2p) (frequency unit), Dnu_1s(n) and the np<->(n-1)s/d rates entering the Raman/two-photon death-probability corrections."""
struct HIPDEAtom
    lyn::Vector{VoigtLine}
    lyn_A21::Vector{Float64}
    nu21::Float64
    Dnu_1s::Vector{Float64}           # by principal quantum number n = 1..3 (independent of l in the native atom)
    A_3p2s::Float64; nu_3p2s::Float64 # Level(3,1).A21(2,0), nu21(2,0)
    A_3s2p::Float64; nu_3s2p::Float64 # Level(3,0).A21(2,1), nu21(2,1)
    A_3d2p::Float64; nu_3d2p::Float64 # Level(3,2).A21(2,1), nu21(2,1)
    lyn_lambda21::Vector{Float64}     # Ly-n wavelengths (tau_S of the correction integrals)
end
const NATIVE_HI_PDE_ATOM = HIPDEAtom([VoigtLine(2466038423768827.0, 626490298.80265749, 1.0), VoigtLine(2922712205948239.0, 189701000.25209922, 1.0)],
                                     [626490298.80265749, 167252720.15153983], 2466038423768827.0, [0.0, 2466038423768827.0, 2922712205948239.0],
                                     22448280.100559387, 456673782179412.38, 6313578.7782823257, 456673782179412.38, 64651046.689611033, 456673782179412.38,
                                     [1.2156844561319915e-05, 1.025733759861368e-05])

"""`Voigtprofile_Dawson::dphi_dx(x, a)` (Voigtprofiles.cpp): derivative of the Mihalas expansion, wing expansion for |x| >= 30."""
function voigt_dphi_dx(x, a)
    if abs(x) >= DPESC_X_WING
        x2 = x * x; a2 = a * a
        return -2.0 * a / (x * x2 * 3.1415926535897931) * (1.0 + (2.0 * (1.5 - a2) + (15.0 * (0.75 - a2) + 105.0 * (0.5 - a2) / x2) / x2) / x2)
    end
    x2 = x * x
    F = dawson_jc(x)
    s = DPESC_SQRT_PI
    d0 = -2.0 * x * exp(-x2)
    d1 = 4.0 / s * (F + x - 2.0 * F * x2)
    d2 = d0 * (3.0 - 2.0 * x2)
    d3 = 4.0 / 3.0 / s * ((5.0 - 2.0 * x2) * x + F * (3.0 + 4.0 * (x2 - 3.0) * x2))
    d4 = d0 / 6.0 * (15.0 - 4.0 * (5.0 - x2) * x2)
    d5 = 2.0 / 15.0 / s * ((33.0 + 4.0 * (x2 - 7.0) * x2) * x + F * (15.0 - (90.0 - (60.0 - 8.0 * x2) * x2) * x2))
    d6 = d0 / 90.0 * (105.0 - (210.0 - (84.0 - 8.0 * x2) * x2) * x2)
    d7 = 2.0 / 315.0 / s * ((279.0 - (370.0 - (108.0 - 8.0 * x2) * x2) * x2) * x + F * (105.0 - (840.0 - (840.0 - (224.0 - 16.0 * x2) * x2) * x2) * x2))
    d8 = d0 / 315.0 * (945.0 / 8.0 - (315.0 - (189.0 - (36.0 - 2.0 * x2) * x2) * x2) * x2)
    r = d0 + (d1 + (d2 + (d3 + (d4 + (d5 + (d6 + (d7 + d8 * a) * a) * a) * a) * a) * a) * a) * a
    return r / s
end

"""Everything the native PDE right-hand side reads: the armed setup (grid, ratios), the population and pd/Dnem splines of the current history, the
explicit Cosmos accessors (loaded H), atom data and constants."""
struct HIPDEModel{P,C,A}
    setup::HIPDESetup
    pops::P
    coef::C
    cosmos::A
    atom::HIPDEAtom
    k::HIPDEConstants
end
HIPDEModel(setup::HIPDESetup, pops, coef, cosmos; atom = NATIVE_HI_PDE_ATOM, k = NATIVE_HI_PDE_CONSTANTS) = HIPDEModel(setup, pops, coef, cosmos, atom, k)

_hi_level(n::Int, l::Int) = n * (n - 1) ÷ 2 + l          # native HI index (1s = 0); equals the HIPDELevels position for 2s..3d
_pde_Dnem(m::HIPDEModel, z, n, l) = hi_Dnem(m.coef, z, _hi_level(n, l))
_pde_pd(m::HIPDEModel, z, n, l) = hi_pd(m.coef, z, _hi_level(n, l))
_pde_Xnl(m::HIPDEModel, z, n, l) = hi_Xi(m.pops, z, _hi_level(n, l))

"""Work arrays and the side products of one `def_PDE_Lyn_and_2s1s` call (native `PDE_funcs` globals) needed by the correction integrals."""
struct HIPDEState{T}
    A::Vector{T}; B::Vector{T}; C::Vector{T}; D::Vector{T}
    exp_x::Vector{T}
    phi::Vector{Vector{T}}        # Ly-n Voigt profiles, n = 2, 3
    hp_phi::Vector{Vector{T}}     # Voigt profiles at x_res(n) - x (two-photon low-frequency side)
    aV::Vector{T}; Dnu::Vector{T}; pd::Vector{T}; pd_eff::Vector{T}; Dnem_eff::Vector{T}
    P2g_ns::Vector{Vector{T}}; P2g_nd::Vector{Vector{T}}; PR_ns::Vector{Vector{T}}; PR_nd::Vector{Vector{T}}
end
function HIPDEState{T}(npts::Int) where {T}
    v() = zeros(T, npts)
    return HIPDEState{T}(v(), v(), v(), v(), v(), [v(), v()], [v(), v()], zeros(T, 2), zeros(T, 2), zeros(T, 2), zeros(T, 2), zeros(T, 2),
                         [v(), v()], [v(), v()], [v(), v()], [v(), v()])
end

# Dp_Dnem_function (define_PDE.cpp:617-715) for the production processes (Raman up to 2s, two-photon up to n = 3)
function _pde_Dp_Dn(m::HIPDEModel, n::Int, z, Tg, pd_np, nem, N1s)
    a = m.atom; h_kb = m.k.h_kb
    NH = cosmos_NH(m.cosmos, z)
    DRp = zero(N1s); DRm = zero(N1s)
    if n > 2       # Raman: ni = min(7, n-1) .. 3 is empty for n <= 3; the 2s term for n = 3
        Nj = _pde_Xnl(m, z, 2, 0) * NH
        ex = exp(-h_kb * a.nu_3p2s / Tg); npij = ex / (1.0 - ex)
        DRm += a.A_3p2s * (1.0 + npij)
        DRp += 3.0 / (2.0 * 0 + 1.0) * Nj * a.A_3p2s * npij
    end
    DRp2 = zero(N1s); DRm2 = zero(N1s)
    if n == 2      # two-photon: ni = min(8, 3) .. n+1
        for (li, A, nu) in ((0, a.A_3s2p, a.nu_3s2p), (2, a.A_3d2p, a.nu_3d2p))
            Nj = _pde_Xnl(m, z, 3, li) * NH
            ex = exp(-h_kb * nu / Tg); npij = ex / (1.0 - ex)
            DRm2 += (2.0 * li + 1.0) / 3.0 * A * npij
            DRp2 += A * (1.0 + npij) * Nj
        end
    end
    DRp += DRp2; DRm += DRm2
    DRm == 0.0 && return (zero(N1s), zero(N1s))
    A21 = a.lyn_A21[n - 1]
    Rm_tot = A21 / (1.0 - pd_np)
    Rm = A21 * pd_np / (1.0 - pd_np)
    Dp = DRm / Rm_tot
    Rp = 3.0 * nem * N1s * Rm
    Dn = ((Rp - DRp) - 3.0 * nem * N1s * (Rm - DRm)) / (Rm - DRm) / 3.0 / N1s
    return (Dp, Dn)
end

_two_g_emission_fac(exp_mx, exp_mx0) = (1.0 + exp_mx0 / (exp_mx - exp_mx0)) * (1.0 + exp_mx / (1.0 - exp_mx))

"""
    hi_pde_rhs_coefficients!(st, m, z)

Native `def_PDE_Lyn_and_2s1s(z, x, A, B, C, D)` for the production configuration: Hubble term, electron scattering (Kompaneets diffusion + recoil), Ly-2
and Ly-3 resonance scattering and emission/absorption with Raman/two-photon corrected death probabilities, the 2s-1s Raman term, the 3s/3d-1s and
2s-1s two-photon terms. Fills `st` (coefficients and the `PDE_funcs` side products) and returns it.
"""
function hi_pde_rhs_coefficients!(st::HIPDEState, m::HIPDEModel, z)
    s = m.setup; a = m.atom; k = m.k; x = s.x; npts = length(x)
    nu21 = a.nu21
    Hz = cosmos_H(m.cosmos, z); NH = cosmos_NH(m.cosmos, z); Tg = cosmos_TCMB(m.cosmos, z)
    Ne = NH * hi_Xe(m.pops, z)
    N1s = NH * hi_Xi(m.pops, z, 0)
    rho = hi_rho(m.pops, z)
    Te = Tg * rho
    @inbounds for i in 1:npts
        st.exp_x[i] = exp(-k.h_kb * nu21 * x[i] / Tg)
        st.A[i] = 0.0; st.C[i] = 0.0; st.D[i] = 0.0
        st.B[i] = -x[i] / (1.0 + z)
    end
    # electron scattering
    xifac = k.h_kb * nu21 / Te
    kap_e = -k.sigT * Ne * (k.kB * Te / k.me_gr / k.cl) / Hz / (1.0 + z)
    kap_e *= 1.0
    @inbounds for i in 1:npts
        D_x = kap_e * x[i] * x[i]
        dlnD = D_x * 4.0 / x[i]
        st.A[i] += D_x; st.B[i] += dlnD; st.B[i] += D_x * xifac; st.C[i] += dlnD * xifac
    end
    # Ly-n scattering + emission/absorption
    Tm = rho * Tg
    for n in 2:3
        L = a.lyn[n - 1]; j = n - 1
        st.pd[j] = _pde_pd(m, z, n, 1)
        st.Dnu[j] = sqrt(2.0 * DPESC_K_MHC2 * Tm / L.AM)
        psc = 1.0 - st.pd[j]
        st.aV[j] = voigt_a(L, Tm) / L.Gamma * a.lyn_A21[j] / psc
        aV = st.aV[j]; DnT = voigt_DnuT(L, Tm)
        @inbounds for i in 1:npts
            st.phi[j][i] = voigt_phi((nu21 * x[i] - L.nu21) / DnT, aV)
            st.hp_phi[j][i] = voigt_phi((nu21 * ((1.0 - 1.0 / n / n) / 0.75 - x[i]) - L.nu21) / DnT, aV)
        end
        # add_Ly_n_sc_em_full
        pd_np = st.pd[j]; psc_np = 1.0 - pd_np
        lambda = k.cl / nu21
        DnuD_nu = st.Dnu[j]
        const_sigr = 1.5 * lambda * lambda * (a.lyn_A21[j] / PDE_FOURPI / (st.Dnu[j] * L.nu21))
        kap_sc = -psc_np * const_sigr * N1s * (k.kB * Te / k.mH_gr / k.cl) / Hz / (1.0 + z)
        nun1 = a.Dnu_1s[n]
        x_c = k.h_kb * nun1 / Tg
        Dnem = _pde_Dnem(m, z, n, 1) * nu21
        kap_abs = -pd_np * const_sigr * N1s * k.cl / Hz / (1.0 + z)
        exp_dum = exp(-x_c)
        Dp, Dn = _pde_Dp_Dn(m, n, z, Tg, pd_np, Dnem / nu21 + exp_dum, N1s)
        kap_abs *= 1.0 - Dp / pd_np
        pd_np *= 1.0 - Dp / pd_np
        Dnem += Dn * nu21
        st.pd_eff[j] = pd_np; st.Dnem_eff[j] = Dnem / nu21
        DnTe = voigt_DnuT(L, Te)
        @inbounds for i in 1:npts
            phi_x = st.phi[j][i]
            xD = (nu21 * x[i] - L.nu21) / DnTe
            dphi = voigt_dphi_dx(xD, aV) / DnuD_nu * nu21 / nun1
            D_x = kap_sc * phi_x
            dlnD = kap_sc * (phi_x * 2.0 / x[i] + dphi)
            st.A[i] += D_x; st.B[i] += dlnD; st.B[i] += D_x * xifac; st.C[i] += dlnD * xifac
            exp_fac = exp_dum / st.exp_x[i]
            D_x = kap_abs * phi_x / x[i] / x[i]
            st.C[i] += -D_x * exp_fac
            st.D[i] += D_x * Dnem
        end
    end
    # Raman 2s-1s profile (update_Raman_profiles_nsd_1s(z, 2, ...)) and term
    sig_fac = (k.cl / nu21)^2
    let ni = 2, P = st.PR_ns[1]
        x_c = k.h_kb * a.Dnu_1s[ni] / Tg; exp_dum = exp(-x_c)
        klow = s.index_emission[1]
        dumf = (1.0 - st.pd[2]) * s.A_npns[1][2] / (st.Dnu[2] * a.lyn[2].nu21)
        @inbounds for i in (klow + 1):npts
            P[i] = 0.0
            P[i] += st.phi[2][i] * dumf
            P[i] *= s.ratio_R_ns[1][i]
            P[i] *= st.exp_x[i] / (exp_dum - st.exp_x[i])
        end
        Dnem_nsd = _pde_Dnem(m, z, 2, 0) * nu21
        x_c = k.h_kb * a.Dnu_1s[2] / Tg; exp_dum = exp(-x_c)
        kap = -((2.0 * 0 + 1.0) / 2.0 / PDE_FOURPI * sig_fac) * N1s * k.cl / Hz / (1.0 + z)
        @inbounds for i in (klow + 1):npts
            npl = st.exp_x[i] / (1.0 - st.exp_x[i])
            DD = kap * P[i] / x[i] / x[i]
            st.C[i] += -DD * (exp_dum / npl)
            st.D[i] += DD * Dnem_nsd
        end
    end
    # two-photon profiles: 2s-1s, then 3s/3d-1s
    let P = st.P2g_ns[1]
        nun1 = a.Dnu_1s[2]; exp_dum = exp(-(k.h_kb * nun1 / Tg))
        @inbounds for i in 1:s.index_emission[1]
            P[i] = s.ratio_2g_ns[1][i] * _two_g_emission_fac(st.exp_x[i], exp_dum) / nun1
        end
    end
    for (P, ratio, Avec) in ((st.P2g_ns[2], s.ratio_2g_ns[2], s.A_npns[2]), (st.P2g_nd[2], s.ratio_2g_nd[2], s.A_npnd[2]))
        nun1 = a.Dnu_1s[3]; exp_dum = exp(-(k.h_kb * nun1 / Tg))
        kmax = s.index_emission[2]
        dumf = (1.0 - st.pd[1]) * Avec[1] / (st.Dnu[1] * a.lyn[1].nu21)
        i2 = s.index_2 + 1
        @inbounds for i in 1:kmax
            P[i] = 0.0
            dum = st.phi[1][i]
            dum += i >= i2 ? st.hp_phi[1][i] : st.hp_phi[1][i2]
            P[i] += dum * dumf
            P[i] *= ratio[i]
            P[i] *= _two_g_emission_fac(st.exp_x[i], exp_dum)
        end
    end
    for (n, l, P) in ((3, 0, st.P2g_ns[2]), (3, 2, st.P2g_nd[2]), (2, 0, st.P2g_ns[1]))
        kap = -((2.0 * l + 1.0) / 2.0 / PDE_FOURPI * sig_fac) * N1s * k.cl / Hz / (1.0 + z)
        Dnem_nsd = _pde_Dnem(m, z, n, l) * nu21
        exp_dum = exp(-(k.h_kb * a.Dnu_1s[n] / Tg))
        @inbounds for i in 1:s.index_emission[n - 1]
            npl = st.exp_x[i] / (1.0 - st.exp_x[i])
            DD = kap * P[i] / x[i] / x[i]
            st.C[i] += -DD * (exp_dum / npl)
            st.D[i] += DD * Dnem_nsd
        end
    end
    return st
end

hi_pde_rhs_coefficients(m::HIPDEModel, z) =
    hi_pde_rhs_coefficients!(HIPDEState{promote_type(typeof(z), typeof(hi_Xi(m.pops, z, 0)), typeof(_pde_Dnem(m, z, 2, 1)), Float64)}(length(m.setup.x)), m, z)
