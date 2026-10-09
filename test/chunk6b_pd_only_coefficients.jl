# Chunk 6b: the PDE coefficient builder computes `pd` with the private pd-only kernel `_hi_pd_only` instead of the full `hi_rp_rm_pd` (whose Rp/Rm the builder discarded).
# The PUBLIC `hi_pde_coefficients` must stay bitwise identical to a verbatim copy of the previous builder (below, built on the unchanged public `hi_rp_rm_pd`), including promoted
# types, input-domain errors and ForwardDiff derivatives of populations / background / Tg / Te. `get_rates_all` and `hi_rp_rm_pd` are not modified (6a tests cover their native values).
using Test
using ForwardDiff
using CosmoRec
@isdefined(FX6A) || include("chunk6_helpers.jl")
import CosmoRec: cosmos_H, cosmos_NH, cosmos_TCMB

function flat6b!(acc::Vector{UInt64}, x)
    if x isa Float64
        push!(acc, reinterpret(UInt64, x))
    elseif x isa ForwardDiff.Dual
        flat6b!(acc, ForwardDiff.value(x)); foreach(p -> flat6b!(acc, p), ForwardDiff.partials(x).values)
    elseif x isa BigFloat
        push!(acc, hash(string(x)))
    elseif x isa AbstractFloat
        push!(acc, hash(x))
    elseif x isa Integer
        push!(acc, x % UInt64)
    elseif x isa AbstractArray || x isa Tuple
        foreach(e -> flat6b!(acc, e), x)
    elseif isstructtype(typeof(x)) && fieldcount(typeof(x)) > 0
        foreach(f -> flat6b!(acc, getfield(x, f)), fieldnames(typeof(x)))
    else
        push!(acc, hash(x))
    end
    return acc
end
bits6b(x) = flat6b!(UInt64[], x)

# verbatim copy of the pre-increment builder: pd from the public hi_rp_rm_pd (all rates, Rp/Rm discarded)
function reference_builder6b(rows::AbstractMatrix, pops::HIPopulationSplines, cosmos, t::AtomicRateTable, lnBitot::AbstractMatrix, lv::HIPDELevels;
                             zs::Real, ze::Real, h_kb::Float64 = 4.7992373449498863e-11, pi_::Float64 = 3.1415926535897931)
    lo = ze / 1.0001; hi = zs * 1.0001
    nz = count(r -> CosmoRec._primal(r) >= lo, view(rows, :, 1))
    zg = CosmoRec.init_xarr_linear(lo, hi, nz)
    n = length(lv.index)
    T = promote_type(typeof(hi_Xi(pops, zg[1], 0)), typeof(cosmos_H(cosmos, zg[1])), typeof(cosmos_NH(cosmos, zg[1])))
    lnpd = [Vector{T}(undef, nz) for _ in 1:n]; Dn = [Vector{T}(undef, nz) for _ in 1:n]
    for k in 1:nz
        z = zg[k]
        Tg = cosmos_TCMB(cosmos, z); Te = Tg * hi_rho(pops, z); NH = cosmos_NH(cosmos, z); Hz = cosmos_H(cosmos, z)
        _, pd = hi_rp_rm_pd(z, pops, t, lnBitot, lv, Tg, Te, NH * hi_Xe(pops, z), 1.0 - hi_Xi(pops, z, 0), NH)
        N1s = NH * hi_Xi(pops, z, 0)
        for m in 1:n
            i = lv.index[m]
            ex = exp(-h_kb * lv.Dnu_1s[m] / Tg)
            Ni = NH * hi_Xi(pops, z, i)
            if i == 1
                Dn[m][k] = Ni / N1s - ex
            elseif lv.l[m] == 1
                nL = Ni / 3.0 / N1s
                tauS = lv.A21[m] * lv.lambda21[m]^3 / (8.0 * pi_ * Hz) * (N1s * 3.0 - Ni)
                PS = (1.0 - exp(-tauS)) / tauS
                Dn[m][k] = (nL - ex) * (1.0 + (1.0 / pd[m] - 1.0) * PS)
            else
                Dn[m][k] = Ni / N1s / (2.0 * lv.l[m] + 1.0) - ex
            end
            lnpd[m][k] = log(pd[m])
        end
    end
    return CosmoRec.HIPDECoefficientSplines(zg, [natural_cubic_spline(zg, lnpd[m]) for m in 1:n], [natural_cubic_spline(zg, Dn[m]) for m in 1:n])
end

struct TScaled6b{A,T}      # scales Tg (hence Te = Tg rho)
    a::A
    s::T
end
cosmos_H(b::TScaled6b, z) = cosmos_H(b.a, z)
cosmos_NH(b::TScaled6b, z) = cosmos_NH(b.a, z)
cosmos_TCMB(b::TScaled6b, z) = b.s * cosmos_TCMB(b.a, z)

build6b(rows, bg) = (p = HIPopulationSplines(rows; zs = ZS6, ze = ZE6); hi_pde_coefficients(rows, p, bg, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6))
refbuild6b(rows, bg) = (p = HIPopulationSplines(rows; zs = ZS6, ze = ZE6); reference_builder6b(rows, p, bg, HTAB6, LNBITOT6, NATIVE_HI_PDE_LEVELS; zs = ZS6, ze = ZE6))
outcome6b(f) = try; ("ok", bits6b(f())); catch e; ("err", string(typeof(e)), first(sprint(showerror, e), 300)); end
outs6b(cs) = reduce(vcat, [vcat([log(hi_pd(cs, z, m)) for m in (1, 2, 4)], [hi_Dnem(cs, z, m) for m in 1:5]) for z in (2300.0, 1700.0, 1300.0, 900.0, 600.0)])

@testset "Chunk 6b: PDE coefficient builder with the pd-only kernel" begin
    nat = NODES5[:, 1:9]
    pops = HIPopulationSplines(nat; zs = ZS6, ze = ZE6)
    lv = NATIVE_HI_PDE_LEVELS

    @testset "public builder vs unchanged-reference builder: all arrays bitwise" begin
        for (name, rows, bg) in (("native", nat, ACC5), ("scaled background", nat, CosmoRec.ScaledBackground(ACC5, 1.07, 0.93)), ("Tg scaled", nat, TScaled6b(ACC5, 1.0021)),
                                 ("perturbed populations", hcat(nat[:, 1], nat[:, 2:9] .* (1 .+ 1.0e-3 .* sin.(reshape(1:(8 * size(nat, 1)), size(nat, 1), 8)))), ACC5))
            a = build6b(rows, bg); b = refbuild6b(rows, bg)
            @test typeof(a) == typeof(b)
            @test bits6b(a) == bits6b(b)
            @test all(m -> bits6b([hi_pd(a, z, m) for z in a.z]) == bits6b([hi_pd(b, z, m) for z in b.z]) && bits6b([hi_Dnem(a, z, m) for z in a.z]) == bits6b([hi_Dnem(b, z, m) for z in b.z]), 1:5)
        end
    end

    @testset "promoted types: Dual populations, Dual background, Dual Tg, BigFloat" begin
        dual(x) = ForwardDiff.Dual(x, 1.0)
        for (rows, bg) in ((map(dual, nat), ACC5), (nat, CosmoRec.ScaledBackground(ACC5, ForwardDiff.Dual(1.0, 1.0, 0.0), ForwardDiff.Dual(1.0, 0.0, 1.0))))
            a = build6b(rows, bg); b = refbuild6b(rows, bg)
            @test typeof(a) == typeof(b)
            @test bits6b(a) == bits6b(b)
        end
        # a Dual Tg alone is unsupported by the builder (its element type excludes Tg): the outcome (error) must be unchanged
        @test outcome6b(() -> build6b(nat, TScaled6b(ACC5, ForwardDiff.Dual(1.0, 1.0))))[1:2] == outcome6b(() -> refbuild6b(nat, TScaled6b(ACC5, ForwardDiff.Dual(1.0, 1.0))))[1:2] == ("err", "MethodError")
        setprecision(BigFloat, 256) do
            a = build6b(BigFloat.(nat), ACC5); b = refbuild6b(BigFloat.(nat), ACC5)
            @test typeof(a) == typeof(b)
            @test bits6b(a) == bits6b(b)
        end
    end

    @testset "pd-only kernel vs hi_rp_rm_pd: values, element types, Jacobians" begin
        for z in range(500.0, 2500.0; length = 11)
            NH = cosmos_NH(ACC5, z); Ne = NH * hi_Xe(pops, z); Xp = 1.0 - hi_Xi(pops, z, 0); Tg = cosmos_TCMB(ACC5, z); Te = Tg * hi_rho(pops, z)
            old = hi_rp_rm_pd(z, pops, HTAB6, LNBITOT6, lv, Tg, Te, Ne, Xp, NH)[2]
            new = CosmoRec._hi_pd_only(z, pops, HTAB6, LNBITOT6, lv, Tg, Te, Ne, NH)
            @test eltype(new) == eltype(old) && bits6b(new) == bits6b(old)
            @test new[3] == 1.0 && new[5] == 1.0
            fo(q) = hi_rp_rm_pd(z, pops, HTAB6, LNBITOT6, lv, q[1], q[2], Ne, Xp, NH)[2]
            fn(q) = CosmoRec._hi_pd_only(z, pops, HTAB6, LNBITOT6, lv, q[1], q[2], Ne, NH)
            @test bits6b(ForwardDiff.jacobian(fn, [Tg, Te])) == bits6b(ForwardDiff.jacobian(fo, [Tg, Te]))
            # Dual only in Tg, only in Te, in both, and a Dual Ne/NH: element type and values follow the original promotion
            D = ForwardDiff.Dual
            for (tg, te, ne, nh) in ((D(Tg, 1.0), Te, Ne, NH), (Tg, D(Te, 1.0), Ne, NH), (D(Tg, 1.0), D(Te, 1.0), Ne, NH), (Tg, Te, D(Ne, 1.0), NH), (Tg, Te, Ne, D(NH, 1.0)))
                o = hi_rp_rm_pd(z, pops, HTAB6, LNBITOT6, lv, tg, te, ne, Xp, nh)[2]; n_ = CosmoRec._hi_pd_only(z, pops, HTAB6, LNBITOT6, lv, tg, te, ne, nh)
                @test eltype(o) == eltype(n_) && bits6b(o) == bits6b(n_)
            end
            setprecision(BigFloat, 256) do
                o = hi_rp_rm_pd(z, pops, HTAB6, LNBITOT6, lv, big(Tg), big(Te), big(Ne), big(Xp), big(NH))[2]
                n_ = CosmoRec._hi_pd_only(z, pops, HTAB6, LNBITOT6, lv, big(Tg), big(Te), big(Ne), big(NH))
                @test eltype(o) == eltype(n_) && bits6b(o) == bits6b(n_)
            end
        end
    end

    @testset "input domain: same exceptions as the full rate path" begin
        z = 1500.0; NH = cosmos_NH(ACC5, z); Ne = NH * hi_Xe(pops, z); Xp = 1.0 - hi_Xi(pops, z, 0); Tg = cosmos_TCMB(ACC5, z)
        for Te in (Tg * 1.0e3, Tg * 1.0e-3, NaN, Inf, -1.0, 0.0), Tg_ in (Tg, NaN, -2.0, 1.0e7)
            o = outcome6b(() -> hi_rp_rm_pd(z, pops, HTAB6, LNBITOT6, lv, Tg_, Te, Ne, Xp, NH)[2])
            n_ = outcome6b(() -> CosmoRec._hi_pd_only(z, pops, HTAB6, LNBITOT6, lv, Tg_, Te, Ne, NH))
            @test o == n_
        end
        hot = copy(nat); hot[:, 9] .*= 1.0e3
        cold = copy(nat); cold[:, 9] .= 1.0 + 1.0e-6
        neg = copy(nat); neg[:, 9] .= 0.5
        nanr = copy(nat); nanr[1500, 9] = NaN
        for r in (hot, cold, neg, nanr)
            @test outcome6b(() -> build6b(r, ACC5)) == outcome6b(() -> refbuild6b(r, ACC5))
        end
        @test_throws CosmoRec.RateTableDomainError build6b(hot, ACC5)
    end

    @testset "ForwardDiff of the public builder: populations, rho (Te), background, Tg" begin
        x0 = vec(nat[:, 2:9]); mk(x) = hcat(nat[:, 1], reshape(x, size(nat, 1), 8))
        v = [sin(0.37 * k) for k in eachindex(x0)] .* abs.(x0)
        vr = zeros(length(x0)); vr[(7 * size(nat, 1) + 1):end] .= nat[:, 9] .- 1.0
        for dir in (v, vr)
            dn = ForwardDiff.derivative(t -> outs6b(build6b(mk(x0 .+ t .* dir), ACC5)), 0.0)
            do_ = ForwardDiff.derivative(t -> outs6b(refbuild6b(mk(x0 .+ t .* dir), ACC5)), 0.0)
            @test all(isfinite, dn) && bits6b(dn) == bits6b(do_)
        end
        for p in ([1.0, 1.0], [1.03, 0.97])
            jn = ForwardDiff.jacobian(q -> outs6b(build6b(nat, CosmoRec.ScaledBackground(ACC5, q[1], q[2]))), p)
            jo = ForwardDiff.jacobian(q -> outs6b(refbuild6b(nat, CosmoRec.ScaledBackground(ACC5, q[1], q[2]))), p)
            @test all(isfinite, jn) && bits6b(jn) == bits6b(jo)
        end
        # Tg/Te derivatives of pd itself are in the kernel testset; a Dual Tg-scale through the builder is an unsupported (erroring) case with unchanged outcome
        @test outcome6b(() -> ForwardDiff.derivative(s -> outs6b(build6b(nat, TScaled6b(ACC5, s))), 1.0))[1:2] == outcome6b(() -> ForwardDiff.derivative(s -> outs6b(refbuild6b(nat, TScaled6b(ACC5, s))), 1.0))[1:2] == ("err", "MethodError")
    end
end
