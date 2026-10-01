using Test
using SHA

# Natural cubic-free 4-point Lagrange interpolation (Base only) on the "grid" rows.
function lagrange4(x, y, xq)
    i = clamp(searchsortedlast(x, xq), 2, length(x) - 2)
    idx = i-1:i+2
    s = 0.0
    for a in idx
        w = 1.0
        for b in idx
            b == a || (w *= (xq - x[b]) / (x[a] - x[b]))
        end
        s += w * y[a]
    end
    return s
end

function read_fixture(path)
    meta = Dict{String,String}()
    cols = String[]
    rows = Vector{Vector{String}}()
    for line in eachline(path)
        if startswith(line, "#")
            m = match(r"^#\s*([^:]+):\s*(.*)$", line)
            m === nothing || (meta[String(m[1])] = String(m[2]))
        elseif isempty(cols)
            cols = split(line)
        else
            push!(rows, split(line))
        end
    end
    return meta, cols, rows
end

@testset "Chunk 1d: native CosmoRec output fixture" begin
    path = joinpath(@__DIR__, "fixtures", "native_camb_cosmorec_thermo.txt")
    meta, cols, rows = read_fixture(path)

    @testset "metadata and hash" begin
        @test meta["source_camb_cosmorec_sha"] == "fa3f097343fbbe427cc04b4f5f0041c22c6ec764"
        @test meta["source_cosmorec_sha"] == "086769055f61ae0c244a53dd381ee65b624d0ac3"
        @test startswith(meta["recombination_model"], "CosmoRec")
        @test haskey(meta, "units")
        @test occursin("km/s/Mpc", meta["units"])
        @test occursin("OUTPUTS ONLY", read(path, String))
        body = join(join.(rows, " ") .* "\n")
        for r in rows
            @test length(r) == length(cols)
        end
        @test parse(Int, meta["n_rows"]) == length(rows)
        @test bytes2hex(sha256(body)) == meta["sha256_rows"]
    end

    @test cols == ["z", "role", "x_e", "T_b", "H", "x_e_noreion", "T_b_noreion"]
    z = parse.(Float64, getindex.(rows, 1))
    role = getindex.(rows, 2)
    col(name) = parse.(Float64, getindex.(rows, findfirst(==(name), cols)))
    xe, Tb, H, xen, Tbn = col("x_e"), col("T_b"), col("H"), col("x_e_noreion"), col("T_b_noreion")

    @testset "grid" begin
        @test issorted(z) && allunique(z)
        @test z[1] == 0.0 && z[end] == 10000.0
        for zz in (50.0, 200.0, 3000.0, 9999.0, 10000.0)
            @test zz in z
        end
        @test all(in(("grid", "holdout")), role)
        @test count(==("holdout"), role) >= 50
        # coverage of the key regions
        for (a, b) in ((0, 50), (50, 200), (550, 1500), (1500, 2500), (2500, 4000), (4000, 7500), (7500, 10000))
            @test count(zi -> a <= zi <= b, z) >= 5
        end
    end

    @testset "finiteness, range, shape" begin
        for v in (xe, Tb, H, xen, Tbn)
            @test all(isfinite, v) && all(>(0), v)
        end
        @test issorted(H)             # H grows with z
        @test issorted(Tb)            # T_b grows with z
        @test H[1] ≈ 67.36 rtol = 1e-12
        @test all(xe .≈ xen .|| z .< 30)   # reionization only affects low z
        @test xe[1] ≈ 1.1640227763583468 rtol = 1e-14
        @test xen[1] < 1e-3
        @test 1.0 < xe[end] < 1.2                 # fully ionised H + He at zmax
        @test maximum(xen[(z .> 50) .& (z .< 200)]) < 1e-3  # freeze-out region
        i1500, i800 = findfirst(==(1500.0), z), findfirst(==(800.0), z)
        @test xe[i800] < 0.2 && xe[i1500] > 0.9    # hydrogen recombination between z~1500 and ~800
        # T_b ≈ T_CMB (1+z) at high z, within 1e-2
        @test Tb[end] ≈ 2.7255 * 10001 rtol = 5e-3
        # standard spot values (original native, this fixture only)
        for (zz, ex, eT, eH) in ((50.0, 0.00023887789169251859, 50.702366072888417, 13865.766597033085),
                                 (200.0, 0.00033834606340520134, 466.4585126318799, 110728.41760713251),
                                 (3000.0, 1.0819197602184856, 8179.2237292156105, 8509774.2573337965))
            i = findfirst(==(zz), z)
            @test xe[i] ≈ ex rtol = 1e-14
            @test Tb[i] ≈ eT rtol = 1e-14
            @test H[i] ≈ eH rtol = 1e-14
        end
    end

    @testset "off-grid interpolation vs independent native queries" begin
        g = role .== "grid"
        h = role .== "holdout"
        # Interpolate in log(1+z) for H and T_b, in x_e (no-reion) linearly spaced nodes
        lz = log.(1 .+ z)
        errs = Dict("H" => 0.0, "Tb" => 0.0, "xe" => 0.0)
        for i in findall(h)
            errs["H"] = max(errs["H"], abs(exp(lagrange4(lz[g], log.(H[g]), lz[i])) / H[i] - 1))
            errs["Tb"] = max(errs["Tb"], abs(exp(lagrange4(lz[g], log.(Tb[g]), lz[i])) / Tb[i] - 1))
            errs["xe"] = max(errs["xe"], abs(lagrange4(z[g], xen[g], z[i]) - xen[i]) / max(xen[i], 1e-3))
        end
        @info "fixture interpolation max errors" errs
        @test errs["H"] < 1e-3
        @test errs["Tb"] < 5e-3
        @test errs["xe"] < 2.5e-3
    end
end
