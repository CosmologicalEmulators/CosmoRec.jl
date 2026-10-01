using Test

@testset "CosmoRec.jl" begin
    include("chunk1a_stiff_ad_probe.jl")
    include("chunk1b_radiation_shaped_probe.jl")
    include("chunk1c_nonlinear_probe.jl")
    include("chunk1d_native_fixture.jl")
    include("chunk2_rate_table_native.jl")
    include("chunk2_rate_table_ad.jl")
    include("chunk3a_rhs_native.jl")
    include("chunk3a_rhs_ad.jl")
    include("chunk3b_rhs_native.jl")
    include("chunk3b_rhs_ad.jl")
    include("chunk3c_hiabs_native.jl")
    include("chunk3c_hiabs_ad.jl")
    include("chunk3d_fcn_native.jl")
    include("chunk3d_fcn_ad.jl")
    include("chunk3e_dpesc_native.jl")
    include("chunk3e_dpesc_ad.jl")
end
