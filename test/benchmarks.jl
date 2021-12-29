using Test, BenchmarkTools
using TripleCorrelations

include("src/helpers.jl")
include("src/old_implementations.jl")

@testset "Sequence classes (N=100,l=7)" begin 
    N=100; l=7
    test_raster = generate_random_raster((N,N))
    test_raster_Int = Matrix{Int}(test_raster)

    contributions_calcd_by_loop = @btime sequence_class_tricorr($test_raster, $l, $l)
    _OLD_contributions_calcd_unrolld = @btime _OLD_sequence_class_tricorr_unrolled($test_raster_Int, $l, $l)
    contributions_calcd_unrolld = @btime sequence_class_tricorr_unrolled($test_raster_Int, $l, $l)

    @test contributions_calcd_by_loop == _OLD_contributions_calcd_unrolld
    @test contributions_calcd_by_loop == contributions_calcd_unrolld
    @test _OLD_contributions_calcd_unrolld == contributions_calcd_unrolld
end


# @testset "Sequence classes (N=500,l=7)" begin 
#     test_raster = generate_random_raster((500,500))

#     contributions_calcd_by_loop = @btime sequence_class_tricorr($test_raster, 7,7)
#     contributions_calcd_unrolld = @btime sequence_class_tricorr_unrolled($test_raster, 7,7)

#     @test contributions_calcd_by_loop == contributions_calcd_unrolld
# end

# @testset "Sequence classes (N=100,l=15)" begin 
#     N=100; l=15
#     test_raster = generate_random_raster((N,N))

#     contributions_calcd_by_loop = @btime sequence_class_tricorr($test_raster, $l, $l)
#     contributions_calcd_unrolld = @btime sequence_class_tricorr_unrolled($test_raster, $l, $l)

#     @test contributions_calcd_by_loop == contributions_calcd_unrolld
# end