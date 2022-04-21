using Test, BenchmarkTools
using TripleCorrelations

include("src/old_implementations.jl")

@testset "Sequence classes (N=100,l=7)" begin 
    N=100; l=7
    test_raster = generate_random_raster((N,N))
    test_raster_Int = Matrix{Int}(test_raster)

    contributions_calcd_by_loop = @btime sequence_class_tricorr($test_raster, ($l, $l))
    _OLD_contributions_calcd_unrolled = @btime _OLD_sequence_class_tricorr_unrolled($test_raster_Int, ($l, $l))

    @test contributions_calcd_by_loop == _OLD_contributions_calcd_unrolled
end