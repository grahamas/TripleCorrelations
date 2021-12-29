using TripleCorrelations
using Test, BenchmarkTools

include("src/helpers.jl")

testdatadir(x...) = joinpath("data", x...)
l = 7
n = 5

@testset "Network Motif I" begin
    raster_I = repeat_padded_motif("I", l, l, n, n)

    contributions_I = sequence_class_tricorr(raster_I, l, l)
    @test contributions_I[1] == sum(raster_I) / prod(size(raster_I)) && all(contributions_I[begin+1:end] .== 0)
end

@testset "Network Motif II" begin
    raster_II = repeat_padded_motif("II", l, l, n, n)

    contributions_II = sequence_class_tricorr(raster_II, l, l)
    @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
end

@testset "Network Motif III" begin
    raster_III = repeat_padded_motif("III", l, l, n, n)

    contributions_III = sequence_class_tricorr(raster_III, l, l)
    @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
end

@testset "Network Motif IV" begin
    raster_IV = repeat_padded_motif("IV", l, l, n, n)

    contributions_IV = sequence_class_tricorr(raster_IV, l, l)
    @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
end

@testset "Network Motif V" begin
    raster_V = repeat_padded_motif("V", l, l, n, n)

    contributions_V = sequence_class_tricorr(raster_V, l, l)
    @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
end

@testset "Network Motif VI" begin
    raster_VI = repeat_padded_motif("VI", l, l, n, n)

    contributions_VI = sequence_class_tricorr(raster_VI, l, l)
    @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
end

@testset "Network Motif VII" begin
    raster_VII = repeat_padded_motif("VII", l, l, n, n)

    contributions_VII = sequence_class_tricorr(raster_VII, l, l)
    @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
end

@testset "Network Motif VIII" begin
    raster_VIII = repeat_padded_motif("VIII", l, l, n, n)

    contributions_VIII = sequence_class_tricorr(raster_VIII, l, l)
    @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
end

@testset "Network Motif IX" begin
    raster_IX = repeat_padded_motif("IX", l, l, n, n)

    contributions_IX = sequence_class_tricorr(raster_IX, l, l)
    @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
end

@testset "Network Motif X" begin
    raster_X = repeat_padded_motif("X", l, l, n, n)

    contributions_X = sequence_class_tricorr(raster_X, l, l)
    @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
end

@testset "Network Motif XI" begin
    raster_XI = repeat_padded_motif("XI", l, l, n, n)

    contributions_XI = sequence_class_tricorr(raster_XI, l, l)
    @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
end

@testset "Network Motif XII" begin
    raster_XII = repeat_padded_motif("XII", l, l, n, n)

    contributions_XII = sequence_class_tricorr(raster_XII, l, l)
    @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
end

@testset "Network Motif XIII" begin
    raster_XIII = repeat_padded_motif("XIII", l, l, n, n)

    contributions_XIII = sequence_class_tricorr(raster_XIII, l, l)
    @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
end

@testset "Network Motif XIV" begin
    raster_XIV = repeat_padded_motif("XIV", l, l, n, n)
    contributions_XIV = @btime sequence_class_tricorr($raster_XIV, l, l)
    @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
end