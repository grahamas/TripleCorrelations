using TripleCorrelations
using Test, BenchmarkTools

include("src/helpers.jl")

testdatadir(x...) = joinpath("data", x...)

@testset "Network Motif I" begin
    fn = testdatadir("I.csv")
    raster_I = load_raster(fn)

    contributions_I = sequence_class_tricorr(raster_I, 7, 7)
    @show contributions_I
    @test contributions_I[1] == sum(raster_I) / prod(size(raster_I)) && all(contributions_I[begin+1:end] .== 0)
end

@testset "Network Motif II" begin
    fn = testdatadir("II.csv")
    raster_II = load_raster(fn)

    contributions_II = sequence_class_tricorr(raster_II, 7, 7)
    @show contributions_II
    @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
end

@testset "Network Motif III" begin
    fn = testdatadir("III.csv")
    raster_III = load_raster(fn)

    contributions_III = sequence_class_tricorr(raster_III, 7, 7)
    @show contributions_III
    @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
end

@testset "Network Motif IV" begin
    fn = testdatadir("IV.csv")
    raster_IV = load_raster(fn)

    contributions_IV = sequence_class_tricorr(raster_IV, 7, 7)
    @show contributions_IV
    @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
end

@testset "Network Motif V" begin
    fn = testdatadir("V.csv")
    raster_V = load_raster(fn)

    contributions_V = sequence_class_tricorr(raster_V, 7, 7)
    @show contributions_V
    @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
end

@testset "Network Motif VI" begin
    fn = testdatadir("VI.csv")
    raster_VI = load_raster(fn)

    contributions_VI = sequence_class_tricorr(raster_VI, 7, 7)
    @show contributions_VI
    @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
end

@testset "Network Motif VII" begin
    fn = testdatadir("VII.csv")
    raster_VII = load_raster(fn)

    contributions_VII = sequence_class_tricorr(raster_VII, 7, 7)
    @show contributions_VII
    @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
end

@testset "Network Motif VIII" begin
    fn = testdatadir("VIII.csv")
    raster_VIII = load_raster(fn)

    contributions_VIII = sequence_class_tricorr(raster_VIII, 7, 7)
    @show contributions_VIII
    @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
end

@testset "Network Motif IX" begin
    fn = testdatadir("IX.csv")
    raster_IX = load_raster(fn)

    contributions_IX = sequence_class_tricorr(raster_IX, 7, 7)
    @show contributions_IX
    @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
end

@testset "Network Motif X" begin
    fn = testdatadir("X.csv")
    raster_X = load_raster(fn)

    contributions_X = sequence_class_tricorr(raster_X, 7, 7)
    @show contributions_X
    @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
end

@testset "Network Motif XI" begin
    fn = testdatadir("XI.csv")
    raster_XI = load_raster(fn)

    contributions_XI = sequence_class_tricorr(raster_XI, 7, 7)
    @show contributions_XI
    @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
end

@testset "Network Motif XII" begin
    fn = testdatadir("XII.csv")
    raster_XII = load_raster(fn)

    contributions_XII = sequence_class_tricorr(raster_XII, 7, 7)
    @show contributions_XII
    @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
end

@testset "Network Motif XIII" begin
    fn = testdatadir("XIII.csv")
    raster_XIII = load_raster(fn)

    contributions_XIII = sequence_class_tricorr(raster_XIII, 7, 7)
    @show contributions_XIII
    @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
end

@testset "Network Motif XIV" begin
    fn = testdatadir("XIV.csv")
    raster_XIV = load_raster(fn)
    contributions_XIV = @btime sequence_class_tricorr($raster_XIV, 7, 7)
    @show contributions_XIV
    @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
end