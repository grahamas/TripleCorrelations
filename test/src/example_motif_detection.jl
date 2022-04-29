using TripleCorrelations
using Test, BenchmarkTools

testdatadir(x...) = joinpath("data", x...)
l = 14
n = 5

@testset "Zeropadding" begin
    @testset "Sequence Motif 0" begin
        motif_idx = 1
        raster = TripleCorrelations.repeat_padded_motif("0", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @show normed_contributions
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] == sum(raster) / prod(size(raster)) && all(contributions[begin+motif_idx:end] .== 0)
        @test contributions[10] == 0
    end

    @testset "Sequence Motif I" begin
        motif_idx = 2
        raster = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @show normed_contributions
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        motif_idx = 3
        raster = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end


    @testset "Sequence Motif III" begin
        motif_idx = 4
        raster = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV TOP" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_top_motif("IV", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV pad >> lag" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_motif("IV", 30, 30, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif V TOP pad >> lag" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_top_motif("IV", 30, 30, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end


    @testset "Sequence Motif V" begin
        motif_idx = 6
        raster = TripleCorrelations.repeat_padded_motif("V", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VI" begin
        motif_idx = 7
        raster = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        motif_idx = 8
        raster = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        motif_idx = 9
        raster = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        motif_idx = 11
        raster = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        motif_idx = 10
        raster = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        motif_idx = 12
        raster = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        motif_idx = 13
        raster = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, ZeroPadded(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        motif_idx = 14
        raster = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)
        @info "Timing zeropad"
        contributions = @btime sequence_class_tricorr($raster, ZeroPadded(), (l, l))
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end
end

@testset "PeriodicExtended" begin
    @testset "Sequence Motif 0" begin
        motif_idx = 1
        raster = TripleCorrelations.repeat_padded_motif("0", l, l, n, n)
        boundary = PeriodicExtended(l)

        contributions = sequence_class_tricorr(raster, boundary, (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @show normed_contributions
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] == sum(raster[:,boundary.boundary:(end-boundary.boundary)]) / TripleCorrelations.calculate_scaling_factor(raster, boundary) && all(contributions[begin+motif_idx:end] .== 0)
        @test contributions[10] == 0
    end

    @testset "Sequence Motif I" begin
        motif_idx = 2
        raster = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @show normed_contributions
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        motif_idx = 3
        raster = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end


    @testset "Sequence Motif III" begin
        motif_idx = 4
        raster = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV TOP" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_top_motif("IV", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV pad >> lag" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_motif("IV", 30, 30, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV TOP pad >> lag" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_top_motif("IV", 30, 30, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end


    @testset "Sequence Motif V" begin
        motif_idx = 6
        raster = TripleCorrelations.repeat_padded_motif("V", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VI" begin
        motif_idx = 7
        raster = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        motif_idx = 8
        raster = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        motif_idx = 9
        raster = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        motif_idx = 11
        raster = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        motif_idx = 10
        raster = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        motif_idx = 12
        raster = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        motif_idx = 13
        raster = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, ZeroPadded(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        motif_idx = 14
        raster = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)
        @info "Timing PeriodicExtended..."
        contributions = @btime sequence_class_tricorr($raster, PeriodicExtended($l ÷ 2), ($l, $l))
        @show contributions
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end
end

@testset "Periodic" begin
    @testset "Sequence Motif 0" begin
        motif_idx = 1
        raster = TripleCorrelations.repeat_padded_motif("0", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[1] == sum(raster) / prod(size(raster)) && all(contributions[begin+motif_idx:end] .== 0)
        @test contributions[10] == 0
    end

    @testset "Sequence Motif I" begin
        motif_idx = 2
        raster = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        motif_idx = 3
        raster = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end


    @testset "Sequence Motif III" begin
        motif_idx = 4
        raster = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV TOP" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_top_motif("IV", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV pad >> lag" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_motif("IV", 30, 30, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IV TOP pad >> lag" begin
        motif_idx = 5
        raster = TripleCorrelations.repeat_padded_top_motif("IV", 30, 30, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @show contributions
        @show normed_contributions
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end


    @testset "Sequence Motif V" begin
        motif_idx = 6
        raster = TripleCorrelations.repeat_padded_motif("V", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VI" begin
        motif_idx = 7
        raster = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        motif_idx = 8
        raster = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        motif_idx = 9
        raster = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        motif_idx = 11
        raster = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        motif_idx = 10
        raster = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        motif_idx = 12
        raster = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        motif_idx = 13
        raster = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions = sequence_class_tricorr(raster, Periodic(), (l, l))
        normed_contributions = rate_normed_sequence_classes(raster, Periodic(), (l, l))
        @test normed_contributions[motif_idx] >= 0 && all(normed_contributions[begin+motif_idx:end] .== -1)
        @test contributions[motif_idx] > 0 && all(contributions[begin+motif_idx:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        motif_idx = 14
        raster = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)
        @info "Timing periodic..."
        contributions = @btime sequence_class_tricorr($raster, Periodic(), (l, l))
        @test contributions[motif_idx] > 0
    end
end