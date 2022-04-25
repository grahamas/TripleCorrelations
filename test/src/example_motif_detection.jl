using TripleCorrelations
using Test, BenchmarkTools

testdatadir(x...) = joinpath("data", x...)
l = 14
n = 5

# @testset "Ambiguous padding" begin
#     @testset "Sequence Motif 0" begin
#         raster_0 = TripleCorrelations.repeat_padded_motif("0", l, l, n, n)

#         contributions_0 = sequence_class_tricorr(raster_0, (l, l))
#         @test contributions_0[1] == sum(raster_0) / prod(size(raster_0)) && all(contributions_0[begin+1:end] .== 0)
#         @test contributions_0[10] == 0
#     end

#     @testset "Sequence Motif I" begin
#         raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

#         contributions_I = sequence_class_tricorr(raster_I, (l, l))
#         @test contributions_I[2] > 0 && all(contributions_I[begin+2:end] .== 0)
#     end

#     @testset "Sequence Motif II" begin
#         raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

#         contributions_II = sequence_class_tricorr(raster_II, (l, l))
#         @test contributions_II[3] > 0 && all(contributions_II[begin+3:end] .== 0)
#     end


#     @testset "Sequence Motif III" begin
#         raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

#         contributions_III = sequence_class_tricorr(raster_III, (l, l))
#         @test contributions_III[4] > 0 && all(contributions_III[begin+4:end] .== 0)
#     end

#     @testset "Sequence Motif IV" begin
#         raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

#         contributions_IV = sequence_class_tricorr(raster_IV, (l, l))
#         @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
#     end

#     @testset "Sequence Motif V TOP" begin
#         raster_IV_top = TripleCorrelations.repeat_padded_top_motif("IV", l, l, n, n)

#         contributions_IV_top = sequence_class_tricorr(raster_IV_top, (l, l))
#         @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
#     end

#     @testset "Sequence Motif V pad >> lag" begin
#         raster_IV = TripleCorrelations.repeat_padded_motif("IV", 30, 30, n, n)

#         contributions_IV = sequence_class_tricorr(raster_IV, (l, l))
#         @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
#     end

#     @testset "Sequence Motif V TOP pad >> lag" begin
#         raster_IV_top = TripleCorrelations.repeat_padded_top_motif("IV", 30, 30, n, n)

#         contributions_IV_top = sequence_class_tricorr(raster_IV_top, (l, l))
#         @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
#     end


#     @testset "Sequence Motif V" begin
#         raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, n, n)

#         contributions_V = sequence_class_tricorr(raster_V, (l, l))
#         @test contributions_V[6] > 0 && all(contributions_V[begin+6:end] .== 0)
#     end

#     @testset "Sequence Motif VI" begin
#         raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

#         contributions_VI = sequence_class_tricorr(raster_VI, (l, l))
#         @test contributions_VI[7] > 0 && all(contributions_VI[begin+7:end] .== 0)
#     end

#     @testset "Sequence Motif VII" begin
#         raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

#         contributions_VII = sequence_class_tricorr(raster_VII, (l, l))
#         @test contributions_VII[8] > 0 && all(contributions_VII[begin+8:end] .== 0)
#     end

#     @testset "Sequence Motif VIII" begin
#         raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

#         contributions_VIII = sequence_class_tricorr(raster_VIII, (l, l))
#         @test contributions_VIII[9] > 0 && all(contributions_VIII[begin+9:end] .== 0)
#     end

#     @testset "Sequence Motif X" begin
#         raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

#         contributions_X = sequence_class_tricorr(raster_X, (l, l))
#         @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
#     end

#     @testset "Sequence Motif IX" begin
#         raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

#         contributions_IX = sequence_class_tricorr(raster_IX, (l, l))
#         @test contributions_IX[11] > 0 && all(contributions_IX[begin+11:end] .== 0)
#     end

#     @testset "Sequence Motif XI" begin
#         raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

#         contributions_XI = sequence_class_tricorr(raster_XI, (l, l))
#         @test contributions_XI[12] > 0 && all(contributions_XI[begin+12:end] .== 0)
#     end

#     @testset "Sequence Motif XII" begin
#         raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

#         contributions_XII = sequence_class_tricorr(raster_XII, (l, l))
#         @test contributions_XII[13] > 0 && all(contributions_XII[begin+13:end] .== 0)
#     end

#     @testset "Sequence Motif XIII" begin
#         raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)
#         @info "Timing default"
#         contributions_XIII = @btime sequence_class_tricorr($raster_XIII, ZeroPadded(), (l, l))
#         @test contributions_XIII[14] > 0 && all(contributions_XIII[begin+14:end] .== 0)
#     end
# end

@testset "Zeropadding" begin
    @testset "Sequence Motif 0" begin
        raster_0 = TripleCorrelations.repeat_padded_motif("0", l, l, n, n)

        contributions_0 = sequence_class_tricorr(raster_0, ZeroPadded(), (l, l))
        @test contributions_0[1] == sum(raster_0) / prod(size(raster_0)) && all(contributions_0[begin+1:end] .== 0)
        @test contributions_0[10] == 0
    end

    @testset "Sequence Motif I" begin
        raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions_I = sequence_class_tricorr(raster_I, ZeroPadded(), (l, l))
        @test contributions_I[2] > 0 && all(contributions_I[begin+2:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions_II = sequence_class_tricorr(raster_II, ZeroPadded(), (l, l))
        @test contributions_II[3] > 0 && all(contributions_II[begin+3:end] .== 0)
    end


    @testset "Sequence Motif III" begin
        raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions_III = sequence_class_tricorr(raster_III, ZeroPadded(), (l, l))
        @test contributions_III[4] > 0 && all(contributions_III[begin+4:end] .== 0)
    end

    @testset "Sequence Motif IV" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions_IV = sequence_class_tricorr(raster_IV, ZeroPadded(), (l, l))
        @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP" begin
        raster_IV_top = TripleCorrelations.repeat_padded_top_motif("IV", l, l, n, n)

        contributions_IV_top = sequence_class_tricorr(raster_IV_top, ZeroPadded(), (l, l))
        @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V pad >> lag" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", 30, 30, n, n)

        contributions_IV = sequence_class_tricorr(raster_IV, ZeroPadded(), (l, l))
        @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP pad >> lag" begin
        raster_IV_top = TripleCorrelations.repeat_padded_top_motif("IV", 30, 30, n, n)

        contributions_IV_top = sequence_class_tricorr(raster_IV_top, ZeroPadded(), (l, l))
        @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
    end


    @testset "Sequence Motif V" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, n, n)

        contributions_V = sequence_class_tricorr(raster_V, ZeroPadded(), (l, l))
        @test contributions_V[6] > 0 && all(contributions_V[begin+6:end] .== 0)
    end

    @testset "Sequence Motif VI" begin
        raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions_VI = sequence_class_tricorr(raster_VI, ZeroPadded(), (l, l))
        @test contributions_VI[7] > 0 && all(contributions_VI[begin+7:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions_VII = sequence_class_tricorr(raster_VII, ZeroPadded(), (l, l))
        @test contributions_VII[8] > 0 && all(contributions_VII[begin+8:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions_VIII = sequence_class_tricorr(raster_VIII, ZeroPadded(), (l, l))
        @test contributions_VIII[9] > 0 && all(contributions_VIII[begin+9:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions_X = sequence_class_tricorr(raster_X, ZeroPadded(), (l, l))
        @test contributions_X[11] > 0 && all(contributions_X[begin+11:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions_IX = sequence_class_tricorr(raster_IX, ZeroPadded(), (l, l))
        @test contributions_IX[10] > 0 && all(contributions_IX[begin+10:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions_XI = sequence_class_tricorr(raster_XI, ZeroPadded(), (l, l))
        @test contributions_XI[12] > 0 && all(contributions_XI[begin+12:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions_XII = sequence_class_tricorr(raster_XII, ZeroPadded(), (l, l))
        @test contributions_XII[13] > 0 && all(contributions_XII[begin+13:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)
        @info "Timing zeropad"
        contributions_XIII = @btime sequence_class_tricorr($raster_XIII, ZeroPadded(), (l, l))
        @test contributions_XIII[14] > 0 && all(contributions_XIII[begin+14:end] .== 0)
    end
end

@testset "PeriodicExtended" begin
    @testset "Sequence Motif 0" begin
        raster = TripleCorrelations.repeat_padded_motif("0", l, l, n, n)
        boundary = PeriodicExtended(l)

        contributions_0 = sequence_class_tricorr(raster, boundary, (l, l))
        @test contributions_0[1] == sum(raster[:,boundary.boundary:(end-boundary.boundary)]) / TripleCorrelations.calculate_scaling_factor(raster, boundary) && all(contributions_0[begin+1:end] .== 0)
        @test contributions_0[10] == 0
    end

    @testset "Sequence Motif I" begin
        raster = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions_I = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_I[2] > 0 && all(contributions_I[begin+2:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        raster = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions_II = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_II[3] > 0 && all(contributions_II[begin+3:end] .== 0)
    end


    @testset "Sequence Motif III" begin
        raster = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions_III = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_III[4] > 0 && all(contributions_III[begin+4:end] .== 0)
    end

    @testset "Sequence Motif IV" begin
        raster = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions_IV = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP" begin
        raster = TripleCorrelations.repeat_padded_top_motif("IV", l, l, n, n)

        contributions_IV_top = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V pad >> lag" begin
        raster = TripleCorrelations.repeat_padded_motif("IV", 30, 30, n, n)

        contributions_IV = sequence_class_tricorr(raster, PeriodicExtended((l+1, size(raster,2)-l)), (l, l))
        @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP pad >> lag" begin
        raster = TripleCorrelations.repeat_padded_top_motif("IV", 30, 30, n, n)

        contributions_IV_top = sequence_class_tricorr(raster, PeriodicExtended((l+1, size(raster,2)-l)), (l, l))
        @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
    end


    @testset "Sequence Motif V" begin
        raster = TripleCorrelations.repeat_padded_motif("V", l, l, n, n)

        contributions_V = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_V[6] > 0 && all(contributions_V[begin+6:end] .== 0)
    end

    @testset "Sequence Motif VI" begin
        raster = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions_VI = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_VI[7] > 0 && all(contributions_VI[begin+7:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        raster = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions_VII = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_VII[8] > 0 && all(contributions_VII[begin+8:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        raster = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions_VIII = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_VIII[9] > 0 && all(contributions_VIII[begin+9:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        raster = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions_X = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_X[11] > 0 && all(contributions_X[begin+11:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        raster = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions_IX = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_IX[10] > 0 && all(contributions_IX[begin+10:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        raster = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions_XI = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_XI[12] > 0 && all(contributions_XI[begin+12:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        raster = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions_XII = sequence_class_tricorr(raster, PeriodicExtended(l), (l, l))
        @test contributions_XII[13] > 0 && all(contributions_XII[begin+13:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        raster = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)
        @info "Timing PeriodicExtended..."
        contributions_XIII = @btime sequence_class_tricorr($raster, PeriodicExtended((1+$l, size($raster,2)-$l)), ($l, $l))
        @test contributions_XIII[14] > 0 && all(contributions_XIII[begin+14:end] .== 0)
    end
end

@testset "Periodic" begin
    @testset "Sequence Motif 0" begin
        raster_0 = TripleCorrelations.repeat_padded_motif("0", l, l, n, n)

        contributions_0 = sequence_class_tricorr(raster_0, Periodic(), (l, l))
        @test contributions_0[1] == sum(raster_0) / prod(size(raster_0)) && all(contributions_0[begin+1:end] .== 0)
        @test contributions_0[10] == 0
    end

    @testset "Sequence Motif I" begin
        raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions_I = sequence_class_tricorr(raster_I, Periodic(), (l, l))
        @test contributions_I[2] > 0 && all(contributions_I[begin+2:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions_II = sequence_class_tricorr(raster_II, Periodic(), (l, l))
        @test contributions_II[3] > 0 && all(contributions_II[begin+3:end] .== 0)
    end


    @testset "Sequence Motif III" begin
        raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions_III = sequence_class_tricorr(raster_III, Periodic(), (l, l))
        @test contributions_III[4] > 0 && all(contributions_III[begin+4:end] .== 0)
    end

    @testset "Sequence Motif IV" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions_IV = sequence_class_tricorr(raster_IV, Periodic(), (l, l))
        @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
    end

    @testset "Sequence Motif IV TOP" begin
        raster_IV_top = TripleCorrelations.repeat_padded_top_motif("IV", l, l, n, n)

        contributions_IV_top = sequence_class_tricorr(raster_IV_top, Periodic(), (l, l))
        @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
    end

    @testset "Sequence Motif IV pad >> lag" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", 30, 30, n, n)

        contributions_IV = sequence_class_tricorr(raster_IV, Periodic(), (l, l))
        @test contributions_IV[5] > 0 && all(contributions_IV[begin+5:end] .== 0)
    end

    @testset "Sequence Motif IV TOP pad >> lag" begin
        raster_IV_top = TripleCorrelations.repeat_padded_top_motif("IV", 30, 30, n, n)

        contributions_IV_top = sequence_class_tricorr(raster_IV_top, Periodic(), (l, l))
        @test contributions_IV_top[5] > 0 && all(contributions_IV_top[begin+5:end] .== 0)
    end


    @testset "Sequence Motif V" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, n, n)

        contributions_V = sequence_class_tricorr(raster_V, Periodic(), (l, l))
        @test contributions_V[6] > 0 && all(contributions_V[begin+6:end] .== 0)
    end

    @testset "Sequence Motif VI" begin
        raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions_VI = sequence_class_tricorr(raster_VI, Periodic(), (l, l))
        @test contributions_VI[7] > 0 && all(contributions_VI[begin+7:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions_VII = sequence_class_tricorr(raster_VII, Periodic(), (l, l))
        @test contributions_VII[8] > 0 && all(contributions_VII[begin+8:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions_VIII = sequence_class_tricorr(raster_VIII, Periodic(), (l, l))
        @test contributions_VIII[9] > 0 && all(contributions_VIII[begin+9:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions_X = sequence_class_tricorr(raster_X, Periodic(), (l, l))
        @test contributions_X[11] > 0 && all(contributions_X[begin+11:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions_IX = sequence_class_tricorr(raster_IX, Periodic(), (l, l))
        @test contributions_IX[10] > 0 && all(contributions_IX[begin+10:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions_XI = sequence_class_tricorr(raster_XI, Periodic(), (l, l))
        @test contributions_XI[12] > 0 && all(contributions_XI[begin+12:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions_XII = sequence_class_tricorr(raster_XII, Periodic(), (l, l))
        @test contributions_XII[13] > 0 && all(contributions_XII[begin+13:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)
        @info "Timing periodic..."
        contributions_XIII = @btime sequence_class_tricorr($raster_XIII, Periodic(), (l, l))
        @test contributions_XIII[14] > 0
    end
end