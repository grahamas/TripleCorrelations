using TripleCorrelations
using Test, BenchmarkTools

testdatadir(x...) = joinpath("data", x...)
l = 7
n = 5

# @testset "Ambiguous padding" begin
#     @testset "Sequence Motif I" begin
#         raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

#         contributions_I = sequence_class_tricorr(raster_I, l, l)
#         @test contributions_I[1] == sum(raster_I) / prod(size(raster_I)) && all(contributions_I[begin+1:end] .== 0)
#     end

#     @testset "Sequence Motif II" begin
#         raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

#         contributions_II = sequence_class_tricorr(raster_II, l, l)
#         @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
#     end

#     @testset "Sequence Motif III" begin
#         raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

#         contributions_III = sequence_class_tricorr(raster_III, l, l)
#         @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
#     end


#     @testset "Sequence Motif IV" begin
#         raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

#         contributions_IV = sequence_class_tricorr(raster_IV, l, l)
#         @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
#     end

#     @testset "Sequence Motif V" begin
#         raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, 1, 1)

#         contributions_V = sequence_class_tricorr(raster_V, l, l)
#         @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
#     end

#     @testset "Sequence Motif V TOP" begin
#         raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", l, l, 1, 1)

#         contributions_V_top = sequence_class_tricorr(raster_V_top, l, l)
#         @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
#     end

#     @testset "Sequence Motif V pad >> lag" begin
#         raster_V = TripleCorrelations.repeat_padded_motif("V", 30, 30, 1, 1)

#         contributions_V = sequence_class_tricorr(raster_V, l, l)
#         @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
#     end

#     @testset "Sequence Motif V TOP pad >> lag" begin
#         raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", 30, 30, 1, 1)

#         contributions_V_top = sequence_class_tricorr(raster_V_top, l, l)
#         @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
#     end


#     @testset "Sequence Motif VI" begin
#         raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

#         contributions_VI = sequence_class_tricorr(raster_VI, l, l)
#         @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
#     end

#     @testset "Sequence Motif VII" begin
#         raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

#         contributions_VII = sequence_class_tricorr(raster_VII, l, l)
#         @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
#     end

#     @testset "Sequence Motif VIII" begin
#         raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

#         contributions_VIII = sequence_class_tricorr(raster_VIII, l, l)
#         @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
#     end

#     @testset "Sequence Motif IX" begin
#         raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

#         contributions_IX = sequence_class_tricorr(raster_IX, l, l)
#         @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
#     end

#     @testset "Sequence Motif X" begin
#         raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

#         contributions_X = sequence_class_tricorr(raster_X, l, l)
#         @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
#     end

#     @testset "Sequence Motif XI" begin
#         raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

#         contributions_XI = sequence_class_tricorr(raster_XI, l, l)
#         @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
#     end

#     @testset "Sequence Motif XII" begin
#         raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

#         contributions_XII = sequence_class_tricorr(raster_XII, l, l)
#         @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
#     end

#     @testset "Sequence Motif XIII" begin
#         raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)

#         contributions_XIII = sequence_class_tricorr(raster_XIII, l, l)
#         @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
#     end

#     @testset "Sequence Motif XIV" begin
#         raster_XIV = TripleCorrelations.repeat_padded_motif("XIV", l, l, n, n)
#         @info "Timing default"
#         contributions_XIV = @btime sequence_class_tricorr($raster_XIV, ZeroPadded(), l, l)
#         @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
#     end
# end

@testset "Zeropadding" begin
    @testset "Sequence Motif I" begin
        raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions_I = sequence_class_tricorr(raster_I, ZeroPadded(), l, l)
        @test contributions_I[1] == sum(raster_I) / prod(size(raster_I)) && all(contributions_I[begin+1:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions_II = sequence_class_tricorr(raster_II, ZeroPadded(), l, l)
        @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
    end

    @testset "Sequence Motif III" begin
        raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions_III = sequence_class_tricorr(raster_III, ZeroPadded(), l, l)
        @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
    end


    @testset "Sequence Motif IV" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions_IV = sequence_class_tricorr(raster_IV, ZeroPadded(), l, l)
        @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
    end

    @testset "Sequence Motif V" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, 1, 1)

        contributions_V = sequence_class_tricorr(raster_V, ZeroPadded(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", l, l, 1, 1)

        contributions_V_top = sequence_class_tricorr(raster_V_top, ZeroPadded(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V pad >> lag" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", 30, 30, 1, 1)

        contributions_V = sequence_class_tricorr(raster_V, ZeroPadded(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP pad >> lag" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", 30, 30, 1, 1)

        contributions_V_top = sequence_class_tricorr(raster_V_top, ZeroPadded(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end


    @testset "Sequence Motif VI" begin
        raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions_VI = sequence_class_tricorr(raster_VI, ZeroPadded(), l, l)
        @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions_VII = sequence_class_tricorr(raster_VII, ZeroPadded(), l, l)
        @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions_VIII = sequence_class_tricorr(raster_VIII, ZeroPadded(), l, l)
        @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions_IX = sequence_class_tricorr(raster_IX, ZeroPadded(), l, l)
        @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions_X = sequence_class_tricorr(raster_X, ZeroPadded(), l, l)
        @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions_XI = sequence_class_tricorr(raster_XI, ZeroPadded(), l, l)
        @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions_XII = sequence_class_tricorr(raster_XII, ZeroPadded(), l, l)
        @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)

        contributions_XIII = sequence_class_tricorr(raster_XIII, ZeroPadded(), l, l)
        @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
    end

    @testset "Sequence Motif XIV" begin
        raster_XIV = TripleCorrelations.repeat_padded_motif("XIV", l, l, n, n)
        @info "Timing zeropad"
        contributions_XIV = @btime sequence_class_tricorr($raster_XIV, ZeroPadded(), l, l)
        @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
    end
end

@testset "Zeropadding, unrolled" begin
    @testset "Sequence Motif I" begin
        raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions_I = sequence_class_tricorr_unrolled(raster_I, ZeroPadded(), l, l)
        @test contributions_I[1] == sum(raster_I) / prod(size(raster_I)) && all(contributions_I[begin+1:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions_II = sequence_class_tricorr_unrolled(raster_II, ZeroPadded(), l, l)
        @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
    end

    @testset "Sequence Motif III" begin
        raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions_III = sequence_class_tricorr_unrolled(raster_III, ZeroPadded(), l, l)
        @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
    end


    @testset "Sequence Motif IV" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions_IV = sequence_class_tricorr_unrolled(raster_IV, ZeroPadded(), l, l)
        @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
    end

    @testset "Sequence Motif V" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, 1, 1)

        contributions_V = sequence_class_tricorr_unrolled(raster_V, ZeroPadded(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", l, l, 1, 1)

        contributions_V_top = sequence_class_tricorr_unrolled(raster_V_top, ZeroPadded(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V pad >> lag" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", 30, 30, 1, 1)

        contributions_V = sequence_class_tricorr_unrolled(raster_V, ZeroPadded(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP pad >> lag" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", 30, 30, 1, 1)

        contributions_V_top = sequence_class_tricorr_unrolled(raster_V_top, ZeroPadded(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end


    @testset "Sequence Motif VI" begin
        raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions_VI = sequence_class_tricorr_unrolled(raster_VI, ZeroPadded(), l, l)
        @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions_VII = sequence_class_tricorr_unrolled(raster_VII, ZeroPadded(), l, l)
        @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions_VIII = sequence_class_tricorr_unrolled(raster_VIII, ZeroPadded(), l, l)
        @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions_IX = sequence_class_tricorr_unrolled(raster_IX, ZeroPadded(), l, l)
        @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions_X = sequence_class_tricorr_unrolled(raster_X, ZeroPadded(), l, l)
        @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions_XI = sequence_class_tricorr_unrolled(raster_XI, ZeroPadded(), l, l)
        @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions_XII = sequence_class_tricorr_unrolled(raster_XII, ZeroPadded(), l, l)
        @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)

        contributions_XIII = sequence_class_tricorr_unrolled(raster_XIII, ZeroPadded(), l, l)
        @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
    end

    @testset "Sequence Motif XIV" begin
        raster_XIV = TripleCorrelations.repeat_padded_motif("XIV", l, l, n, n)
        @info "Timing zeropad unrolled"
        contributions_XIV = @btime sequence_class_tricorr_unrolled($raster_XIV, ZeroPadded(), l, l)
        @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
    end
end

@testset "Periodic" begin
    @testset "Sequence Motif I" begin
        raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions_I = sequence_class_tricorr(raster_I, Periodic(), l, l)
        @test contributions_I[1] == sum(raster_I) / prod(size(raster_I)) && all(contributions_I[begin+1:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions_II = sequence_class_tricorr(raster_II, Periodic(), l, l)
        @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
    end

    @testset "Sequence Motif III" begin
        raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions_III = sequence_class_tricorr(raster_III, Periodic(), l, l)
        @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
    end


    @testset "Sequence Motif IV" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions_IV = sequence_class_tricorr(raster_IV, Periodic(), l, l)
        @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
    end

    @testset "Sequence Motif V" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, 1, 1)

        contributions_V = sequence_class_tricorr(raster_V, Periodic(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", l, l, 1, 1)

        contributions_V_top = sequence_class_tricorr(raster_V_top, Periodic(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V pad >> lag" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", 30, 30, 1, 1)

        contributions_V = sequence_class_tricorr(raster_V, Periodic(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP pad >> lag" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", 30, 30, 1, 1)

        contributions_V_top = sequence_class_tricorr(raster_V_top, Periodic(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end


    @testset "Sequence Motif VI" begin
        raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions_VI = sequence_class_tricorr(raster_VI, Periodic(), l, l)
        @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions_VII = sequence_class_tricorr(raster_VII, Periodic(), l, l)
        @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions_VIII = sequence_class_tricorr(raster_VIII, Periodic(), l, l)
        @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions_IX = sequence_class_tricorr(raster_IX, Periodic(), l, l)
        @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions_X = sequence_class_tricorr(raster_X, Periodic(), l, l)
        @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions_XI = sequence_class_tricorr(raster_XI, Periodic(), l, l)
        @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions_XII = sequence_class_tricorr(raster_XII, Periodic(), l, l)
        @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)

        contributions_XIII = sequence_class_tricorr(raster_XIII, Periodic(), l, l)
        @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
    end

    @testset "Sequence Motif XIV" begin
        raster_XIV = TripleCorrelations.repeat_padded_motif("XIV", l, l, n, n)
        @info "Timing periodic..."
        contributions_XIV = @btime sequence_class_tricorr($raster_XIV, Periodic(), l, l)
        @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
    end
end


@testset "Periodic, unrolled" begin
    @testset "Sequence Motif I" begin
        raster_I = TripleCorrelations.repeat_padded_motif("I", l, l, n, n)

        contributions_I = sequence_class_tricorr_unrolled(raster_I, Periodic(), l, l)
        @test contributions_I[1] == sum(raster_I) / prod(size(raster_I)) && all(contributions_I[begin+1:end] .== 0)
    end

    @testset "Sequence Motif II" begin
        raster_II = TripleCorrelations.repeat_padded_motif("II", l, l, n, n)

        contributions_II = sequence_class_tricorr_unrolled(raster_II, Periodic(), l, l)
        @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
    end

    @testset "Sequence Motif III" begin
        raster_III = TripleCorrelations.repeat_padded_motif("III", l, l, n, n)

        contributions_III = sequence_class_tricorr_unrolled(raster_III, Periodic(), l, l)
        @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
    end


    @testset "Sequence Motif IV" begin
        raster_IV = TripleCorrelations.repeat_padded_motif("IV", l, l, n, n)

        contributions_IV = sequence_class_tricorr_unrolled(raster_IV, Periodic(), l, l)
        @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
    end

    @testset "Sequence Motif V" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", l, l, 1, 1)

        contributions_V = sequence_class_tricorr_unrolled(raster_V, Periodic(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", l, l, 1, 1)

        contributions_V_top = sequence_class_tricorr_unrolled(raster_V_top, Periodic(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V pad >> lag" begin
        raster_V = TripleCorrelations.repeat_padded_motif("V", 30, 30, 1, 1)

        contributions_V = sequence_class_tricorr_unrolled(raster_V, Periodic(), l, l)
        @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
    end

    @testset "Sequence Motif V TOP pad >> lag" begin
        raster_V_top = TripleCorrelations.repeat_padded_top_motif("V", 30, 30, 1, 1)

        contributions_V_top = sequence_class_tricorr_unrolled(raster_V_top, Periodic(), l, l)
        @test contributions_V_top[5] > 0 && all(contributions_V_top[begin+5:end] .== 0)
    end


    @testset "Sequence Motif VI" begin
        raster_VI = TripleCorrelations.repeat_padded_motif("VI", l, l, n, n)

        contributions_VI = sequence_class_tricorr_unrolled(raster_VI, Periodic(), l, l)
        @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
    end

    @testset "Sequence Motif VII" begin
        raster_VII = TripleCorrelations.repeat_padded_motif("VII", l, l, n, n)

        contributions_VII = sequence_class_tricorr_unrolled(raster_VII, Periodic(), l, l)
        @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
    end

    @testset "Sequence Motif VIII" begin
        raster_VIII = TripleCorrelations.repeat_padded_motif("VIII", l, l, n, n)

        contributions_VIII = sequence_class_tricorr_unrolled(raster_VIII, Periodic(), l, l)
        @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
    end

    @testset "Sequence Motif IX" begin
        raster_IX = TripleCorrelations.repeat_padded_motif("IX", l, l, n, n)

        contributions_IX = sequence_class_tricorr_unrolled(raster_IX, Periodic(), l, l)
        @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
    end

    @testset "Sequence Motif X" begin
        raster_X = TripleCorrelations.repeat_padded_motif("X", l, l, n, n)

        contributions_X = sequence_class_tricorr_unrolled(raster_X, Periodic(), l, l)
        @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
    end

    @testset "Sequence Motif XI" begin
        raster_XI = TripleCorrelations.repeat_padded_motif("XI", l, l, n, n)

        contributions_XI = sequence_class_tricorr_unrolled(raster_XI, Periodic(), l, l)
        @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
    end

    @testset "Sequence Motif XII" begin
        raster_XII = TripleCorrelations.repeat_padded_motif("XII", l, l, n, n)

        contributions_XII = sequence_class_tricorr_unrolled(raster_XII, Periodic(), l, l)
        @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
    end

    @testset "Sequence Motif XIII" begin
        raster_XIII = TripleCorrelations.repeat_padded_motif("XIII", l, l, n, n)

        contributions_XIII = sequence_class_tricorr_unrolled(raster_XIII, Periodic(), l, l)
        @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
    end

    @testset "Sequence Motif XIV" begin
        raster_XIV = TripleCorrelations.repeat_padded_motif("XIV", l, l, n, n)
        @info "Timing periodic, unrolled..."
        contributions_XIV = @btime sequence_class_tricorr_unrolled($raster_XIV, Periodic(), l, l)
        @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
    end
end