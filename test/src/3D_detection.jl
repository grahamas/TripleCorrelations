using TripleCorrelations
using Test, BenchmarkTools, MAT

testdatadir(x...) = joinpath("data", x...)

mat = matread(testdatadir("test_cases_3D_rasters_10x10x100.mat"))

lag_extents = (4, 4, 24)

@testset "Periodic tests in 3D (space x space x time)" begin

@testset "3D Periodic --- Motif 0" begin
    raster = mat["raster_3D_motif_0"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[1] == sum(raster) / prod(size(raster)) && all(contributions[begin+1:end] .== 0)
end

@testset "3D Periodic --- Motif I" begin
    raster = mat["raster_3D_motif_I"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[2] > 0 && all(contributions[begin+2:end] .== 0)
end

@testset "3D Periodic --- Motif II" begin
    raster = mat["raster_3D_motif_II"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[3] > 0 && all(contributions[begin+3:end] .== 0)
end

@testset "3D Periodic --- Motif III" begin
    raster = mat["raster_3D_motif_III"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[4] > 0 && all(contributions[begin+4:end] .== 0)
end

@testset "3D Periodic --- Motif IV" begin
    raster = mat["raster_3D_motif_IV"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[5] > 0 && all(contributions[begin+5:end] .== 0)
end

@testset "3D Periodic --- Motif V" begin
    raster = mat["raster_3D_motif_V"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[6] > 0 && all(contributions[begin+6:end] .== 0)
end

@testset "3D Periodic --- Motif VI" begin
    raster = mat["raster_3D_motif_VI"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[7] > 0 && all(contributions[begin+7:end] .== 0)
end

@testset "3D Periodic --- Motif VII" begin
    raster = mat["raster_3D_motif_VII"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[8] > 0 && all(contributions[begin+8:end] .== 0)
end

@testset "3D Periodic --- Motif VIII" begin
    raster = mat["raster_3D_motif_VIII"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[9] > 0 && all(contributions[begin+9:end] .== 0)
end

@testset "3D Periodic --- Motif IX" begin
    raster = mat["raster_3D_motif_IX"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[10] > 0 && all(contributions[begin+10:end] .== 0)
end

@testset "3D Periodic --- Motif X" begin
    raster = mat["raster_3D_motif_X"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[11] > 0 && all(contributions[begin+11:end] .== 0)
end

@testset "3D Periodic --- Motif XI" begin
    raster = mat["raster_3D_motif_XI"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[12] > 0 && all(contributions[begin+12:end] .== 0)
end

@testset "3D Periodic --- Motif XII" begin
    raster = mat["raster_3D_motif_XII"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[13] > 0 && all(contributions[begin+13:end] .== 0)
end

@testset "3D Periodic --- Motif XIII" begin
    raster = mat["raster_3D_motif_XIII"]

    contributions = sequence_class_tricorr(raster, Periodic(), lag_extents)
    @test contributions[14] > 0 && all(contributions[begin+14:end] .== 0)
end

end

@testset "ZeroPadded tests in 3D (space x space x time)" begin

    @testset "3D ZeroPadded --- Motif 0" begin
        raster = mat["raster_3D_motif_0"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[1] == sum(raster) / prod(size(raster)) && all(contributions[begin+1:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif I" begin
        raster = mat["raster_3D_motif_I"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[2] > 0 && all(contributions[begin+2:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif II" begin
        raster = mat["raster_3D_motif_II"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[3] > 0 && all(contributions[begin+3:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif III" begin
        raster = mat["raster_3D_motif_III"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[4] > 0 && all(contributions[begin+4:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif IV" begin
        raster = mat["raster_3D_motif_IV"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[5] > 0 && all(contributions[begin+5:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif V" begin
        raster = mat["raster_3D_motif_V"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[6] > 0 && all(contributions[begin+6:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif VI" begin
        raster = mat["raster_3D_motif_VI"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[7] > 0 && all(contributions[begin+7:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif VII" begin
        raster = mat["raster_3D_motif_VII"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[8] > 0 && all(contributions[begin+8:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif VIII" begin
        raster = mat["raster_3D_motif_VIII"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[9] > 0 && all(contributions[begin+9:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif IX" begin
        raster = mat["raster_3D_motif_IX"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[10] > 0 && all(contributions[begin+10:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif X" begin
        raster = mat["raster_3D_motif_X"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[11] > 0 && all(contributions[begin+11:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif XI" begin
        raster = mat["raster_3D_motif_XI"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[12] > 0 && all(contributions[begin+12:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif XII" begin
        raster = mat["raster_3D_motif_XII"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[13] > 0 && all(contributions[begin+13:end] .== 0)
    end
    
    @testset "3D ZeroPadded --- Motif XIII" begin
        raster = mat["raster_3D_motif_XIII"]
    
        contributions = sequence_class_tricorr(raster, ZeroPadded(), lag_extents)
        @test contributions[14] > 0 && all(contributions[begin+14:end] .== 0)
    end
    
    end