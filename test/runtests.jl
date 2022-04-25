using SafeTestsets

@time begin
    #@time @safetestset "CSV Motif Detection" begin include("src/csv_motif_detection.jl") end
    @time @safetestset "Example Motif Detection" begin include("src/example_motif_detection.jl") end
    @time @safetestset "3D Motif Detection" begin include("src/3D_detection.jl") end
    #@time @safetestset "Benchmarks" begin include("src/benchmarks.jl") end
end