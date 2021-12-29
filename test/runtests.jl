using SafeTestsets

@time begin
    #@time @safetestset "CSV Motif Detection" begin include("csv_motif_detection.jl") end
    @time @safetestset "Example Motif Detection" begin include("example_motif_detection.jl") end

    #@time @safetestset "Benchmarks" begin include("benchmarks.jl") end
end