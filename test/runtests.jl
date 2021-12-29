using SafeTestsets

@time begin
    #@time @safetestset "Motif Detection" begin include("motif_detection.jl") end
    @time @safetestset "Benchmarks" begin include("benchmarks.jl") end
end