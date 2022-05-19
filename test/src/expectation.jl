@quickactivate "test"

using Test, BenchmarkTools
using TripleCorrelations

using CairoMakie, Random

boundary = Periodic()#PeriodicExtended(9)
N_trials = 10
N_neurons = 50; N_times = 70
lag_extents = (20, 17)
firing_rate = 0.1
#raster = Matrix{Bool}(rand(N_neurons, N_times) .<= firing_rate)
raster = zeros(Bool, N_neurons, N_times)
raster[1:30] .= 1
bootstrapped = mapreduce(vcat, 1:N_trials) do _
    raster .= rand(N_neurons, N_times) .<= firing_rate
    expectation = expectation_of_independent_spiking_conditioned_on_rate(raster, boundary, lag_extents)
    sequence_class_tricorr(raster, boundary, lag_extents) .- expectation
end
motif_class = repeat(collect(0:13), outer=(N_trials,))
boxplot(motif_class, bootstrapped)

boundary = Periodic()#PeriodicExtended(9)
N_trials = 100
N_neurons = 50; N_times = 70
lag_extents = (20, 17)
firing_rate = 0.1
#raster = Matrix{Bool}(rand(N_neurons, N_times) .<= firing_rate)
raster = randn(N_neurons, N_times)
bootstrapped = mapreduce(vcat, 1:N_trials) do _
    randn!(raster)
    variance = variance_of_standard_normals(boundary, lag_extents)
    sequence_class_tricorr(raster, boundary, lag_extents) ./ sqrt.(variance)
end
motif_class = repeat(collect(0:13), outer=(N_trials,))
boxplot(motif_class, bootstrapped)