@quickactivate "test"

using Test, BenchmarkTools
using TripleCorrelations

using CairoMakie, Random

N_trials = 100
N_neurons = 50; N_times = 70
lag_extents = (20, 17)
firing_rate = 0.1
raster = Matrix{Bool}(rand(N_neurons, N_times) .<= firing_rate)
bootstrapped = mapreduce(vcat, 1:N_trials) do _
    raster .= rand(N_neurons, N_times) .<= firing_rate
    expectation = expectation_of_independent_spiking_conditioned_on_rate(raster, Periodic(), lag_extents)
    sequence_class_tricorr(raster, Periodic(), lag_extents) .- expectation
end
motif_class = repeat(collect(0:13), outer=(N_trials,))
boxplot(motif_class, bootstrapped)

