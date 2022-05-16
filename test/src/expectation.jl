using Test, BenchmarkTools
using TripleCorrelations

using CairoMakie, Random

# @testset "Test expectation for Periodic" begin
    N_trials = 1000
    N_neurons = 50; N_times = 27
    lag_neuron = 6; lag_time = 7
    firing_rate = 0.1
    raster = zeros(Bool, N_neurons, N_times)
    raster[1:round(Int, N_neurons*N_times*firing_rate)] .= 1
    expectation = expectation_conditioned_on_spike_count(raster, Periodic(), (lag_neuron, lag_time))
    bootstrapped = mapreduce(vcat, 1:N_trials) do _
        sequence_class_tricorr(shuffle!(raster), Periodic(), (lag_neuron, lag_time)) .- expectation
    end
    motif_class = repeat(collect(0:13), outer=(N_trials,))
    boxplot(motif_class, bootstrapped)
# end

