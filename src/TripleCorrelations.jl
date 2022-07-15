module TripleCorrelations

using LoopVectorization, OffsetArrays, PaddedViews, Base.Iterators
using CSV, Tables
using Memoize
using Random
using Base.Threads, ThreadsX
using Statistics, StatsBase

include("motif_examples.jl")
export rand_motif

include("boundaries.jl")
export AbstractBoundaryCondition, Periodic, ZeroPadded, PeriodicExtended

include("triple_correlation.jl")
export TripleCorrelation, triple_correlation

include("sequence_class.jl")
export sequence_class_tricorr, sequence_class_tricorr!

include("bootstrap.jl")
export bootstrap_sequence_classes, bootstrap_sequence_classes_nonzero,
    bootstrap_normed_sequence_classes

# FIXME switch marginals to new boundary system
include("marginal.jl")
export time_tricorr_zeropad, space_tricorr_zeropad, space_time_tricorr_zeropad,
    marginal_tricorr_zeropad,
    time_tricorr_zeropad!, space_time_tricorr_zeropad!, space_tricorr_zeropad!,
    marginal_tricorr_zeropad!

include("expectation.jl")
export expectation_of_independent_spiking_conditioned_on_rate, 
    sequence_classes_divide_E_given_constituents,
    sequence_classes_divide_E_given_rate,
    variance_of_standard_normals,
    estimate_std_of_standard_normals,
    estimate_μ, estimate_σ

include("snippet_processing.jl")
export extrema_01!, zscore!

include("contributions_processing.jl")
export precalculate

end # module
