module TripleCorrelations

using LoopVectorization, OffsetArrays, PaddedViews, Base.Iterators
using CSV, Tables
using Memoize
using Random
using Base.Threads, ThreadsX

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
export rate_normed_sequence_classes, expectation_conditioned_on_spike_count, expectation_conditioned_on_lower_orders

end # module
