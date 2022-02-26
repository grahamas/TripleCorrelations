module TripleCorrelations

using LoopVectorization, OffsetArrays, PaddedViews, IterTools
using CSV, Tables
using Memoize
using Random

include("motif_examples.jl")
export rand_motif

include("boundaries.jl")
export AbstractBoundaryCondition, Periodic, ZeroPadded

include("triple_correlation.jl")
export TripleCorrelation, triple_correlation

include("sequence_class.jl")
export sequence_class_tricorr, sequence_class_tricorr!
export sequence_class_tricorr_unrolled,
    sequence_class_tricorr_unrolled

include("bootstrap.jl")
export bootstrap_sequence_classes, bootstrap_sequence_classes_nonzero,
    bootstrap_normed_sequence_classes

# FIXME switch marginals to new boundary system
include("marginal.jl")
export time_tricorr_zeropad, space_tricorr_zeropad, space_time_tricorr_zeropad,
    marginal_tricorr_zeropad,
    time_tricorr_zeropad!, space_time_tricorr_zeropad!, space_tricorr_zeropad!,
    marginal_tricorr_zeropad!

end # module
