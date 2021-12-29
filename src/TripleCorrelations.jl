module TripleCorrelations

using LoopVectorization, OffsetArrays, PaddedViews

include("triple_correlation.jl")
export TripleCorrelation, triple_correlation

include("sequence_class.jl")
export sequence_class_tricorr, sequence_class_tricorr_zeropad,
    sequence_class_tricorr!, sequence_class_tricorr_zeropad!,
    sequence_class_tricorr_unrolled

include("marginal.jl")
export time_tricorr_zeropad, space_tricorr_zeropad, space_time_tricorr_zeropad,
    marginal_tricorr_zeropad,
    time_tricorr_zeropad!, space_time_tricorr_zeropad!, space_tricorr_zeropad!,
    marginal_tricorr_zeropad!
end # module
