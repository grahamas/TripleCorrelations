module TripleCorrelations

using LoopVectorization, OffsetArrays, PaddedViews

include("triple_correlation.jl")
export TripleCorrelation, triple_correlation

include("sequence_class.jl")
export sequence_class_tricorr, sequence_class_tricorr_zeropad

end # module
