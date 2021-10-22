module TripleCorrelations

using LoopVectorization, OffsetArrays

include("triple_correlation.jl")
export TripleCorrelation, triple_correlation

include("sequence_class.jl")
export sequence_class_tricorr

end # module
