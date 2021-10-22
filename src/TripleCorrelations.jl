module TripleCorrelations

using LoopVectorization, OffsetArrays

include("triple_correlation.jl")
export TripleCorrelation, triple_correlation

end # module
