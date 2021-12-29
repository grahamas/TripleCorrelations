using CSV, OffsetArrays, Tables

function load_raster(fn)
    arr = CSV.read(fn, BitArray ∘ Tables.matrix, delim=",", header=false)
    return OffsetArray(arr, 1:size(arr,1), 0:(size(arr,2)-1))
end

function generate_random_raster(raster_size::Tuple, spike_prob=0.1)
    putative_inputs = rand(Float64, raster_size)
    BitArray(putative_inputs .<= spike_prob)
end