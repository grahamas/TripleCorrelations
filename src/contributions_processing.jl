

function precalculate(func::Function, args...)
    func
end

function precalculate(::typeof(zscore!), assumption::IndStdNormal, condition, all_rasters, boundary, lag_extents)
    σ = estimate_σ(assumption, condition, first(all_rasters), boundary, lag_extents)
    function zscore_stdnorm!(o, i, raster)
        o .= i ./ σ
    end
end
function precalculate(::typeof(zscore!), assumption::IndBernoulli, condition, all_rasters::AbstractArray, boundary, lag_extents)
    spike_counts = count.(all_rasters)
    observed_spike_counts = unique(spike_counts)
    raster_size = size(first(all_rasters))
    min_spike_counts = max(0, minimum(observed_spike_counts)-3) #-3 if sig overlaps with noise
    max_spike_counts = min(prod(raster_size), maximum(observed_spike_counts)+3) #+3 if sig overlapped with noise in example
    possible_spike_counts = union(observed_spike_counts, Set(min_spike_counts:minimum(observed_spike_counts)), Set(maximum(observed_spike_counts):max_spike_counts))
    μ_from_count = Dict(count => estimate_μ(assumption, condition, count, raster_size, boundary, lag_extents) for count in possible_spike_counts)
    σ_from_count = Dict(count => estimate_σ(assumption, condition, count, raster_size, boundary, lag_extents) for count in possible_spike_counts)
    function zscore_bernoulli!(o, i, raster)
        μ = μ_from_count[count(raster)]
        σ = σ_from_count[count(raster)]
        @. o = (i - μ) / σ
    end
end

# prep: Function -> 