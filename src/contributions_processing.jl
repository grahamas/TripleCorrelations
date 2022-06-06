

function precalculate(func::Function, args...)
    func
end

function precalculate(::typeof(zscore!), assumption::IndStdNormal, condition, all_rasters, boundary, lag_extents)
    σ = estimate_σ(assumption, condition, example_raster, boundary, lag_extents)
    function zscore_stdnorm!(o, i, raster)
        zscore!(o, i, 0, σ)
    end
end
function precalculate(::typeof(zscore!), assumption::IndBernoulli, condition, all_rasters::AbstractArray, boundary, lag_extents)
    spike_counts = count.(all_rasters)
    raster_size = size(first(all_rasters))
    μ_from_count = Dict(count => estimate_μ(assumption, condition, count, raster_size, boundary, lag_extents) for count in unique(spike_counts))
    σ_from_count = Dict(count => estimate_σ(assumption, condition, count, raster_size, boundary, lag_extents) for count in unique(spike_counts))
    function zscore_bernoulli!(o, i, raster)
        μ = μ_from_count[count(raster)]
        σ = σ_from_count[count(raster)]
        zscore!(o, i, μ, σ)
    end
end

# prep: Function -> 