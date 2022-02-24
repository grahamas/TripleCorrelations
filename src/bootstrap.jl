
function bootstrap_normed_sequence_classes(raster, boundary, n_lag, t_lag; n_bootstraps, bootstraps_step=max(n_bootstraps÷5,1))
    raw_sequence_classes = sequence_class_tricorr(raster, boundary, n_lag, t_lag)
    raw_sequence_classes ./ bootstrap_sequence_classes_nonzero(raster, boundary, n_lag, t_lag, n_bootstraps, bootstraps_step)
end

# Don't memoize because raster is variable
function bootstrap_sequence_classes(raster::BitArray, boundary, n_lag::Int, t_lag::Int, n_bootstraps::Int)
    bootstrap_sequence_classes(boundary, size(raster)..., sum(raster), n_lag, t_lag,  n_bootstraps)
end

@memoize function bootstrap_sequence_classes(boundary, n::Int, t::Int, count_ones::Int, n_lag::Int, t_lag::Int, n_bootstraps::Int)
    unshuffled_raster = zeros(Int, n, t)
    unshuffled_raster[1:count_ones] .= 1
    bootstrap_sequence_classes!(unshuffled_raster, boundary, n_lag, t_lag, n_bootstraps)
end

function bootstrap_sequence_classes!(inplace_src_raster, boundary, n_lag::Int, t_lag::Int, n_bootstraps::Int)
    bootstrap_tricorr = sum(
        sequence_class_tricorr(shuffle!(inplace_src_raster), boundary, n_lag, t_lag) 
            for _ ∈ 1:n_bootstraps
    )
    bootstrap_tricorr ./= n_bootstraps
    return bootstrap_tricorr
end

# Don't memoize because raster is variable
function bootstrap_sequence_classes_nonzero(raster::BitArray, boundary, n_lag::Int, t_lag::Int, n_bootstraps::Int, bootstraps_step::Int)
    bootstrap_sequence_classes_nonzero(boundary, size(raster)..., count(raster), n_lag, t_lag, n_bootstraps, bootstraps_step)
end

@memoize function bootstrap_sequence_classes_nonzero(boundary, n::Int, t::Int, count_ones::Int, n_lag::Int, t_lag::Int, n_bootstraps::Int, bootstraps_step::Int)
    unshuffled_raster = zeros(Int, n, t)
    unshuffled_raster[1:count_ones] .= 1
    sequence_class_bootstrapped = bootstrap_sequence_classes!(unshuffled_raster, boundary, n_lag, t_lag, n_bootstraps) .* n_bootstraps
    while any(sequence_class_bootstrapped .== 0)
        @warn "insufficient n_bootstraps = $n_bootstraps [(n,t) = $((n,t)); lag = $((n_lag,t_lag))]"
        n_bootstraps += bootstraps_step
        sequence_class_bootstrapped += sum(
            sequence_class_tricorr(shuffle!(unshuffled_raster), boundary, n_lag, t_lag) 
                for _ ∈ 1:bootstraps_step
        )
    end
    sequence_class_bootstrapped ./= n_bootstraps
    return sequence_class_bootstrapped
end