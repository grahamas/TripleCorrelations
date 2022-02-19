# Don't memoize because raster is variable
function bootstrap_sequence_classes(raster::BitArray, n_lag::Int, t_lag::Int, n_bootstrap::Int)
    bootstrap_sequence_classes(size(raster)..., sum(raster), n_lag, t_lag,  n_bootstrap)
end

@memoize function bootstrap_sequence_classes(n::Int, t::Int, count_ones::Int, n_lag::Int, t_lag::Int, n_bootstrap::Int=10)
    unshuffled_raster = zeros(Int, n, t)
    unshuffled_raster[1:count_ones] .= 1
    bootstrap_sequence_classes!(unshuffled_raster, n_lag, t_lag, n_bootstrap)
end

function bootstrap_sequence_classes!(inplace_src_raster, n_lag::Int, t_lag::Int, n_bootstrap::Int=10)
    bootstrap_tricorr = sum(
        sequence_class_tricorr(shuffle!(inplace_src_raster), n_lag, t_lag) 
            for _ ∈ 1:n_bootstrap
    )
    bootstrap_tricorr ./= n_bootstrap
    return bootstrap_tricorr
end

# Don't memoize because raster is variable
function bootstrap_sequence_classes_nonzero(raster::BitArray, n_lag::Int, t_lag::Int, n_bootstrap::Int=10, bootstrap_step::Int=5)
    bootstrap_sequence_classes_nonzero(size(raster)..., count(raster), n_lag, t_lag, n_bootstrap, bootstrap_step)
end

@memoize function bootstrap_sequence_classes_nonzero(n::Int, t::Int, count_ones::Int, n_lag::Int, t_lag::Int, n_bootstrap::Int=10, bootstrap_step::Int=5)
    unshuffled_raster = zeros(Int, n, t)
    unshuffled_raster[1:count_ones] .= 1
    sequence_class_bootstrapped = bootstrap_sequence_classes!(unshuffled_raster, n_lag, t_lag, n_bootstrap) .* n_bootstrap
    while any(sequence_class_bootstrapped .== 0)
        @warn "insufficient n_bootstrap = $n_bootstrap [(n,t) = $((n,t)); lag = ($(n_lag,t_lag))]"
        n_bootstrap += bootstrap_step
        sequence_class_bootstrapped += sum(
            sequence_class_tricorr(shuffle!(unshuffled_raster), n_lag, t_lag) 
                for _ ∈ 1:bootstrap_step
        )
    end
    sequence_class_bootstrapped ./= n_bootstrap
    return sequence_class_bootstrapped
end