# FIXME revert to _unrolled
function bootstrap_normed_sequence_classes(raster, boundary, max_lags; n_bootstraps, bootstraps_step=max(n_bootstraps÷5,1))
    raw_sequence_classes = sequence_class_tricorr(raster, boundary, max_lags)
    if n_bootstraps > 0
        raw_sequence_classes ./ bootstrap_sequence_classes_nonzero(raster, boundary, max_lags, n_bootstraps, bootstraps_step)
    else
        raw_sequence_classes
    end
end

# Don't memoize because raster is variable
function bootstrap_sequence_classes(raster::Array{Bool}, boundary, max_lags, n_bootstraps::Int)
    bootstrap_sequence_classes(boundary, size(raster)..., sum(raster), max_lags,  n_bootstraps)
end
function bootstrap_sequence_classes(raster::Array{Bool}, boundary::PeriodicExtended, max_lags, n_bootstraps::Int)
    @error "Unsupported--- needs extension shuffling."
    bootstrap_sequence_classes(boundary, size(raster)..., sum(raster), max_lags,  n_bootstraps)
end

@memoize function bootstrap_sequence_classes(boundary, n::Int, t::Int, count_ones::Int, max_lags, n_bootstraps::Int)
    unshuffled_raster = zeros(Bool, n, t)
    unshuffled_raster[1:count_ones] .= 1
    bootstrap_sequence_classes!(unshuffled_raster, boundary, max_lags, n_bootstraps)
end


function bootstrap_sequence_classes!(inplace_src, boundary, max_lags, n_bootstraps::Int)
    bootstrap_tricorr = sum(
        sequence_class_tricorr(shuffle!(inplace_src), boundary, max_lags) 
            for _ ∈ 1:n_bootstraps
    )
    bootstrap_tricorr ./= n_bootstraps
    return bootstrap_tricorr
end
function bootstrap_sequence_classes!(inplace_src, boundary::PeriodicExtended, max_lags, n_bootstraps::Int)
    t0, t1 = boundary.t_bounds
    bootstrap_tricorr = mapreduce(+, 1:n_bootstraps) do _
        shuffle!(@view inplace_src[:, 1:(t0-1)])
        shuffle!(@view inplace_src[:, t0:t1])
        shuffle!(@view inplace_src[:, (t1+1):end])
        sequence_class_tricorr(inplace_src, boundary, max_lags)
    end
    bootstrap_tricorr ./= n_bootstraps
    return bootstrap_tricorr
end


# Don't memoize because raster is variable
function bootstrap_sequence_classes_nonzero(raster::Array{Bool}, boundary, max_lags, n_bootstraps::Int, bootstraps_step::Int)
    bootstrap_sequence_classes_nonzero(boundary, size(raster)..., count(raster), max_lags, n_bootstraps, bootstraps_step)
end


@memoize function bootstrap_sequence_classes_nonzero(boundary, n::Int, t::Int, count_ones::Int, max_lags, n_bootstraps::Int, bootstraps_step::Int)
    max_bootstraps = n_bootstraps * 20
    unshuffled_raster = zeros(Bool, n, t)
    unshuffled_raster[1:count_ones] .= 1
    sequence_class_bootstrapped = bootstrap_sequence_classes!(unshuffled_raster, boundary, max_lags, n_bootstraps) .* n_bootstraps
    while any(sequence_class_bootstrapped .== 0) && (n_bootstraps < max_bootstraps)
        @warn "insufficient n_bootstraps = $n_bootstraps [(n,t) = $((n,t)); lag = $((n_lag,t_lag))]"
        n_bootstraps += bootstraps_step
        sequence_class_bootstrapped += sum(
            sequence_class_tricorr(shuffle!(unshuffled_raster), boundary, max_lags) 
                for _ ∈ 1:bootstraps_step
        )
    end
    if n_bootstraps >= max_bootstraps
        @warn "BAD"
        @show boundary n t count_ones n_lag t_lag
    end
    sequence_class_bootstrapped ./= n_bootstraps
    return sequence_class_bootstrapped
end

@memoize function bootstrap_sequence_classes_nonzero(boundary::PeriodicExtended, n::Int, t::Int, count_ones::Int, max_lags, n_bootstraps::Int, bootstraps_step::Int)
    max_bootstraps = n_bootstraps * 20
    inplace_arr = zeros(Bool, n, t)
    inplace_arr[1:count_ones] .= 1
    sequence_class_bootstrapped = bootstrap_sequence_classes!(inplace_arr, boundary, max_lags, n_bootstraps) .* n_bootstraps
    t0, t1 = boundary.t_bounds
    while any(sequence_class_bootstrapped .== 0) && (n_bootstraps < max_bootstraps)
        @warn "insufficient n_bootstraps = $n_bootstraps [(n,t) = $((n,t)); lag = $((n_lag,t_lag))]"
        n_bootstraps += bootstraps_step
        
        sequence_class_bootstrapped += mapreduce(+, 1:bootstraps_step) do _
            shuffle!(@view inplace_arr[:, 1:(t0-1)])
            shuffle!(@view inplace_arr[:, t0:t1])
            shuffle!(@view inplace_arr[:, (t1+1):end])
            sequence_class_tricorr(inplace_arr, boundary, max_lags)
        end
    end
    if n_bootstraps >= max_bootstraps
        @warn "BAD"
        @show boundary n t count_ones n_lag t_lag
    end
    sequence_class_bootstrapped ./= n_bootstraps
    return sequence_class_bootstrapped
end

function bootstrap_sequence_classes_nonzero(arr::AbstractArray{<:AbstractFloat}, boundary::PeriodicExtended, max_lags, n_bootstraps::Int, bootstraps_step::Int)
    inplace_arr = copy(arr)
    max_bootstraps = n_bootstraps * 20
    sequence_class_bootstrapped = bootstrap_sequence_classes!(inplace_arr, boundary, max_lags, n_bootstraps) .* n_bootstraps
    t0, t1 = boundary.t_bounds
    while any(sequence_class_bootstrapped .== 0) && (n_bootstraps < max_bootstraps)
        @warn "insufficient n_bootstraps = $n_bootstraps [(n,t) = $((n,t)); lag = $((n_lag,t_lag))]"
        n_bootstraps += bootstraps_step
        sequence_class_bootstrapped += mapreduce(+, 1:bootstraps_step) do _
            shuffle!(@view inplace_arr[:, 1:(t0-1)])
            shuffle!(@view inplace_arr[:, t0:t1])
            shuffle!(@view inplace_arr[:, (t1+1):end])
            sequence_class_tricorr(inplace_arr, boundary, max_lags)
        end
    end
    if n_bootstraps >= max_bootstraps
        @warn "Degenerate array or something (all zeros)."
    end
    sequence_class_bootstrapped ./= n_bootstraps
    return sequence_class_bootstrapped
end

function bootstrap_sequence_classes_nonzero(arr::AbstractArray{<:AbstractFloat}, boundary::Union{Periodic,ZeroPadded}, max_lags, n_bootstraps::Int, bootstraps_step::Int)
    inplace_arr = copy(arr)
    max_bootstraps = n_bootstraps * 20
    sequence_class_bootstrapped = bootstrap_sequence_classes!(inplace_arr, boundary, max_lags, n_bootstraps) .* n_bootstraps
    while any(sequence_class_bootstrapped .== 0) && (n_bootstraps < max_bootstraps)
        @warn "insufficient n_bootstraps = $n_bootstraps [(n,t) = $((n,t)); lag = $((n_lag,t_lag))]"
        n_bootstraps += bootstraps_step
        sequence_class_bootstrapped += sum(
            sequence_class_tricorr(shuffle!(inplace_arr), boundary, max_lags) 
                for _ ∈ 1:bootstraps_step
        )
    end
    if n_bootstraps >= max_bootstraps
        @warn "Degenerate array or something (all zeros)."
    end
    sequence_class_bootstrapped ./= n_bootstraps
    return sequence_class_bootstrapped
end