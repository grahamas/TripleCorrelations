# Truncating calculation

function lag_contribution(data::D, boundary::Periodic, p1::NTuple{N,Int}, p2::NTuple{N,Int}, data_λ₁=D(undef, size(data)), data_λ₂=D(undef, size(data))) where {N, T, D <: AbstractArray{T,N}}
    # Assume ns, ts < size(data)
    contribution = 0
    circshift!(data_λ₁, data, .-p1)
    circshift!(data_λ₂, data, .-p2)

    @tturbo for p0 ∈ CartesianIndices(data)
        contribution += data[p0] * data_λ₁[p0] * data_λ₂[p0]
    end
    return contribution
end

function lag_contribution_pre_shifted(data::D, boundary::Periodic, data_λ₁, data_λ₂) where {N, T, D <: AbstractArray{T,N}}
    # Assume ns, ts < size(data)
    contribution = 0

    @tturbo for p0 ∈ CartesianIndices(data)
        contribution += data[p0] * data_λ₁[p0] * data_λ₂[p0]
    end
    return contribution
end

function lag_contribution_pre_shifted(data::D, boundary::PeriodicExtended, data_λ₁, data_λ₂) where {N, T, D <: AbstractArray{T}}
    # Periodic in n; extended in t
    # FIXME should validate extension holds lags
    t_start, t_end = boundary.t_bounds
    contribution = 0

    data_view = view_slice_last(data, t_start:t_end)
    data_λ₁_view = view_slice_last(data_λ₁, t_start:t_end)
    data_λ₂_view = view_slice_last(data_λ₂, t_start:t_end)

    @tturbo for p ∈ CartesianIndices(data_view)
        contribution += data_view[p] * data_λ₁_view[p] * data_λ₂_view[p]
    end
    return contribution
end

# if n > 0: n; else 0
function if_positive(n)
    if n > 0
        n
    else
        0
    end
end

function lag_contribution(data::D, boundary::ZeroPadded, λ₁::NTuple{N}, λ₂::NTuple{N}) where {T,N, D <: AbstractArray{T,N}}
    λ_starts = map(λ₁, λ₂) do p1, p2
        max(1 - min(p1, p2), 1)
    end
    λ_stops = map(λ₁, λ₂, size(data)) do p1, p2, dim_max
        min(dim_max - max(p1, p2), dim_max)
    end

    base_ranges = map(λ_starts, λ_stops) do start, stop
        UnitRange(start, stop)
    end
    λ₁_ranges = map(base_ranges, λ₁) do range, p
        range .+ p
    end
    λ₂_ranges = map(base_ranges, λ₂) do range, p
        range .+ p
    end

    A = @view data[base_ranges...]
    A_λ₁ = @view data[λ₁_ranges...]
    A_λ₂ = @view data[λ₂_ranges...]

    contribution = 0
    @tturbo for p ∈ CartesianIndices(A)
        contribution += A[p] * A_λ₁[p] * A_λ₂[p]
    end
    return contribution
end

@inline function view_slice_last(arr::AbstractArray{T,N}, dx) where {T,N}
    view(arr, ntuple(_ -> Colon(), N - 1)..., dx)
end

function lag_contribution(data::D, boundary::PeriodicExtended, λ₁::NTuple{N,Int}, λ₂::NTuple{N,Int}, data_λ₁=D(undef, size(data)), data_λ₂=D(undef, size(data))) where {N, T, D <: AbstractArray{T}}
    # Periodic in n; extended in t
    # FIXME should validate extension holds lags
    t_start, t_end = boundary.t_bounds
    contribution = 0
    
    circshift!(data_λ₁, data, .-λ₁)
    circshift!(data_λ₂, data, .-λ₂)

    data_view = view_slice_last(data, t_start:t_end)
    data_λ₁_view = view_slice_last(data_λ₁, t_start:t_end)
    data_λ₂_view = view_slice_last(data_λ₂, t_start:t_end)

    @tturbo for p ∈ CartesianIndices(data_view)
        contribution += data_view[p] * data_λ₁_view[p] * data_λ₂_view[p]
    end
    return contribution
end

# FIXME 
# instead of filtering the element, track the el_idx,
# then split inner loop to go begin:el_idx-1 and el_idx+1:end
function filter_element(arr, el)
    filter(≠(el), arr)
end

function lag_sequence_class_contributions(tricorr::TripleCorrelation)
    contributions = zeros(14)
    cartesian_indices = CartesianIndices(tricorr.arr)
    for i ∈ 1:length(cartesian_indices)
        idx = cartesian_indices[i]
        class = TripleCorrelations.lag_motif_sequence_class(idx.I)
        contributions[class] += tricorr.arr[i]
    end
    return contributions
end

# # FIXME is this zeropadding?
# function sequence_class_tricorr!(class_contribution::AbstractVector, src, boundary:: lag_extents, lags_classifier::Function)
#     src = parent(src)

#     space_lag_range = -(space_max_lag):(space_max_lag)        
#     time_lag_range = -(time_max_lag):(time_max_lag)

#     (N_n, N_t) = size(src)

#     class_contribution .= 0
#     for n1 ∈ space_lag_range, n2 ∈ space_lag_range, 
#         t1 ∈ time_lag_range, t2 ∈ time_lag_range
#         class = lags_classifier(n1, n2, t1, t2)
#         n_start = max(1 - min(n1, n2), 1); t_start = max(1 - min(t1, t2), 1)
#         n_end = min(N_n - max(n1,n2), N_n)
#         t_end = min(N_t - max(t1,t2), N_t)
        
#         contribution = 0
#         # tturbo missing
#         for n ∈ n_start:n_end, t ∈ t_start:t_end
#             contribution += src[n, t] * src[n+n1,t+t1] * src[n+n2,t+t2]
#         end
#         class_contribution[class] += contribution
#     end
#     class_contribution ./= calculate_scaling_factor_zeropad(src)
# end

function sequence_class_tricorr(src, boundary::AbstractBoundaryCondition, lag_extents)
    N_network_classifications = 14
    network_class_contributions = Array{Float64}(undef, N_network_classifications)
    lags_classifier = lag_motif_sequence_class

    sequence_class_tricorr!(network_class_contributions, src, boundary, lag_extents, lags_classifier)
end

function sequence_class_tricorr!(class_contribution::AbstractVector, src::AbstractArray{T}, boundary::AbstractBoundaryCondition, lag_extents, lags_classifier::Function) where T
    src = parent(src)

    lag_ranges = map(((start,stop),) -> UnitRange(start,stop), zip(.-floor.(Int, lag_extents ./ 2), ceil.(Int, lag_extents ./ 2)))

    class_contribution .= 0
    for λ₁ ∈ Iterators.product(lag_ranges...), λ₂ ∈ Iterators.product(lag_ranges...)
        class = lags_classifier(λ₁, λ₂) #FIXME function call
        contribution = lag_contribution(src, boundary, λ₁, λ₂)
        class_contribution[class] += contribution
    end
    class_contribution ./= calculate_scaling_factor(src, boundary)
end

function sequence_class_tricorr(src::OffsetArray, args...)
    sequence_class_tricorr(parent(src), args...)
end

# Periodic calculation

function sequence_class_tricorr!(class_contribution::AbstractVector, src::SRC, boundary::Periodic, lag_extents, lags_classifier::Function) where {T, SRC<:AbstractArray{T}}
    src = parent(src)
    lag1_cache = typeof(src)(undef, size(src))
    lag2_cache = typeof(src)(undef, size(src))

    lag_ranges = map(((start,stop),) -> UnitRange(start,stop), zip(.-floor.(Int, lag_extents ./ 2), ceil.(Int, lag_extents ./ 2)))

    class_contribution .= 0
    for λ₁ ∈ Iterators.product(lag_ranges...)
        circshift!(lag1_cache, src, .-λ₁)
        for λ₂ ∈ Iterators.product(lag_ranges...)
            circshift!(lag2_cache, src, .-λ₂)
            class = lags_classifier(λ₁, λ₂)
            contribution = lag_contribution_pre_shifted(src, boundary, lag1_cache, lag2_cache)
            class_contribution[class] += contribution
        end
    end
    class_contribution ./= calculate_scaling_factor(src, boundary)
end

function sequence_class_tricorr!(class_contribution::AbstractVector, src::SRC, boundary::PeriodicExtended, lag_extents, lags_classifier::Function) where {T, SRC<:AbstractArray{T}}
    src = parent(src)
    lag1_cache = typeof(src)(undef, size(src))
    lag2_cache = typeof(src)(undef, size(src))

    periodic_lag_extents = lag_extents[1:end-1]
    periodic_lag_ranges = map(((start,stop),) -> UnitRange(start,stop), zip(.-floor.(Int, lag_extents ./ 2), ceil.(Int, periodic_lag_extents ./ 2)))

    extended_dim_lag_range = UnitRange(-floor(Int, lag_extents[end] / 2), ceil(Int, lag_extents[end] / 2))

    class_contribution .= 0
    for λ₁_periodic ∈ Iterators.product(periodic_lag_ranges...)
        circshift!(lag1_cache, src, (.-λ₁_periodic..., 0))
        for λ₂_periodic ∈ Iterators.product(periodic_lag_ranges...)
            circshift!(lag2_cache, src, (.-λ₂_periodic..., 0))
            for t₁ ∈ extended_dim_lag_range, t₂ ∈ extended_dim_lag_range
                λ₁ = (λ₁_periodic..., t₁)
                λ₂ = (λ₂_periodic..., t₂)
                class = lags_classifier(λ₁, λ₂)
                contribution = lag_contribution_pre_shifted(src, boundary, lag1_cache, lag2_cache)
                class_contribution[class] += contribution
            end
        end
    end
    class_contribution ./= calculate_scaling_factor(src, boundary)
end


allzeros(t) = mapreduce(iszero, &, t)
neuron_lags(t::NTuple{N}) where N = ntuple(i -> t[i], N-1)
time_lag(t::Tuple) = t[end]

function count_distinct_neurons(λ₁, λ₂)
    if allzeros(neuron_lags(λ₁)) && allzeros(neuron_lags(λ₂))
        return 1
    elseif allzeros(neuron_lags(λ₁)) || allzeros(neuron_lags(λ₂)) || neuron_lags(λ₁) == neuron_lags(λ₂)
        return 2
    else
        return 3
    end
end

function count_distinct_times(t1, t2)
    if t1 == t2 == 0
        return 1
    elseif t1 == 0 || t2 == 0 || t1 == t2
        return 2
    else
        @assert t1 != t2 != 0
        return 3
    end
end

function _1_channel_seq_class(λ₁, λ₂)
    # All neurons are the same
    # Assume t1 <= t2
    t1 = time_lag(λ₁); t2 = time_lag(λ₂)
    n_distinct_times = count_distinct_times(t1, t2)

    if n_distinct_times == 1
        # All neurons and times are the same
        return 1
    elseif n_distinct_times == 2
        # All neurons are the same, two times
        return 2
    elseif n_distinct_times == 3
        # All neurons are the same, all times distinct
        return 3
    else
        error("Invalid number of distinct times: $n_distinct_times")
    end
end

function _2_channel_seq_class(λ₁, λ₂)
    # Two neurons are the same (two distinct neurons involved)
    # Assume t1 <= t2
    t1 = time_lag(λ₁); t2 = time_lag(λ₂)
    n_distinct_times = count_distinct_times(t1, t2)
    n1 = neuron_lags(λ₁); n2 = neuron_lags(λ₂) 

    if n_distinct_times == 1
        # Two neurons are the same, all times are the same
        return 4
    elseif n_distinct_times == 2
        if (0 == t1 && allzeros(n1)) || (t1 == t2 && n1 == n2) || (0 == t2 && allzeros(n2))
            # Base and first nodes are same; or first and second
            return 6
        elseif (t1 == 0) || (t2 < 0)
            # Synchrony is first
            return 7
        elseif (t2 == 0) || (t1 > 0)
            # Synchrony is second
            return 8
        else
            error("Shouldn't be here")
        end
    elseif n_distinct_times == 3
        if allzeros(n1)
            if (0 < t1)
                return 11
            elseif (t1 < 0)
                if t2 < 0
                    return 10
                elseif t2 > 0
                    return 11
                else
                    error("Shouldn't be here")
                end
            else
                error("Shouldn't be here")
            end
        elseif allzeros(n2)
            if (t2 < 0)
                return 9
            elseif (t2 > 0)
                if t1 > 0 # in between
                    return 10
                elseif t1 < 0
                    return 9
                else
                    error("Shouldn't be here")
                end
            else
                error("Shouldn't be here")
            end
        elseif (n1 == n2)
            if (0 < t1)
                return 9
            elseif (t1 < 0 < t2)
                return 10
            elseif (t2 < 0)
                return 11
            else
                error("Shouldn't be here")
            end
        else
            error("Shouldn't get here.")
        end
    else
        error("Invalid number of same times: $n_distinct_times")
    end
end

function _3_channel_seq_class(λ₁, λ₂)
    # All neurons are distinct
    # Assume t1 <= t2
    t1 = time_lag(λ₁); t2 = time_lag(λ₂)
    n_distinct_times = count_distinct_times(t1, t2)

    if n_distinct_times == 1
        # All neurons are distinct, times are the same
        return 5
    elseif n_distinct_times == 2
        if t1 == 0
            return 13
        elseif t2 == 0
            return 12
        elseif t1 == t2
            if t1 > 0
                return 12
            else
                return 13
            end
        else
            error("Shouldn't get here.")
        end 
    elseif n_distinct_times == 3
        return 14 
    else
        error("Invalid number of same times: $n_distinct_times")
    end
end

function lag_motif_sequence_class(p0::T, p1::T, p2::T) where {T <: NTuple{2}}
    λ₁ = p1 .- p0
    λ₂ = p2 .- p0
    lag_motif_sequence_class(λ₁, λ₂)
end

function (n1, t1, n2, t2)
    n1, t1, n2, t2 = if t1 <= t2
        n1, t1, n2, t2
    else
        n2, t2, n1, t2
    end
end

function lag_motif_sequence_class(λ₁, λ₂) where {N, T, NT <: NTuple{N,T}}

    # order in time
    λ₁, λ₂ = if time_lag(λ₁) <= time_lag(λ₂)
        λ₁, λ₂
    else
        λ₂, λ₁
    end
    n_distinct_neurons = count_distinct_neurons(λ₁, λ₂)

    if n_distinct_neurons == 1
        # All neurons are the same
        return _1_channel_seq_class(λ₁, λ₂)
    elseif n_distinct_neurons == 2
        return _2_channel_seq_class(λ₁, λ₂)
    elseif n_distinct_neurons == 3
        return _3_channel_seq_class(λ₁, λ₂)
    else
        error("Invalid number of same neurons: $n_distinct_neurons")
    end
end