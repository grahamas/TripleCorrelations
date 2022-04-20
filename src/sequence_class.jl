# Truncating calculation
using Base.Threads, ThreadsX

function lag_contribution(data::D, boundary::Periodic, p1::NTuple{N,Int}, p2::NTuple{N,Int}, data_λ₁=D(undef, size(data)), data_λ₂=D(undef, size(data))) where {N, T, D <: AbstractArray{T}}
    # Assume ns, ts < size(data)
    contribution = 0
    circshift!(data_λ₁, data, .-p1)
    circshift!(data_λ₂, data, .-p2)

    @tturbo for p0 ∈ CartesianIndices(data)
        contribution += data[p0] * data_λ₁[p0] * data_λ₂[p0]
    end
    return contribution
end

function lag_contribution(data::D, boundary::ZeroPadded, n1::Int, t1::Int, n2::Int, t2::Int) where {T, D <: AbstractArray{T,2}}
    contribution = 0

    n_start = max(1 - min(n1, n2), 1)
    t_start = max(1 - min(t1, t2), 1)
    n_end = min(size(data,1) - max(n1,n2), size(data,1))
    t_end = min(size(data,2) - max(t1,t2), size(data,2))

    base_n_range = n_start:n_end
    base_t_range = t_start:t_end

    λ₁_n_range = base_n_range .+ n1
    λ₁_t_range = base_t_range .+ t1

    λ₂_n_range = base_n_range .+ n2
    λ₂_t_range = base_t_range .+ t2

    A = @view data[base_n_range, base_t_range]
    A_λ₁ = @view data[λ₁_n_range, λ₁_t_range]
    A_λ₂ = @view data[λ₂_n_range, λ₂_t_range]

    @tturbo for n ∈ axes(A, 1), t ∈ axes(A, 2)
        contribution += A[n,t] * A_λ₁[n,t] * A_λ₂[n,t]
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
# function sequence_class_tricorr!(class_contribution::AbstractVector, src, boundary:: space_max_lag, time_max_lag, lags_classifier::Function)
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

function sequence_class_tricorr(src, boundary::AbstractBoundaryCondition, space_max_lag, time_max_lag)
    N_network_classifications = 14
    network_class_contributions = Array{Float64}(undef, N_network_classifications)
    lags_classifier = lag_motif_sequence_class

    sequence_class_tricorr!(network_class_contributions, src, boundary, space_max_lag, time_max_lag, lags_classifier)
end

function sequence_class_tricorr!(class_contribution::AbstractVector, src::AbstractArray{T}, boundary::AbstractBoundaryCondition, max_lags, lags_classifier::Function) where T
    src = parent(src)

    lag_ranges = map(((start,stop),) -> UnitRange(start,stop), zip(.-max_lags, max_lags))

    class_contribution .= 0
    for λ₁ ∈ IterTools.product(lag_ranges...), λ₂ ∈ IterTools.product(lag_ranges...)
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

function sequence_class_tricorr!(class_contribution::AbstractVector, src::SRC, boundary::Union{Periodic,PeriodicExtended}, max_lags, lags_classifier::Function) where {T, SRC<:AbstractArray{T}}
    src = parent(src)
    lag1_cache = typeof(src)(undef, size(src))
    lag2_cache = typeof(src)(undef, size(src))

    lag_ranges = map(((start,stop),) -> UnitRange(start,stop), zip(.-max_lags, max_lags))

    class_contribution .= 0
    for λ₁ ∈ IterTools.product(lag_ranges...), λ₂ ∈ IterTools.product(lag_ranges...)
        class = lags_classifier(λ₁, λ₂)
        contribution = lag_contribution(src, boundary, λ₁, λ₂, lag1_cache, lag2_cache)
        class_contribution[class] += contribution
    end
    class_contribution ./= calculate_scaling_factor(src, boundary)
end




function count_distinct(args...)
    return length(unique(args))
end

function _1_channel_seq_class(n1, n2, t1, t2)
    # All neurons are the same
    # Assume t1 <= t2
    n_distinct_times = count_distinct(0, t1, t2)

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

function _2_channel_seq_class(n1, n2, t1, t2)
    # Two neurons are the same (two distinct neurons involved)
    # Assume t1 <= t2
    n_distinct_times = count_distinct(0, t1, t2)

    if n_distinct_times == 1
        # Two neurons are the same, all times are the same
        return 4
    elseif n_distinct_times == 2
        if (0 == t1 && all(iszero.(n1))) || (t1 == t2 && n1 == n2) || (0 == t2 && all(iszero.(n2)))
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
        if all(iszero.(n1))
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
        elseif all(iszero.(n2))
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

function _3_channel_seq_class(n1, n2, t1, t2)
    # All neurons are distinct
    # Assume t1 <= t2
    n_distinct_times = count_distinct(0, t1, t2)

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

function lag_motif_sequence_class(λ₁::NT, λ₂::NT) where {N, T, NT <: NTuple{N,T}}
    n1 = λ₁[1:end-1]; t1 = λ₁[end]
    n2 = λ₂[1:end-1]; t2 = λ₂[end]
    n_distinct_neurons = count_distinct(zero.(λ₁), n1, n2)

    n1, n2, t1, t2 = if t1 < t2
        (n1, n2, t1, t2)
    else
        (n2, n1, t2, t1)
    end
    # Assume below that t1 <= t2

    if n_distinct_neurons == 1
        # All neurons are the same
        return _1_channel_seq_class(n1, n2, t1, t2)
    elseif n_distinct_neurons == 2
        return _2_channel_seq_class(n1, n2, t1, t2)
    elseif n_distinct_neurons == 3
        return _3_channel_seq_class(n1, n2, t1, t2)
    else
        error("Invalid number of same neurons: $n_distinct_neurons")
    end
end