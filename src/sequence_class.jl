# Truncating calculation
using Base.Threads, ThreadsX

function circshift2d(arr, d1, d2)
    @views [arr[1+d1:end, 1+d2:end] arr[1+d1:end, 1:d2];
            arr[1:d1, 1+d2:end] arr[1:d1, 1:d2]]
end


function circshift2d_no_views(arr, d1, d2)
    [arr[1+d1:end, 1+d2:end] arr[1+d1:end, 1:d2];
            arr[1:d1, 1+d2:end] arr[1:d1, 1:d2]]
end

function lag_contribution(data::D, boundary::Periodic, n1::Int, t1::Int, n2::Int, t2::Int, data_λ₁=D(undef, size(data)), data_λ₂=D(undef, size(data))) where {T, D <: AbstractArray{T,2}}
    # Assume ns, ts < size(data)
    contribution = 0
    circshift!(data_λ₁, data, (-n1, -t1))
    circshift!(data_λ₂, data, (-n2, -t2))

    @tturbo for n ∈ axes(data, 1), t ∈ axes(data, 2)
        contribution += data[n,t] * data_λ₁[n,t] * data_λ₂[n,t]
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

function lag_contribution(data::D, boundary::PeriodicExtended, n1::Int, t1::Int, n2::Int, t2::Int, data_λ₁=D(undef, size(data)), data_λ₂=D(undef, size(data))) where {T, D <: Matrix{T}}
    # Periodic in n; extended in t
    # FIXME should validate extension holds lags
    t_start, t_end = boundary.t_bounds
    contribution = 0
    
    circshift!(data_λ₁, data, (-n1, -t1))
    circshift!(data_λ₂, data, (-n2, -t2))

    @tturbo for n ∈ axes(data,1), t ∈ t_start:t_end
        contribution += data[n,t] * data_λ₁[n,t] * data_λ₂[n,t]
    end
    return contribution
end

function sequence_class_tricorr_unrolled(data::AbstractArray, boundary::ZeroPadded, n_max_lag, t_max_lag)
    contributions = zeros(14)
    data = parent(data)

    n_lag_axis = -n_max_lag:n_max_lag
    t_lag_axis = -t_max_lag:t_max_lag

    N_space, N_times = size(data)
    @assert 2n_max_lag < N_space && 2t_max_lag < N_times

    negative_n = n_lag_axis[begin:end÷2]
    positive_n = n_lag_axis[end÷2+2:end]
    negative_t = t_lag_axis[begin:end÷2]
    positive_t = t_lag_axis[end÷2+2:end]

    nonzero_n = [negative_n; positive_n]
    nonzero_t = [negative_t; positive_t]

    # Class I
    contributions[1] = lag_contribution(data, boundary,  0,0,0,0)

    # Class II
    # n1, n2 == 0, 0
    # t1 == 0 or t2 == 0
    @inbounds for t ∈ nonzero_t
        contributions[2] += lag_contribution(data, boundary,  0,t,0,0) + lag_contribution(data, boundary,  0,0,0,t)
    end
    # t1 == t2
    @inbounds for t ∈ nonzero_t
        contributions[2] += lag_contribution(data, boundary,  0,t,0,t)
    end

    # Class III
    # n1, n2 == 0, 0
    # t1 ≠ t2 ≠ 0
    @inbounds for t1 ∈ nonzero_t
        nonzero_nont1_t = filter_element(nonzero_t, t1)
        for t2 ∈ nonzero_nont1_t
            contributions[3] += lag_contribution(data, boundary,  0,t1,0,t2)
        end
    end

    # Class IV
    # t1, t2 == 0, 0
    # n1 == 0 or n2 == 0
    @inbounds for n ∈ nonzero_n
        contributions[4] += lag_contribution(data, boundary,  n,0,0,0) + lag_contribution(data, boundary,  0,0,n,0)
    end
    # t1 == t2
    @inbounds for n ∈ nonzero_n
        contributions[4] += lag_contribution(data, boundary,  n,0,n,0)
    end

    # Class V
    # t1, t2 == 0, 0
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_nonn1_n = filter_element(nonzero_n, n1)
        for n2 ∈ nonzero_nonn1_n
            contributions[5] += lag_contribution(data, boundary,  n1,0,n2,0)
        end
    end

    # Class VI
    # n1, t1 == 0, 0 && n2, t2 ≠ 0, 0
    @inbounds for n2 ∈ nonzero_n, t2 ∈ nonzero_t
        contributions[6] += lag_contribution(data, boundary,  0, 0, n2, t2)
    end
    # n2, t2 == 0, 0 && n1, t1 ≠ 0, 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ nonzero_t
        contributions[6] += lag_contribution(data, boundary,  n1, t1, 0, 0)
    end
    # n1, t1 == n2, t2 ≠ 0, 0
    @inbounds for n ∈ nonzero_n, t ∈ nonzero_t
        contributions[6] += lag_contribution(data, boundary,  n, t, n, t)
    end

    # Class VII
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n2 ≠ 0 && n1 == 0
    @inbounds for t ∈ negative_t, n2 ∈ nonzero_n
        contributions[7] += lag_contribution(data, boundary,  0, t, n2, t)
    end
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n1 ≠ 0 && n2 == 0
    @inbounds for t ∈ negative_t, n1 ∈ nonzero_n
        contributions[7] += lag_contribution(data, boundary,  n1, t, 0, t)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 > 0 && t1 == 0
    @inbounds for n ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary,  n, 0, n, t2)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 > 0 && t2 == 0
    @inbounds for n ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary,  n, t1, n, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 > 0 && t2 == 0
    @inbounds for n2 ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary,  0, t1, n2, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 > 0 && t1 == 0
    @inbounds for n1 ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary,  n1, 0, 0, t2)
    end

    # Class VIII
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n2 ≠ 0 && n1 == 0
    @inbounds for t ∈ positive_t, n2 ∈ nonzero_n
        contributions[8] += lag_contribution(data, boundary,  0, t, n2, t)
    end
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n1 ≠ 0 && n2 == 0
    @inbounds for t ∈ positive_t, n1 ∈ nonzero_n
        contributions[8] += lag_contribution(data, boundary,  n1, t, 0, t)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 < 0 && t1 == 0
    @inbounds for n ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary,  n, 0, n, t2)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 < 0 && t2 == 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary,  n, t1, n, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 < 0 && t2 == 0
    @inbounds for n2 ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary,  0, t1, n2, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 < 0 && t1 == 0
    @inbounds for n1 ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary,  n1, 0, 0, t2)
    end

    # Class IX
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 ≠ t2 > 0
    @inbounds for t1 ∈ positive_t
        positive_nont1_t = filter_element(positive_t, t1)
        # n in inner loop bc filter_element allocates
        @inbounds for t2 ∈ positive_nont1_t, n ∈ nonzero_n
            contributions[9] += lag_contribution(data, boundary,  n, t1, n, t2)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < 0 && t1 > 0
    @inbounds for t1 ∈ positive_t, n2 ∈ nonzero_n, t2 ∈ negative_t
        contributions[9] += lag_contribution(data, boundary,  0, t1, n2, t2)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < 0 && t2 > 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[9] += lag_contribution(data, boundary,  n1, t1, 0, t2)
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < t1 < 0
    @inbounds for t1 ∈ negative_t[begin+1:end], n2 ∈ nonzero_n
        @inbounds for t2 ∈ negative_t[begin]:(t1-1)
            contributions[9] += lag_contribution(data, boundary,  0, t1, n2, t2)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < t2 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        for t2 ∈ (t1+1):-1
            contributions[9] += lag_contribution(data, boundary,  n1, t1, 0, t2)
        end
    end

    # Class X
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 < t2 < 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        for t2 ∈ (t1+1):-1
            contributions[10] += lag_contribution(data, boundary,  n, t1, n, t2)
        end
    end
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t2 < t1 < 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[10] += lag_contribution(data, boundary,  n, t1, n, t2)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 > 0 && t1 < 0
    @inbounds for t1 ∈ negative_t, n2 ∈ nonzero_n, t2 ∈ positive_t
        contributions[10] += lag_contribution(data, boundary,  0, t1, n2, t2)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 > 0 && t2 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[10] += lag_contribution(data, boundary,  n1, t1, 0, t2)
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && 0 < t1 < t2
    @inbounds for t1 ∈ positive_t[begin:end-1], n2 ∈ nonzero_n
        for t2 ∈ (t1+1):positive_t[end]
            contributions[10] += lag_contribution(data, boundary,  0, t1, n2, t2)
        end
    end
    # Assume left point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && 0 < t2 < t1
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t[begin+1:end]
        for t2 ∈ positive_t[begin]:(t1-1)
            contributions[10] += lag_contribution(data, boundary,  n1, t1, 0, t2)
        end
    end

    # Class XI
    # Assume left point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; 0 < t1 < t2
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t[begin:end-1]
        for t2 ∈ (t1+1):positive_t[end]
            contributions[11] += lag_contribution(data, boundary,  n1, t1, 0, t2)
        end
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; 0 < t2 < t1
    @inbounds for t1 ∈ positive_t[begin+1:end], n2 ∈ nonzero_n 
        for t2 ∈ 1:(t1-1)
            contributions[11] += lag_contribution(data, boundary,  0, t1, n2, t2)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; t2 < t1 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[11] += lag_contribution(data, boundary,  n1, t1, 0, t2)
        end
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; t1 < t2 < 0
    @inbounds for t1 ∈ negative_t[begin:end-1], n2 ∈ nonzero_n 
        for t2 ∈ (t1+1):-1
            contributions[11] += lag_contribution(data, boundary,  0, t1, n2, t2)
        end
    end
    # Assume middle point (0,0); n1 left
    # n1 == n2 ≠ 0; t1 < 0 < t2
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[11] += lag_contribution(data, boundary,  n, t1, n, t2)
    end
    # Assume middle point (0,0); n2 left
    # n1 == n2 ≠ 0; t2 < 0 < t1
    @inbounds for n ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[11] += lag_contribution(data, boundary,  n, t1, n, t2)
    end

    # Class XII
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 > 0
        @inbounds for t ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, boundary,  n1, t, n2, t)
        end
        # Assume n1 odd: t1 < 0; t2 == 0
        @inbounds for t1 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, boundary,  n1, t1, n2, 0)
        end
        # Assume n2 odd: t2 < 0; t1 == 0
        @inbounds for t2 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, boundary,  n1, 0, n2, t2)
        end
    end

    # Class XIII
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 < 0
        @inbounds for t ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, boundary,  n1, t, n2, t)
        end
        # Assume n1 odd: t1 > 0; t2 == 0
        @inbounds for t1 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, boundary,  n1, t1, n2, 0)
        end
        # Assume n2 odd: t2 > 0; t1 == 0
        @inbounds for t2 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, boundary,  n1, 0, n2, t2)
        end
    end

    # Class XIV
    # n1 ≠ n2 ≠ 0
    cont = 0
    n1_with_nonzero_filtered = [(n1, filter_element(nonzero_n, n1)) for n1 ∈ nonzero_n]
    t1_with_nonzero_filtered = [(t1, filter_element(nonzero_t, t1)) for t1 ∈ nonzero_t]
    @inbounds for (n1, nonzero_notn1) ∈ n1_with_nonzero_filtered, (t1, nonzero_nott1) ∈ t1_with_nonzero_filtered
        for n2 ∈ nonzero_notn1, t2 ∈ nonzero_nott1
            contributions[14] += lag_contribution(data, boundary,  n1, t1, n2, t2)
        end
    end

    return contributions ./ calculate_scaling_factor(data, boundary)

end

function sequence_class_tricorr_unrolled(data::AbstractArray,boundary::Union{Periodic,PeriodicExtended}, n_max_lag, t_max_lag)
    contributions = zeros(14)
    data = parent(data)
    lag1_cache = typeof(data)(undef, size(data))
    lag2_cache = typeof(data)(undef, size(data))

    n_lag_axis = -n_max_lag:n_max_lag
    t_lag_axis = -t_max_lag:t_max_lag

    N_space, N_times = size(data)
    @assert (2n_max_lag < N_space && 2t_max_lag < N_times) "Need 2$(n_max_lag) < $N_space and 2$t_max_lag < $N_times"

    negative_n = n_lag_axis[begin:end÷2]
    positive_n = n_lag_axis[end÷2+2:end]
    negative_t = t_lag_axis[begin:end÷2]
    positive_t = t_lag_axis[end÷2+2:end]

    nonzero_n = [negative_n; positive_n]
    nonzero_t = [negative_t; positive_t]

    # Class I
    contributions[1] = lag_contribution(data, boundary, 0,0,0,0,lag1_cache,lag2_cache)

    # Class II
    # n1, n2 == 0, 0
    # t1 == 0 or t2 == 0
    @inbounds for t ∈ nonzero_t
        contributions[2] += lag_contribution(data, boundary, 0,t,0,0,lag1_cache,lag2_cache) + lag_contribution(data, boundary, 0,0,0,t,lag1_cache,lag2_cache)
    end
    # t1 == t2
    @inbounds for t ∈ nonzero_t
        contributions[2] += lag_contribution(data, boundary, 0,t,0,t,lag1_cache,lag2_cache)
    end

    # Class III
    # n1, n2 == 0, 0
    # t1 ≠ t2 ≠ 0
    @inbounds for t1 ∈ nonzero_t
        nonzero_nont1_t = filter_element(nonzero_t, t1)
        for t2 ∈ nonzero_nont1_t
            contributions[3] += lag_contribution(data, boundary, 0,t1,0,t2,lag1_cache,lag2_cache)
        end
    end

    # Class IV
    # t1, t2 == 0, 0
    # n1 == 0 or n2 == 0
    @inbounds for n ∈ nonzero_n
        contributions[4] += lag_contribution(data, boundary, n,0,0,0,lag1_cache,lag2_cache) + lag_contribution(data, boundary, 0,0,n,0,lag1_cache,lag2_cache)
    end
    # t1 == t2
    @inbounds for n ∈ nonzero_n
        contributions[4] += lag_contribution(data, boundary, n,0,n,0,lag1_cache,lag2_cache)
    end

    # Class V
    # t1, t2 == 0, 0
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_nonn1_n = filter_element(nonzero_n, n1)
        for n2 ∈ nonzero_nonn1_n
            contributions[5] += lag_contribution(data, boundary, n1,0,n2,0,lag1_cache,lag2_cache)
        end
    end

    # Class VI
    # n1, t1 == 0, 0 && n2, t2 ≠ 0, 0
    @inbounds for n2 ∈ nonzero_n, t2 ∈ nonzero_t
        contributions[6] += lag_contribution(data, boundary, 0, 0, n2, t2,lag1_cache,lag2_cache)
    end
    # n2, t2 == 0, 0 && n1, t1 ≠ 0, 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ nonzero_t
        contributions[6] += lag_contribution(data, boundary, n1, t1, 0, 0,lag1_cache,lag2_cache)
    end
    # n1, t1 == n2, t2 ≠ 0, 0
    @inbounds for n ∈ nonzero_n, t ∈ nonzero_t
        contributions[6] += lag_contribution(data, boundary, n, t, n, t,lag1_cache,lag2_cache)
    end

    # Class VII
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n2 ≠ 0 && n1 == 0
    @inbounds for t ∈ negative_t, n2 ∈ nonzero_n
        contributions[7] += lag_contribution(data, boundary, 0, t, n2, t,lag1_cache,lag2_cache)
    end
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n1 ≠ 0 && n2 == 0
    @inbounds for t ∈ negative_t, n1 ∈ nonzero_n
        contributions[7] += lag_contribution(data, boundary, n1, t, 0, t,lag1_cache,lag2_cache)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 > 0 && t1 == 0
    @inbounds for n ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary, n, 0, n, t2,lag1_cache,lag2_cache)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 > 0 && t2 == 0
    @inbounds for n ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary, n, t1, n, 0,lag1_cache,lag2_cache)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 > 0 && t2 == 0
    @inbounds for n2 ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary, 0, t1, n2, 0,lag1_cache,lag2_cache)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 > 0 && t1 == 0
    @inbounds for n1 ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += lag_contribution(data, boundary, n1, 0, 0, t2,lag1_cache,lag2_cache)
    end

    # Class VIII
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n2 ≠ 0 && n1 == 0
    @inbounds for t ∈ positive_t, n2 ∈ nonzero_n
        contributions[8] += lag_contribution(data, boundary, 0, t, n2, t,lag1_cache,lag2_cache)
    end
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n1 ≠ 0 && n2 == 0
    @inbounds for t ∈ positive_t, n1 ∈ nonzero_n
        contributions[8] += lag_contribution(data, boundary, n1, t, 0, t,lag1_cache,lag2_cache)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 < 0 && t1 == 0
    @inbounds for n ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary, n, 0, n, t2,lag1_cache,lag2_cache)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 < 0 && t2 == 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary, n, t1, n, 0,lag1_cache,lag2_cache)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 < 0 && t2 == 0
    @inbounds for n2 ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary, 0, t1, n2, 0,lag1_cache,lag2_cache)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 < 0 && t1 == 0
    @inbounds for n1 ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += lag_contribution(data, boundary, n1, 0, 0, t2,lag1_cache,lag2_cache)
    end

    # Class IX
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 ≠ t2 > 0
    @inbounds for t1 ∈ positive_t
        positive_nont1_t = filter_element(positive_t, t1)
        # n in inner loop bc filter_element allocates
        @inbounds for t2 ∈ positive_nont1_t, n ∈ nonzero_n
            contributions[9] += lag_contribution(data, boundary, n, t1, n, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < 0 && t1 > 0
    @inbounds for t1 ∈ positive_t, n2 ∈ nonzero_n, t2 ∈ negative_t
        contributions[9] += lag_contribution(data, boundary, 0, t1, n2, t2,lag1_cache,lag2_cache)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < 0 && t2 > 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[9] += lag_contribution(data, boundary, n1, t1, 0, t2,lag1_cache,lag2_cache)
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < t1 < 0
    @inbounds for t1 ∈ negative_t[begin+1:end], n2 ∈ nonzero_n
        @inbounds for t2 ∈ negative_t[begin]:(t1-1)
            contributions[9] += lag_contribution(data, boundary, 0, t1, n2, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < t2 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        for t2 ∈ (t1+1):-1
            contributions[9] += lag_contribution(data, boundary, n1, t1, 0, t2,lag1_cache,lag2_cache)
        end
    end

    # Class X
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 < t2 < 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        for t2 ∈ (t1+1):-1
            contributions[10] += lag_contribution(data, boundary, n, t1, n, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t2 < t1 < 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[10] += lag_contribution(data, boundary, n, t1, n, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 > 0 && t1 < 0
    @inbounds for t1 ∈ negative_t, n2 ∈ nonzero_n, t2 ∈ positive_t
        contributions[10] += lag_contribution(data, boundary, 0, t1, n2, t2,lag1_cache,lag2_cache)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 > 0 && t2 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[10] += lag_contribution(data, boundary, n1, t1, 0, t2,lag1_cache,lag2_cache)
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && 0 < t1 < t2
    @inbounds for t1 ∈ positive_t[begin:end-1], n2 ∈ nonzero_n
        for t2 ∈ (t1+1):positive_t[end]
            contributions[10] += lag_contribution(data, boundary, 0, t1, n2, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume left point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && 0 < t2 < t1
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t[begin+1:end]
        for t2 ∈ positive_t[begin]:(t1-1)
            contributions[10] += lag_contribution(data, boundary, n1, t1, 0, t2,lag1_cache,lag2_cache)
        end
    end

    # Class XI
    # Assume left point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; 0 < t1 < t2
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t[begin:end-1]
        for t2 ∈ (t1+1):positive_t[end]
            contributions[11] += lag_contribution(data, boundary, n1, t1, 0, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; 0 < t2 < t1
    @inbounds for t1 ∈ positive_t[begin+1:end], n2 ∈ nonzero_n 
        for t2 ∈ 1:(t1-1)
            contributions[11] += lag_contribution(data, boundary, 0, t1, n2, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; t2 < t1 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[11] += lag_contribution(data, boundary, n1, t1, 0, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; t1 < t2 < 0
    @inbounds for t1 ∈ negative_t[begin:end-1], n2 ∈ nonzero_n 
        for t2 ∈ (t1+1):-1
            contributions[11] += lag_contribution(data, boundary, 0, t1, n2, t2,lag1_cache,lag2_cache)
        end
    end
    # Assume middle point (0,0); n1 left
    # n1 == n2 ≠ 0; t1 < 0 < t2
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[11] += lag_contribution(data, boundary, n, t1, n, t2,lag1_cache,lag2_cache)
    end
    # Assume middle point (0,0); n2 left
    # n1 == n2 ≠ 0; t2 < 0 < t1
    @inbounds for n ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[11] += lag_contribution(data, boundary, n, t1, n, t2,lag1_cache,lag2_cache)
    end

    # Class XII
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 > 0
        @inbounds for t ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, boundary, n1, t, n2, t,lag1_cache,lag2_cache)
        end
        # Assume n1 odd: t1 < 0; t2 == 0
        @inbounds for t1 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, boundary, n1, t1, n2, 0,lag1_cache,lag2_cache)
        end
        # Assume n2 odd: t2 < 0; t1 == 0
        @inbounds for t2 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, boundary, n1, 0, n2, t2,lag1_cache,lag2_cache)
        end
    end

    # Class XIII
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 < 0
        @inbounds for t ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, boundary, n1, t, n2, t,lag1_cache,lag2_cache)
        end
        # Assume n1 odd: t1 > 0; t2 == 0
        @inbounds for t1 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, boundary, n1, t1, n2, 0,lag1_cache,lag2_cache)
        end
        # Assume n2 odd: t2 > 0; t1 == 0
        @inbounds for t2 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, boundary, n1, 0, n2, t2,lag1_cache,lag2_cache)
        end
    end

    # Class XIV
    # n1 ≠ n2 ≠ 0
    cont = 0
    n1_with_nonzero_filtered = [(n1, filter_element(nonzero_n, n1)) for n1 ∈ nonzero_n]
    t1_with_nonzero_filtered = [(t1, filter_element(nonzero_t, t1)) for t1 ∈ nonzero_t]
    @inbounds for (n1, nonzero_notn1) ∈ n1_with_nonzero_filtered, (t1, nonzero_nott1) ∈ t1_with_nonzero_filtered
        for n2 ∈ nonzero_notn1, t2 ∈ nonzero_nott1
            contributions[14] += lag_contribution(data, boundary, n1, t1, n2, t2,lag1_cache,lag2_cache)
        end
    end

    return contributions ./ calculate_scaling_factor(data, boundary)

end


# FIXME 
# instead of filtering the element, track the el_idx,
# then split inner loop to go begin:el_idx-1 and el_idx+1:end
function filter_element(arr, el)
    filter(≠(el), arr)
end

function sequence_class_tricorr_unrolled(tricorr::TripleCorrelation)
    arr = tricorr.arr
    contributions = zeros(14)

    n_lag_axis = axes(arr, 1)
    t_lag_axis = axes(arr, 2)
    @assert length(n_lag_axis) % 2 == length(t_lag_axis) % 2 == 1

    negative_n = n_lag_axis[begin:end÷2]
    positive_n = n_lag_axis[end÷2+2:end]
    negative_t = t_lag_axis[begin:end÷2]
    positive_t = t_lag_axis[end÷2+2:end]

    nonzero_n = [negative_n; positive_n]
    nonzero_t = [negative_t; positive_t]

    # Class I
    contributions[1] = arr[0,0,0,0]

    # Class II
    # n1, n2 == 0, 0
    # t1 == 0 or t2 == 0
    contributions[2] += sum(arr[0,nonzero_t,0,0]) + sum(arr[0,0,0,nonzero_t])
    # t1 == t2
    for t ∈ nonzero_t
        contributions[2] += arr[0,t,0,t]
    end

    # Class III
    # n1, n2 == 0, 0
    # t1 ≠ t2 ≠ 0
    for t1 ∈ nonzero_t
        nonzero_nont1_t = filter_element(nonzero_t, t1)
        for t2 ∈ nonzero_nont1_t
            contributions[3] += arr[0,t1,0,t2]
        end
    end

    # Class IV
    # t1, t2 == 0, 0
    # n1 == 0 or n2 == 0
    contributions[4] += sum(arr[nonzero_n,0,0,0]) + sum(arr[0,0,nonzero_n,0])
    # t1 == t2
    for n ∈ nonzero_n
        contributions[4] += arr[n,0,n,0]
    end

    # Class V
    # t1, t2 == 0, 0
    # n1 ≠ n2 ≠ 0
    for n1 ∈ nonzero_n
        nonzero_nonn1_n = filter_element(nonzero_n, n1)
        for n2 ∈ nonzero_nonn1_n
            contributions[5] += arr[n1,0,n2,0]
        end
    end

    # Class VI
    # n1, t1 == 0, 0 && n2, t2 ≠ 0, 0
    for n2 ∈ nonzero_n, t2 ∈ nonzero_t
        contributions[6] += arr[0, 0, n2, t2]
    end
    # n2, t2 == 0, 0 && n1, t1 ≠ 0, 0
    for n1 ∈ nonzero_n, t1 ∈ nonzero_t
        contributions[6] += arr[n1, t1, 0, 0]
    end
    # n1, t1 == n2, t2 ≠ 0, 0
    for n ∈ nonzero_n, t ∈ nonzero_t
        contributions[6] += arr[n, t, n, t]
    end

    # Class VII
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n2 ≠ 0 && n1 == 0
    for t ∈ negative_t, n2 ∈ nonzero_n
        contributions[7] += arr[0, t, n2, t]
    end
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n1 ≠ 0 && n2 == 0
    for t ∈ negative_t, n1 ∈ nonzero_n
        contributions[7] += arr[n1, t, 0, t]
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 > 0 && t1 == 0
    for n ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += arr[n, 0, n, t2]
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 > 0 && t2 == 0
    for n ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += arr[n, t1, n, 0]
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 > 0 && t2 == 0
    for n2 ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += arr[0, t1, n2, 0]
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 > 0 && t1 == 0
    for n1 ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += arr[n1, 0, 0, t2]
    end

    # Class VIII
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n2 ≠ 0 && n1 == 0
    for t ∈ positive_t, n2 ∈ nonzero_n
        contributions[8] += arr[0, t, n2, t]
    end
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n1 ≠ 0 && n2 == 0
    for t ∈ positive_t, n1 ∈ nonzero_n
        contributions[8] += arr[n1, t, 0, t]
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 < 0 && t1 == 0
    for n ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += arr[n, 0, n, t2]
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 < 0 && t2 == 0
    for n ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += arr[n, t1, n, 0]
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 < 0 && t2 == 0
    for n2 ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += arr[0, t1, n2, 0]
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 < 0 && t1 == 0
    for n1 ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += arr[n1, 0, 0, t2]
    end

    # Class IX
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 ≠ t2 > 0
    for t1 ∈ positive_t
        positive_nont1_t = filter_element(positive_t, t1)
        # n in inner loop bc filter_element allocates
        for t2 ∈ positive_nont1_t, n ∈ nonzero_n
            contributions[9] += arr[n, t1, n, t2]
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < 0 && t1 > 0
    for t1 ∈ positive_t, n2 ∈ nonzero_n, t2 ∈ negative_t
        contributions[9] += arr[0, t1, n2, t2]
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < 0 && t2 > 0
    for n1 ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[9] += arr[n1, t1, 0, t2]
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < t1 < 0
    for t1 ∈ negative_t[begin+1:end], n2 ∈ nonzero_n
        contributions[9] += sum(arr[0, t1, n2, begin:(t1-1)])
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < t2 < 0
    for n1 ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        contributions[9] += sum(arr[n1, t1, 0, (t1+1):-1])
    end

    # Class X
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 < t2 < 0
    for n ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        for t2 ∈ (t1+1):-1
            contributions[10] += arr[n, t1, n, t2]
        end
    end
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t2 < t1 < 0
    for n ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[10] += arr[n, t1, n, t2]
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 > 0 && t1 < 0
    for t1 ∈ negative_t, n2 ∈ nonzero_n, t2 ∈ positive_t
        contributions[10] += arr[0, t1, n2, t2]
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 > 0 && t2 < 0
    for n1 ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[10] += arr[n1, t1, 0, t2]
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && 0 < t1 < t2
    for t1 ∈ positive_t[begin:end-1], n2 ∈ nonzero_n
        for t2 ∈ (t1+1):positive_t[end]
            contributions[10] += arr[0, t1, n2, t2]
        end
    end
    # Assume left point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && 0 < t2 < t1
    for n1 ∈ nonzero_n, t1 ∈ positive_t[begin+1:end]
        for t2 ∈ positive_t[begin]:(t1-1)
            contributions[10] += arr[n1, t1, 0, t2]
        end
    end

    # Class XI
    # Assume left point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; 0 < t1 < t2
    for n1 ∈ nonzero_n, t1 ∈ positive_t[begin:end-1]
        for t2 ∈ (t1+1):positive_t[end]
            contributions[11] += arr[n1, t1, 0, t2]
        end
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; 0 < t2 < t1
    for t1 ∈ positive_t[begin+1:end], n2 ∈ nonzero_n 
        for t2 ∈ 1:(t1-1)
            contributions[11] += arr[0, t1, n2, t2]
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; t2 < t1 < 0
    for n1 ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[11] += arr[n1, t1, 0, t2]
        end
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; t1 < t2 < 0
    for t1 ∈ negative_t[begin:end-1], n2 ∈ nonzero_n 
        for t2 ∈ (t1+1):-1
            contributions[11] += arr[0, t1, n2, t2]
        end
    end
    # Assume middle point (0,0); n1 left
    # n1 == n2 ≠ 0; t1 < 0 < t2
    for n ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[11] += arr[n, t1, n, t2]
    end
    # Assume middle point (0,0); n2 left
    # n1 == n2 ≠ 0; t2 < 0 < t1
    for n ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[11] += arr[n, t1, n, t2]
    end

    # Class XII
    # n1 ≠ n2 ≠ 0
    for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 > 0
        for t ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[12] += arr[n1, t, n2, t]
        end
        # Assume n1 odd: t1 < 0; t2 == 0
        for t1 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += arr[n1, t1, n2, 0]
        end
        # Assume n2 odd: t2 < 0; t1 == 0
        for t2 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += arr[n1, 0, n2, t2]
        end
    end

    # Class XIII
    # n1 ≠ n2 ≠ 0
    for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 < 0
        for t ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[13] += arr[n1, t, n2, t]
        end
        # Assume n1 odd: t1 > 0; t2 == 0
        for t1 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += arr[n1, t1, n2, 0]
        end
        # Assume n2 odd: t2 > 0; t1 == 0
        for t2 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += arr[n1, 0, n2, t2]
        end
    end

    # Class XIV
    # n1 ≠ n2 ≠ 0
    n1_with_nonzero_filtered = [(n1, filter_element(nonzero_n, n1)) for n1 ∈ nonzero_n]
    t1_with_nonzero_filtered = [(t1, filter_element(nonzero_t, t1)) for t1 ∈ nonzero_t]
    for (n1, nonzero_notn1) ∈ n1_with_nonzero_filtered, (t1, nonzero_nott1) ∈ t1_with_nonzero_filtered
        for n2 ∈ nonzero_notn1, t2 ∈ nonzero_nott1
            contributions[14] += arr[n1, t1, n2, t2]
        end
    end

    return contributions

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

function sequence_class_tricorr!(class_contribution::AbstractVector, src::AbstractArray{T}, boundary::AbstractBoundaryCondition, space_max_lag, time_max_lag, lags_classifier::Function) where T
    src = parent(src)
    space_lag_range = -(space_max_lag):(space_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    class_contribution .= 0
    for n1 ∈ space_lag_range, n2 ∈ space_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        contribution = lag_contribution(src, boundary, n1, t1, n2, t2)
        class_contribution[class] += contribution
    end
    class_contribution ./= calculate_scaling_factor(src, boundary)
end

function sequence_class_tricorr(src::OffsetArray, args...)
    sequence_class_tricorr(parent(src), args...)
end

# Periodic calculation

function sequence_class_tricorr!(class_contribution::AbstractVector, src::SRC, boundary::Periodic, space_max_lag, time_max_lag, lags_classifier::Function) where {T, SRC<:AbstractArray{T}}
    src = parent(src)
    lag1_cache = typeof(src)(undef, size(src))
    lag2_cache = typeof(src)(undef, size(src))

    space_lag_range = -(space_max_lag):(space_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    class_contribution .= 0
    for n1 ∈ space_lag_range, n2 ∈ space_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        contribution = lag_contribution(src, boundary, n1, t1, n2, t2, lag1_cache, lag2_cache)
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
        if (0 == t1 && 0 == n1) || (t1 == t2 && n1 == n2) || (0 == t2 && 0 == n2)
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
        if (n1 == 0)
            if (0 < t1)
                return 10
            elseif (t1 < 0)
                if t2 < 0
                    return 11
                elseif t2 > 0
                    return 10
                else
                    error("Shouldn't be here")
                end
            else
                error("Shouldn't be here")
            end
        elseif (n2 == 0)
            if (t2 < 0)
                return 9
            elseif (t2 > 0)
                if t1 > 0 # in between
                    return 11
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
                return 11
            elseif (t2 < 0)
                return 10
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

function lag_motif_sequence_class(tup::NTuple{4})
    lag_motif_sequence_class(tup[1], tup[2], tup[3], tup[4])
end

function lag_motif_sequence_class(n1, n2, t1, t2)
    n_distinct_neurons = count_distinct(0, n1, n2)

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