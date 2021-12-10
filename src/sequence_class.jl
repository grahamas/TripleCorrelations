# Truncating calculation

function lag_contribution(data, n_offset, t_offset, n1::Int, t1::Int, n2::Int, t2::Int)
    contribution = 0
    @inbounds for data_n ∈ n_offset, data_t ∈ t_offset
        contribution += data[data_n, data_t] * 
            data[data_n+n1, data_t+t1] *
            data[data_n+n2, data_t+t2]
    end
    return contribution
end

function lag_sequence_class_contribution_unrolled(data::AbstractArray, λ_max::NTuple{2})
    contributions = zeros(14)
    data = parent(data)

    n_max_lag, t_max_lag = λ_max
    n_lag_axis = -n_max_lag:n_max_lag
    t_lag_axis = -t_max_lag:t_max_lag

    N_space, N_times = size(data)
    t_offset = (1-minimum(t_lag_axis)):(N_times-maximum(t_lag_axis))
    n_offset = (1-minimum(n_lag_axis)):(N_space-maximum(n_lag_axis))
    @show data t_offset n_offset n_lag_axis 

    negative_n = n_lag_axis[begin:end÷2]
    positive_n = n_lag_axis[end÷2+2:end]
    negative_t = t_lag_axis[begin:end÷2]
    positive_t = t_lag_axis[end÷2+2:end]

    nonzero_n = [negative_n; positive_n]
    nonzero_t = [negative_t; positive_t]

    # Class I
    contributions[1] = lag_contribution(data, n_offset, t_offset, 0,0,0,0)

    # Class II
    # n1, n2 == 0, 0
    # t1 == 0 or t2 == 0
    for t ∈ nonzero_t
        contributions[2] += lag_contribution(data, n_offset, t_offset, 0,t,0,0) + lag_contribution(data, n_offset, t_offset, 0,0,0,t)
    end
    # t1 == t2
    for t ∈ nonzero_t
        contributions[2] += lag_contribution(data, n_offset, t_offset, 0,t,0,t)
    end

    # Class III
    # n1, n2 == 0, 0
    # t1 ≠ t2 ≠ 0
    for t1 ∈ nonzero_t
        nonzero_nont1_t = filter_element(nonzero_t, t1)
        for t2 ∈ nonzero_nont1_t
            contributions[3] += lag_contribution(data, n_offset, t_offset, 0,t1,0,t2)
        end
    end

    # Class IV
    # t1, t2 == 0, 0
    # n1 == 0 or n2 == 0
    for n ∈ nonzero_n
        contributions[4] += lag_contribution(data, n_offset, t_offset, n,0,0,0) + lag_contribution(data, n_offset, t_offset, 0,0,n,0)
    end
    # t1 == t2
    for n ∈ nonzero_n
        contributions[4] += lag_contribution(data, n_offset, t_offset, n,0,n,0)
    end

    # Class V
    # t1, t2 == 0, 0
    # n1 ≠ n2 ≠ 0
    for n1 ∈ nonzero_n
        nonzero_nonn1_n = filter_element(nonzero_n, n1)
        for n2 ∈ nonzero_nonn1_n
            contributions[5] += lag_contribution(data, n_offset, t_offset, n1,0,n2,0)
        end
    end

    # Class VI
    # n1, t1 == 0, 0 && n2, t2 ≠ 0, 0
    for n2 ∈ nonzero_n, t2 ∈ nonzero_t
        contributions[6] += lag_contribution(data, n_offset, t_offset, 0, 0, n2, t2)
    end
    # n2, t2 == 0, 0 && n1, t1 ≠ 0, 0
    for n1 ∈ nonzero_n, t1 ∈ nonzero_t
        contributions[6] += lag_contribution(data, n_offset, t_offset, n1, t1, 0, 0)
    end
    # n1, t1 == n2, t2 ≠ 0, 0
    for n ∈ nonzero_n, t ∈ nonzero_t
        contributions[6] += lag_contribution(data, n_offset, t_offset, n, t, n, t)
    end

    # Class VII
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n2 ≠ 0 && n1 == 0
    for t ∈ negative_t, n2 ∈ nonzero_n
        contributions[7] += lag_contribution(data, n_offset, t_offset, 0, t, n2, t)
    end
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n1 ≠ 0 && n2 == 0
    for t ∈ negative_t, n1 ∈ nonzero_n
        contributions[7] += lag_contribution(data, n_offset, t_offset, n1, t, 0, t)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 > 0 && t1 == 0
    for n ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += lag_contribution(data, n_offset, t_offset, n, 0, n, t2)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 > 0 && t2 == 0
    for n ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += lag_contribution(data, n_offset, t_offset, n, t1, n, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 > 0 && t2 == 0
    for n2 ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 > 0 && t1 == 0
    for n1 ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += lag_contribution(data, n_offset, t_offset, n1, 0, 0, t2)
    end

    # Class VIII
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n2 ≠ 0 && n1 == 0
    for t ∈ positive_t, n2 ∈ nonzero_n
        contributions[8] += lag_contribution(data, n_offset, t_offset, 0, t, n2, t)
    end
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n1 ≠ 0 && n2 == 0
    for t ∈ positive_t, n1 ∈ nonzero_n
        contributions[8] += lag_contribution(data, n_offset, t_offset, n1, t, 0, t)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 < 0 && t1 == 0
    for n ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += lag_contribution(data, n_offset, t_offset, n, 0, n, t2)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 < 0 && t2 == 0
    for n ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += lag_contribution(data, n_offset, t_offset, n, t1, n, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 < 0 && t2 == 0
    for n2 ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 < 0 && t1 == 0
    for n1 ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += lag_contribution(data, n_offset, t_offset, n1, 0, 0, t2)
    end

    # Class IX
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 ≠ t2 > 0
    for t1 ∈ positive_t
        positive_nont1_t = filter_element(positive_t, t1)
        # n in inner loop bc filter_element allocates
        for t2 ∈ positive_nont1_t, n ∈ nonzero_n
            contributions[9] += lag_contribution(data, n_offset, t_offset, n, t1, n, t2)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < 0 && t1 > 0
    for t1 ∈ positive_t, n2 ∈ nonzero_n, t2 ∈ negative_t
        contributions[9] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, t2)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < 0 && t2 > 0
    for n1 ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[9] += lag_contribution(data, n_offset, t_offset, n1, t1, 0, t2)
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < t1 < 0
    for t1 ∈ negative_t[begin+1:end], n2 ∈ nonzero_n
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[9] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, t2)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < t2 < 0
    for n1 ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        for t2 ∈ (t1+1):-1
            contributions[9] += lag_contribution(data, n_offset, t_offset, n1, t1, 0, t2)
        end
    end

    # Class X
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 < t2 < 0
    for n ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        for t2 ∈ (t1+1):-1
            contributions[10] += lag_contribution(data, n_offset, t_offset, n, t1, n, t2)
        end
    end
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t2 < t1 < 0
    for n ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[10] += lag_contribution(data, n_offset, t_offset, n, t1, n, t2)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 > 0 && t1 < 0
    for t1 ∈ negative_t, n2 ∈ nonzero_n, t2 ∈ positive_t
        contributions[10] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, t2)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 > 0 && t2 < 0
    for n1 ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[10] += lag_contribution(data, n_offset, t_offset, n1, t1, 0, t2)
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && 0 < t1 < t2
    for t1 ∈ positive_t[begin:end-1], n2 ∈ nonzero_n
        for t2 ∈ (t1+1):positive_t[end]
            contributions[10] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, t2)
        end
    end
    # Assume left point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && 0 < t2 < t1
    for n1 ∈ nonzero_n, t1 ∈ positive_t[begin+1:end]
        for t2 ∈ positive_t[begin]:(t1-1)
            contributions[10] += lag_contribution(data, n_offset, t_offset, n1, t1, 0, t2)
        end
    end

    # Class XI
    # Assume left point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; 0 < t1 < t2
    for n1 ∈ nonzero_n, t1 ∈ positive_t[begin:end-1]
        for t2 ∈ (t1+1):positive_t[end]
            contributions[11] += lag_contribution(data, n_offset, t_offset, n1, t1, 0, t2)
        end
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; 0 < t2 < t1
    for t1 ∈ positive_t[begin+1:end], n2 ∈ nonzero_n 
        for t2 ∈ 1:(t1-1)
            contributions[11] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, t2)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; t2 < t1 < 0
    for n1 ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        for t2 ∈ negative_t[begin]:(t1-1)
            contributions[11] += lag_contribution(data, n_offset, t_offset, n1, t1, 0, t2)
        end
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; t1 < t2 < 0
    for t1 ∈ negative_t[begin+1:end], n2 ∈ nonzero_n 
        for t2 ∈ (t1+1):-1
            contributions[11] += lag_contribution(data, n_offset, t_offset, 0, t1, n2, t2)
        end
    end
    # Assume middle point (0,0); n1 left
    # n1 == n2 ≠ 0; t1 < 0 < t2
    for n ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[11] += lag_contribution(data, n_offset, t_offset, n, t1, n, t2)
    end
    # Assume middle point (0,0); n2 left
    # n1 == n2 ≠ 0; t2 < 0 < t1
    for n ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[11] += lag_contribution(data, n_offset, t_offset, n, t1, n, t2)
    end

    # Class XII
    # n1 ≠ n2 ≠ 0
    for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 > 0
        for t ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, n_offset, t_offset, n1, t, n2, t)
        end
        # Assume n1 odd: t1 < 0; t2 == 0
        for t1 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, n_offset, t_offset, n1, t1, n2, 0)
        end
        # Assume n2 odd: t2 < 0; t1 == 0
        for t2 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += lag_contribution(data, n_offset, t_offset, n1, 0, n2, t2)
        end
    end

    # Class XIII
    # n1 ≠ n2 ≠ 0
    for n1 ∈ nonzero_n
        nonzero_notn1 = filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 < 0
        for t ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, n_offset, t_offset, n1, t, n2, t)
        end
        # Assume n1 odd: t1 > 0; t2 == 0
        for t1 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, n_offset, t_offset, n1, t1, n2, 0)
        end
        # Assume n2 odd: t2 > 0; t1 == 0
        for t2 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += lag_contribution(data, n_offset, t_offset, n1, 0, n2, t2)
        end
    end

    # Class XIV
    # n1 ≠ n2 ≠ 0
    n1_with_nonzero_filtered = [(n1, filter_element(nonzero_n, n1)) for n1 ∈ nonzero_n]
    t1_with_nonzero_filtered = [(t1, filter_element(nonzero_n, t1)) for t1 ∈ nonzero_n]
    for (n1, nonzero_notn1) ∈ n1_with_nonzero_filtered, (t1, nonzero_nott1) ∈ t1_with_nonzero_filtered
        for n2 ∈ nonzero_notn1, t2 ∈ nonzero_nott1
            contributions[14] += lag_contribution(data, n_offset, t_offset, n1, t1, n2, t2)
        end
    end

    return contributions

end

# FIXME 
# instead of filtering the element, track the el_idx,
# then split inner loop to go begin:el_idx-1 and el_idx+1:end
function filter_element(arr, el)
    filter(≠(el), arr)
end

function lag_sequence_class_contribution_unrolled(tricorr::TripleCorrelation)
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
    for t1 ∈ negative_t[begin+1:end], n2 ∈ nonzero_n 
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
    t1_with_nonzero_filtered = [(t1, filter_element(nonzero_n, t1)) for t1 ∈ nonzero_n]
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

function sequence_class_tricorr!(class_contribution::AbstractVector, src, space_max_lag, time_max_lag, lags_classifier::Function)
    src = parent(src)

    space_lag_range = -(space_max_lag):(space_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    (N_space, N_times) = size(src)
    time_range = (1-minimum(time_lag_range)):(N_times-maximum(time_lag_range))
    space_range = (1-minimum(space_lag_range)):(N_space-maximum(space_lag_range))

    class_contribution .= 0
    for n1 ∈ space_lag_range, n2 ∈ space_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        
        contribution = 0
        # tturbo missing
        for i_neuron ∈ space_range, i_time ∈ time_range
            contribution += src[i_neuron, i_time] * src[i_neuron+n1,i_time+t1] * src[i_neuron+n2,i_time+t2]
        end
        class_contribution[class] += contribution
    end
    class_contribution ./= calculate_scaling_factor(src, (space_max_lag, time_max_lag))
end

function sequence_class_tricorr(src, space_max_lag, time_max_lag)
    N_network_classifications = 14
    network_class_contributions = Array{Float64}(undef, N_network_classifications)
    lags_classifier = lag_motif_sequence_class

    sequence_class_tricorr!(network_class_contributions, src, space_max_lag, time_max_lag, lags_classifier)
end

function sequence_class_tricorr(src::OffsetArray, args...)
    sequence_class_tricorr(parent(src), args...)
end

# Zero padding calculation

function sequence_class_tricorr_zeropad!(class_contribution::AbstractVector, src::AbstractArray{T}, space_max_lag, time_max_lag, lags_classifier::Function) where T
    src = parent(src)
    space_lag_range = -(space_max_lag):(space_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    space_axis, time_axis = axes(src)

    class_contribution .= 0
    for n1 ∈ space_lag_range, n2 ∈ space_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        contribution = 0

        n_min = first(space_axis)-min(0,min(n1,n2))
        n_max = last(space_axis)-max(0,max(n1,n2))
        n_range = n_min:n_max

        t_min = first(time_axis)-min(0,min(t1,t2))
        t_max = last(time_axis)-max(0,max(t1,t2))
        t_range = t_min:t_max

        # tturbo missing
        for i_space ∈ n_range, i_time ∈ t_range
            contribution += src[i_space, i_time] * src[i_space+n1,i_time+t1] * src[i_space+n2,i_time+t2]
        end
        class_contribution[class] += contribution
    end
    class_contribution ./= calculate_scaling_factor_zeropad(src, (space_max_lag, time_max_lag))
end

function sequence_class_tricorr_zeropad(src, space_max_lag, time_max_lag)
    N_network_classifications = 14
    network_class_contributions = Array{Float64}(undef, N_network_classifications)
    lags_classifier = lag_motif_sequence_class

    sequence_class_tricorr_zeropad!(network_class_contributions, src, space_max_lag, time_max_lag, lags_classifier)
    return network_class_contributions
end

function sequence_class_tricorr_zeropad(src::OffsetArray, args...)
    sequence_class_tricorr_zeropad(parent(src), args...)
end


# Helpers



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