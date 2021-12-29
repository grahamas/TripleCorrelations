using LoopVectorization

function _OLD_lag_contribution(data::Matrix, n1::Int, t1::Int, n2::Int, t2::Int)
    contribution = 0
    n_start = max(1 - min(n1, n2), 1); t_start = max(1 - min(t1, t2), 1)
    n_end = min(size(data,1) - max(n1,n2), size(data,1))
    t_end = min(size(data,2) - max(t1,t2), size(data,2))
    @turbo for data_n ∈ n_start:n_end, data_t ∈ t_start:t_end
        contribution += data[data_n, data_t] * 
            data[data_n+n1, data_t+t1] *
            data[data_n+n2, data_t+t2]
    end
    return contribution
end

function _OLD_sequence_class_tricorr_unrolled(data::AbstractArray, n_max_lag, t_max_lag)
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
    contributions[1] = _OLD_lag_contribution(data, 0,0,0,0)

    # Class II
    # n1, n2 == 0, 0
    # t1 == 0 or t2 == 0
    @inbounds for t ∈ nonzero_t
        contributions[2] += _OLD_lag_contribution(data, 0,t,0,0) + _OLD_lag_contribution(data, 0,0,0,t)
    end
    # t1 == t2
    @inbounds for t ∈ nonzero_t
        contributions[2] += _OLD_lag_contribution(data, 0,t,0,t)
    end

    # Class III
    # n1, n2 == 0, 0
    # t1 ≠ t2 ≠ 0
    @inbounds for t1 ∈ nonzero_t
        nonzero_nont1_t = TripleCorrelations.filter_element(nonzero_t, t1)
        @inbounds for t2 ∈ nonzero_nont1_t
            contributions[3] += _OLD_lag_contribution(data, 0,t1,0,t2)
        end
    end

    # Class IV
    # t1, t2 == 0, 0
    # n1 == 0 or n2 == 0
    @inbounds for n ∈ nonzero_n
        contributions[4] += _OLD_lag_contribution(data, n,0,0,0) + _OLD_lag_contribution(data, 0,0,n,0)
    end
    # t1 == t2
    @inbounds for n ∈ nonzero_n
        contributions[4] += _OLD_lag_contribution(data, n,0,n,0)
    end

    # Class V
    # t1, t2 == 0, 0
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_nonn1_n = TripleCorrelations.filter_element(nonzero_n, n1)
        @inbounds for n2 ∈ nonzero_nonn1_n
            contributions[5] += _OLD_lag_contribution(data, n1,0,n2,0)
        end
    end

    # Class VI
    # n1, t1 == 0, 0 && n2, t2 ≠ 0, 0
    @inbounds for n2 ∈ nonzero_n, t2 ∈ nonzero_t
        contributions[6] += _OLD_lag_contribution(data, 0, 0, n2, t2)
    end
    # n2, t2 == 0, 0 && n1, t1 ≠ 0, 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ nonzero_t
        contributions[6] += _OLD_lag_contribution(data, n1, t1, 0, 0)
    end
    # n1, t1 == n2, t2 ≠ 0, 0
    @inbounds for n ∈ nonzero_n, t ∈ nonzero_t
        contributions[6] += _OLD_lag_contribution(data, n, t, n, t)
    end

    # Class VII
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n2 ≠ 0 && n1 == 0
    @inbounds for t ∈ negative_t, n2 ∈ nonzero_n
        contributions[7] += _OLD_lag_contribution(data, 0, t, n2, t)
    end
    # Assuming horz arm (rh) point (0,0)
    # t1 == t2 < 0 && n1 ≠ 0 && n2 == 0
    @inbounds for t ∈ negative_t, n1 ∈ nonzero_n
        contributions[7] += _OLD_lag_contribution(data, n1, t, 0, t)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 > 0 && t1 == 0
    @inbounds for n ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += _OLD_lag_contribution(data, n, 0, n, t2)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 > 0 && t2 == 0
    @inbounds for n ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += _OLD_lag_contribution(data, n, t1, n, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 > 0 && t2 == 0
    @inbounds for n2 ∈ nonzero_n, t1 ∈ positive_t
        contributions[7] += _OLD_lag_contribution(data, 0, t1, n2, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 > 0 && t1 == 0
    @inbounds for n1 ∈ nonzero_n, t2 ∈ positive_t
        contributions[7] += _OLD_lag_contribution(data, n1, 0, 0, t2)
    end

    # Class VIII
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n2 ≠ 0 && n1 == 0
    @inbounds for t ∈ positive_t, n2 ∈ nonzero_n
        contributions[8] += _OLD_lag_contribution(data, 0, t, n2, t)
    end
    # Assuming horz arm (lh) point (0,0)
    # t1 == t2 > 0 && n1 ≠ 0 && n2 == 0
    @inbounds for t ∈ positive_t, n1 ∈ nonzero_n
        contributions[8] += _OLD_lag_contribution(data, n1, t, 0, t)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t2 < 0 && t1 == 0
    @inbounds for n ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += _OLD_lag_contribution(data, n, 0, n, t2)
    end
    # Assuming vert arm (up or down) point (0,0)
    # n1 == n2 ≠ 0 && t1 < 0 && t2 == 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += _OLD_lag_contribution(data, n, t1, n, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n1 == 0 && t1 < 0 && t2 == 0
    @inbounds for n2 ∈ nonzero_n, t1 ∈ negative_t
        contributions[8] += _OLD_lag_contribution(data, 0, t1, n2, 0)
    end
    # Assuming corner point (0,0)
    # n1 ≠ n2 && n2 == 0 && t2 < 0 && t1 == 0
    @inbounds for n1 ∈ nonzero_n, t2 ∈ negative_t
        contributions[8] += _OLD_lag_contribution(data, n1, 0, 0, t2)
    end

    # Class IX
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 ≠ t2 > 0
    @inbounds for t1 ∈ positive_t
        positive_nont1_t = TripleCorrelations.filter_element(positive_t, t1)
        # n in inner loop bc filter_element allocates
        @inbounds for t2 ∈ positive_nont1_t, n ∈ nonzero_n
            contributions[9] += _OLD_lag_contribution(data, n, t1, n, t2)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < 0 && t1 > 0
    @inbounds for t1 ∈ positive_t, n2 ∈ nonzero_n, t2 ∈ negative_t
        contributions[9] += _OLD_lag_contribution(data, 0, t1, n2, t2)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < 0 && t2 > 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[9] += _OLD_lag_contribution(data, n1, t1, 0, t2)
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 < t1 < 0
    @inbounds for t1 ∈ negative_t[begin+1:end], n2 ∈ nonzero_n
        @inbounds for t2 ∈ negative_t[begin]:(t1-1)
            contributions[9] += _OLD_lag_contribution(data, 0, t1, n2, t2)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 < t2 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        @inbounds for t2 ∈ (t1+1):-1
            contributions[9] += _OLD_lag_contribution(data, n1, t1, 0, t2)
        end
    end

    # Class X
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t1 < t2 < 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t[begin:end-1]
        @inbounds for t2 ∈ (t1+1):-1
            contributions[10] += _OLD_lag_contribution(data, n, t1, n, t2)
        end
    end
    # Assume odd point (0,0)
    # n1 == n2 ≠ 0 && t2 < t1 < 0
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        @inbounds for t2 ∈ negative_t[begin]:(t1-1)
            contributions[10] += _OLD_lag_contribution(data, n, t1, n, t2)
        end
    end
    # Assume middle point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && t2 > 0 && t1 < 0
    @inbounds for t1 ∈ negative_t, n2 ∈ nonzero_n, t2 ∈ positive_t
        contributions[10] += _OLD_lag_contribution(data, 0, t1, n2, t2)
    end
    # Assume middle point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && t1 > 0 && t2 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[10] += _OLD_lag_contribution(data, n1, t1, 0, t2)
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0 && n2 ≠ 0 && 0 < t1 < t2
    @inbounds for t1 ∈ positive_t[begin:end-1], n2 ∈ nonzero_n
        @inbounds for t2 ∈ (t1+1):positive_t[end]
            contributions[10] += _OLD_lag_contribution(data, 0, t1, n2, t2)
        end
    end
    # Assume left point (0,0); n1 odd
    # n2 == 0 && n1 ≠ 0 && 0 < t2 < t1
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t[begin+1:end]
        @inbounds for t2 ∈ positive_t[begin]:(t1-1)
            contributions[10] += _OLD_lag_contribution(data, n1, t1, 0, t2)
        end
    end

    # Class XI
    # Assume left point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; 0 < t1 < t2
    @inbounds for n1 ∈ nonzero_n, t1 ∈ positive_t[begin:end-1]
        @inbounds for t2 ∈ (t1+1):positive_t[end]
            contributions[11] += _OLD_lag_contribution(data, n1, t1, 0, t2)
        end
    end
    # Assume left point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; 0 < t2 < t1
    @inbounds for t1 ∈ positive_t[begin+1:end], n2 ∈ nonzero_n 
        @inbounds for t2 ∈ 1:(t1-1)
            contributions[11] += _OLD_lag_contribution(data, 0, t1, n2, t2)
        end
    end
    # Assume right point (0,0); n1 odd
    # n2 == 0; n1 ≠ 0; t2 < t1 < 0
    @inbounds for n1 ∈ nonzero_n, t1 ∈ negative_t[begin+1:end]
        @inbounds for t2 ∈ negative_t[begin]:(t1-1)
            contributions[11] += _OLD_lag_contribution(data, n1, t1, 0, t2)
        end
    end
    # Assume right point (0,0); n2 odd
    # n1 == 0; n2 ≠ 0; t1 < t2 < 0
    @inbounds for t1 ∈ negative_t[begin:end-1], n2 ∈ nonzero_n 
        @inbounds for t2 ∈ (t1+1):-1
            contributions[11] += _OLD_lag_contribution(data, 0, t1, n2, t2)
        end
    end
    # Assume middle point (0,0); n1 left
    # n1 == n2 ≠ 0; t1 < 0 < t2
    @inbounds for n ∈ nonzero_n, t1 ∈ negative_t, t2 ∈ positive_t
        contributions[11] += _OLD_lag_contribution(data, n, t1, n, t2)
    end
    # Assume middle point (0,0); n2 left
    # n1 == n2 ≠ 0; t2 < 0 < t1
    @inbounds for n ∈ nonzero_n, t1 ∈ positive_t, t2 ∈ negative_t
        contributions[11] += _OLD_lag_contribution(data, n, t1, n, t2)
    end

    # Class XII
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_notn1 = TripleCorrelations.filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 > 0
        @inbounds for t ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[12] += _OLD_lag_contribution(data, n1, t, n2, t)
        end
        # Assume n1 odd: t1 < 0; t2 == 0
        @inbounds for t1 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += _OLD_lag_contribution(data, n1, t1, n2, 0)
        end
        # Assume n2 odd: t2 < 0; t1 == 0
        @inbounds for t2 ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[12] += _OLD_lag_contribution(data, n1, 0, n2, t2)
        end
    end

    # Class XIII
    # n1 ≠ n2 ≠ 0
    @inbounds for n1 ∈ nonzero_n
        nonzero_notn1 = TripleCorrelations.filter_element(nonzero_n, n1)
        # Assume (0,0) odd: t1 == t2 < 0
        @inbounds for t ∈ negative_t, n2 ∈ nonzero_notn1
            contributions[13] += _OLD_lag_contribution(data, n1, t, n2, t)
        end
        # Assume n1 odd: t1 > 0; t2 == 0
        @inbounds for t1 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += _OLD_lag_contribution(data, n1, t1, n2, 0)
        end
        # Assume n2 odd: t2 > 0; t1 == 0
        @inbounds for t2 ∈ positive_t, n2 ∈ nonzero_notn1
            contributions[13] += _OLD_lag_contribution(data, n1, 0, n2, t2)
        end
    end

    # Class XIV
    # n1 ≠ n2 ≠ 0
    n1_with_nonzero_filtered = [(n1, TripleCorrelations.filter_element(nonzero_n, n1)) for n1 ∈ nonzero_n]
    t1_with_nonzero_filtered = [(t1, TripleCorrelations.filter_element(nonzero_t, t1)) for t1 ∈ nonzero_t]
    @inbounds for (n1, nonzero_notn1) ∈ n1_with_nonzero_filtered, (t1, nonzero_nott1) ∈ t1_with_nonzero_filtered
        @inbounds for n2 ∈ nonzero_notn1, t2 ∈ nonzero_nott1
            contributions[14] += _OLD_lag_contribution(data, n1, t1, n2, t2)
        end
    end

    return contributions ./ TripleCorrelations.calculate_scaling_factor_zeropad(data)

end
