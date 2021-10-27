# Truncating calculation

function class_tricorr!(class_contribution::Vector, src, space_max_lag, time_max_lag, lags_classifier::Function)
    space_lag_range = -(space_max_lag):(space_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    (N_space, N_times) = size(src)
    time_range = (1-minimum(time_lag_range)):(N_times-maximum(time_lag_range))
    space_range = (1-minimum(space_lag_range)):(N_space-maximum(space_lag_range))

    class_contribution .= 0
    @tturbo for n1 ∈ space_lag_range, n2 ∈ space_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        contribution = 0
        for i_neuron ∈ space_range, i_time ∈ time_range
            contribution += src[i_neuron, i_time] * src[i_neuron+n1,i_time+t1] * src[i_neuron+n2,i_time+t2]
        end
        class_contribution[class] += contribution
    end
    return class_contribution ./ calculate_scaling_factor(src, (space_max_lag, time_max_lag))
end

function sequence_class_tricorr(src, space_max_lag, time_max_lag)
    N_network_classifications = 14
    network_class_contributions = Array{Float64}(undef, N_network_classifications)
    lags_classifier = lag_motif_sequence_class

    class_tricorr!(network_class_contributions, src, space_max_lag, time_max_lag, lags_classifier)
end

function sequence_class_tricorr(src::OffsetArray, args...)
    sequence_class_tricorr(parent(src), args...)
end

# Zero padding calculation

function class_tricorr_zeropad!(class_contribution::Vector, src::AbstractArray{T}, space_max_lag, time_max_lag, lags_classifier::Function) where T
    space_lag_range = -(space_max_lag):(space_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    space_axis, time_axis = axes(src)
    
    padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag)) |> collect

    padded_space_axis = space_axis .+ space_max_lag
    padded_time_axis = time_axis .+ time_max_lag

    class_contribution .= 0
    for n1 ∈ space_lag_range, n2 ∈ space_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        contribution = 0
        @tturbo for i_space ∈ padded_space_axis, i_time ∈ padded_time_axis
            contribution += padded_src[i_space, i_time] * padded_src[i_space+n1,i_time+t1] * padded_src[i_space+n2,i_time+t2]
        end
        class_contribution[class] += contribution
    end
    return class_contribution ./ calculate_scaling_factor(src, (space_max_lag, time_max_lag))
end

function sequence_class_tricorr_zeropad(src, space_max_lag, time_max_lag)
    N_network_classifications = 14
    network_class_contributions = Array{Float64}(undef, N_network_classifications)
    lags_classifier = lag_motif_sequence_class

    class_tricorr_zeropad!(network_class_contributions, src, space_max_lag, time_max_lag, lags_classifier)
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