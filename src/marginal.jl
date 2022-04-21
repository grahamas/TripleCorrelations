

function time_tricorr_zeropad!(time_contribution_offset::OffsetMatrix, src::AbstractArray{T}, max_lags) where T
    time_contribution = parent(time_contribution_offset)
    space_lag_range = -(space_max_lag):(space_max_lag)    
    time_lag_range = -(time_max_lag):(time_max_lag)

    space_axis, time_axis = axes(src)
    
    padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag)) |> collect

    padded_space_axis = space_axis .+ space_max_lag
    padded_time_axis = time_axis .+ time_max_lag

    # padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag))

    # padded_space_axis = space_axis
    # padded_time_axis = time_axis

    time_contribution .= 0.
    #@warn "No turbo."
    @tturbo for n1i ∈ eachindex(space_lag_range), 
            n2i ∈ eachindex(space_lag_range), 
            t1i ∈ eachindex(time_lag_range), 
            t2i ∈ eachindex(time_lag_range)
        n1 = space_lag_range[n1i]; n2 = space_lag_range[n2i] 
        t1 = time_lag_range[t1i]; t2 = time_lag_range[t2i]
        contribution = time_contribution[t1i, t2i]
        for i_neuron ∈ padded_space_axis, i_time ∈ padded_time_axis
            contribution += padded_src[i_neuron, i_time] * padded_src[i_neuron+n1,i_time+t1] * padded_src[i_neuron+n2,i_time+t2]
        end
        time_contribution[t1i, t2i] = contribution
    end
    time_contribution ./ calculate_scaling_factor(src, (max_lags))
    return
end

function time_tricorr_zeropad(src, max_lags)
    time_lag_range = -(time_max_lag):time_max_lag
    time_contributions = OffsetArray{Float64}(undef, time_lag_range, time_lag_range)
    time_tricorr_zeropad!(time_contributions, src, max_lags)
    return time_contributions
end

function time_tricorr_zeropad(src::OffsetArray, args...)
    time_tricorr_zeropad(parent(src), args...)
end

#### Space time

function space_time_tricorr_zeropad!(space_time_contribution_offset::OffsetMatrix, src::AbstractArray{T}, max_lags) where T
    space_time_contribution = parent(space_time_contribution_offset)
    space_lag_range = -(space_max_lag):(space_max_lag)    
    time_lag_range = -(time_max_lag):(time_max_lag)

    space_axis, time_axis = axes(src)
    
    padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag)) |> collect

    padded_space_axis = space_axis .+ space_max_lag
    padded_time_axis = time_axis .+ time_max_lag

    # padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag))

    # padded_space_axis = space_axis
    # padded_time_axis = time_axis

    all_contribution = 0.
    space_time_contribution .= 0.
    #@warn "No turbo."
    @tturbo for n1i ∈ eachindex(space_lag_range), 
            n2i ∈ eachindex(space_lag_range), 
            t1i ∈ eachindex(time_lag_range), 
            t2i ∈ eachindex(time_lag_range)
        n1 = space_lag_range[n1i]; n2 = space_lag_range[n2i] 
        t1 = time_lag_range[t1i]; t2 = time_lag_range[t2i]
        contribution = space_time_contribution[n1i, t1i]
        for i_neuron ∈ padded_space_axis, i_time ∈ padded_time_axis
            contribution += padded_src[i_neuron, i_time] * padded_src[i_neuron+n1,i_time+t1] * padded_src[i_neuron+n2,i_time+t2]
        end
        space_time_contribution[n1i, t1i] = contribution
    end
    space_time_contribution ./= calculate_scaling_factor(src, (max_lags))
    return
end

function space_time_tricorr_zeropad(src, max_lags)
    time_lag_range = -(time_max_lag):time_max_lag
    space_lag_range = -(space_max_lag):space_max_lag
    space_time_contributions = OffsetArray{Float64}(undef, space_lag_range, time_lag_range)
    space_time_tricorr_zeropad!(space_time_contributions, src, max_lags)
    return space_time_contributions
end

function space_time_tricorr_zeropad(src::OffsetArray, args...)
    space_time_tricorr_zeropad(parent(src), args...)
end


###### Space

function space_tricorr_zeropad!(space_contribution_offset::OffsetMatrix, src::AbstractArray{T}, max_lags) where T
    space_contribution = parent(space_contribution_offset)
    space_lag_range = -(space_max_lag):(space_max_lag)    
    time_lag_range = -(time_max_lag):(time_max_lag)

    space_axis, time_axis = axes(src)
    
    padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag)) |> collect

    padded_space_axis = space_axis .+ space_max_lag
    padded_time_axis = time_axis .+ time_max_lag

    # padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag))

    # padded_space_axis = space_axis
    # padded_time_axis = time_axis

    all_contribution = 0.
    space_contribution .= 0.
    #@warn "No turbo."
    @tturbo for n1i ∈ eachindex(space_lag_range), 
            n2i ∈ eachindex(space_lag_range), 
            t1i ∈ eachindex(time_lag_range), 
            t2i ∈ eachindex(time_lag_range)
        n1 = space_lag_range[n1i]; n2 = space_lag_range[n2i] 
        t1 = time_lag_range[t1i]; t2 = time_lag_range[t2i]
        contribution = space_contribution[n1i, n2i]
        for i_neuron ∈ padded_space_axis, i_time ∈ padded_time_axis
            contribution += padded_src[i_neuron, i_time] * padded_src[i_neuron+n1,i_time+t1] * padded_src[i_neuron+n2,i_time+t2]
        end
        space_contribution[n1i, n2i] = contribution
    end
    space_contribution ./= calculate_scaling_factor(src, (max_lags))
    return
end

function space_tricorr_zeropad(src, max_lags)
    space_lag_range = -(space_max_lag):space_max_lag
    space_contributions = OffsetArray{Float64}(undef, space_lag_range, space_lag_range)
    space_tricorr_zeropad!(space_contributions, src, max_lags)
    return space_contributions
end

function space_tricorr_zeropad(src::OffsetArray, args...)
    space_tricorr_zeropad(parent(src), args...)
end


#### Marginal

function marginal_tricorr_zeropad!(marginal_contributions_offset::NamedTuple, src::AbstractArray{T}, max_lags) where T
    
    time_contribution_offset = marginal_contributions_offset.time
    space_time_contribution_offset = marginal_contributions_offset.space_time
    space_contribution_offset = marginal_contributions_offset.space
    
    time_contribution = parent(time_contribution_offset)
    space_time_contribution = parent(space_time_contribution_offset)
    space_contribution = parent(space_contribution_offset)

    space_lag_range = -(space_max_lag):(space_max_lag)    
    time_lag_range = -(time_max_lag):(time_max_lag)

    space_axis, time_axis = axes(src)
    
    padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag)) |> collect

    padded_space_axis = space_axis .+ space_max_lag
    padded_time_axis = time_axis .+ time_max_lag

    # padded_src = PaddedView(zero(T), src, (first(space_axis)-space_max_lag:last(space_axis)+space_max_lag, first(time_axis)-time_max_lag:last(time_axis)+time_max_lag))

    # padded_space_axis = space_axis
    # padded_time_axis = time_axis

    time_contribution .= 0.
    space_time_contribution .= 0.
    space_contribution .= 0.
    #@warn "No turbo."
    @tturbo for n1i ∈ eachindex(space_lag_range), 
            n2i ∈ eachindex(space_lag_range), 
            t1i ∈ eachindex(time_lag_range), 
            t2i ∈ eachindex(time_lag_range)
        n1 = space_lag_range[n1i]; n2 = space_lag_range[n2i] 
        t1 = time_lag_range[t1i]; t2 = time_lag_range[t2i]
        contribution = 0
        for i_space ∈ padded_space_axis, i_time ∈ padded_time_axis
            contribution += padded_src[i_space, i_time] * padded_src[i_space+n1,i_time+t1] * padded_src[i_space+n2,i_time+t2]
        end
        time_contribution[t1i, t2i] += contribution
        space_time_contribution[n1i, t1i] += contribution
        space_contribution[n1i, n2i] += contribution
    end
    scaling = calculate_scaling_factor(src, (max_lags))
    time_contribution ./= scaling
    space_time_contribution ./= scaling
    space_contribution ./= scaling
    return
end

function marginal_tricorr_zeropad(src, max_lags)
    time_lag_range = -(time_max_lag):time_max_lag
    space_lag_range = -(space_max_lag):space_max_lag
    space_time_contributions = OffsetArray{Float64}(undef, space_lag_range, time_lag_range)
    time_contributions = OffsetArray{Float64}(undef, time_lag_range, time_lag_range)
    space_contributions = OffsetArray{Float64}(undef, space_lag_range, space_lag_range)
    marginal_contributions = (time=time_contributions, space=space_contributions, space_time=space_time_contributions)
    marginal_tricorr_zeropad!(marginal_contributions, src, max_lags)
    return marginal_contributions
end

# function marginal_tricorr_zeropad(src::OffsetArray, args...)
#     marginal_tricorr_zeropad(parent(src), args...)
# end