function expectation_conditioned_on_spike_count(raster::AbstractArray, boundary::Periodic, lag_extents::NTuple{2})
    expectation_conditioned_on_spike_count(count(raster), size(raster), lag_extents)
end

function expectation_conditioned_on_spike_count(raster::AbstractArray, boundary::PeriodicExtended, lag_extents)
    expectation_conditioned_on_spike_count(count_raster, (size(raster[1:end-1])..., raster[end]-2 * boundary.boundary), lag_extents)
end

function expectation_conditioned_on_spike_count(count::Number, raster_size, lag_extents::NTuple{2})
    n_lag_extent = lag_extents[1]
    t_lag_extent = lag_extents[end]
    NT = prod(raster_size)
    R = count
    DT1 = (t_lag_extent)  # within t distance with 1 spike
    DT2 = (t_lag_extent - 1)  # within t distance with 2 spikes
    DN1 = (n_lag_extent)
    DN2 = (n_lag_extent - 1)
    DOWN_T = floor(Int, t_lag_extent / 2)
    UP_T = ceil(Int, t_lag_extent / 2)
    PS1 = ((R - 1) / (NT - 1))  # probability of spike given prior spike
    PS2 = ((R - 2) / (NT - 2))  # probability of spike given two prior spikes
    C = 3
    [
        R,  # 0
        C * R * DT1 * PS1,  # I
        R * DT1 * PS1 * DT2 * PS2,   # II
        C * R * DN1 * PS1,  # III
        R * DN1 * PS1 * DN2 * PS2,  # IV
        C * R * DT1 * DN1 * PS1,  # V
        C * R * DN1 * PS1 * DT1 * PS2,  # VI
        C * R * DN1 * PS1 * DT1 * PS2,  # VII
        C^2 / 2 * R * DN1 * UP_T * (UP_T - 1) * PS1 * PS2, # VIII
        3.6 * R * DN1 * UP_T * DOWN_T * PS1 * PS2,  # IX; FIXME what
        C^2 / 2 * R * DN1 * DOWN_T * (DOWN_T - 1) * PS1 * PS2, # X; local dynamics precede
        C * R * UP_T * DN1 * PS1 * DN2 * PS2,  # XI
        C * R * DOWN_T * DN1 * PS1 * DN2 * PS2,  # XII
        R * DT1 * DN1 * PS1 * DT2 * DN2 * PS2  # XIII
    ]
end

function expectation_conditioned_on_constituent_parts(actual, count::Number, raster_size, lag_extents::NTuple{2})
    expected = expectation_conditioned_on_spike_count(count, raster_size, lag_extents)
    [
        expected[1],
        expected[2],  # I
        (expected[3] / expected[2]) * actual[2],   # II
        expected[4],  # III
        (expected[5] / expected[4]) * actual[4],  # IV
        expected[6],  # V
        (expected[7] / (expected[2] * expected[4] * expected[6])) * (actual[2] * actual[4] * actual[6]) ,  # VI
        (expected[8] / (expected[2] * expected[4] * expected[6])) * (actual[2] * actual[4] * actual[6]),  # VII
        (expected[9] / (expected[2] * expected[6]^2)) * (actual[2] * actual[6]^2), # VIII
        (expected[10] / (expected[2] * expected[6]^2)) * (actual[2] * actual[6]^2),  # IX; FIXME what
        (expected[11] / (expected[2] * expected[6]^2)) * (actual[2] * actual[6]^2), # X; local dynamics precede
        (expected[12] / (expected[4] * expected[6])) * (actual[4] * actual[6]),  # XI
        (expected[13] / (expected[4] * expected[6])) * (actual[4] * actual[6]),  # XII
        (expected[14] / (expected[6]^3)) * (actual[6]^3)  # XIII
    ]
end
function rate_normed_sequence_classes(raster, boundary, lag_extents)
    # 0 means same as noise
    raw_sequence_classes = sequence_class_tricorr(raster, boundary, lag_extents)
    raw_sequence_classes ./ (expectation_conditioned_on_spike_count(raster, boundary, lag_extents) ./ calculate_scaling_factor(raster, boundary))
end
function constituent_normed_sequence_classes(raster, boundary, lag_extents)
    # 0 means same as noise
    raster_size = get_raster_size(raster, boundary)
    raw_sequence_classes = sequence_class_tricorr(raster, boundary, lag_extents)
    @show raw_sequence_classes
    # @show expectation_conditioned_on_constituent_parts(raw_sequence_classes, count(raster), raster_size, lag_extents)
    # raw_sequence_classes ./ (expectation_conditioned_on_constituent_parts(raw_sequence_classes, count(raster), raster_size, lag_extents) ./ calculate_scaling_factor(raster, boundary))
    @show expectation_conditioned_on_spike_count(count(raster), raster_size, lag_extents)
    raw_sequence_classes ./ (expectation_conditioned_on_spike_count(count(raster), raster_size, lag_extents) ./ calculate_scaling_factor(raster, boundary))
end