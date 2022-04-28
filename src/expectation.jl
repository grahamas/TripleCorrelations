
function expectation_conditioned_on_spike_count(count, raster_size, lag_extents::NTuple{2})
    n_lag_extent = lag_extents[1]
    t_lag_extent = lag_extents[end]
    NT = prod(raster_size)
    R = count
    DT1 = (t_lag_extent)  # within t distance with 1 spike
    DT2 = (t_lag_extent - 1)  # within t distance with 2 spikes
    DN1 = (n_lag_extent)
    DN2 = (n_lag_extent - 1)
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
        C^2 / 2 * R * DN1 * ceil(Int, t_lag_extent / 2) * (ceil(Int, t_lag_extent / 2) - 1) * PS1 * PS2, # VIII
        3.6 * R * DN1 * ceil(Int, t_lag_extent / 2) * floor(Int, t_lag_extent / 2) * PS1 * PS2,  # IX; FIXME what
        C^2 / 2 * R * DN1 * floor(Int, t_lag_extent / 2) * (floor(Int, t_lag_extent / 2) - 1) * PS1 * PS2, # X; local dynamics precede
        C * R * ceil(Int, t_lag_extent / 2) * DN1 * PS1 * DN2 * PS2,  # XI
        C * R * floor(Int, t_lag_extent / 2) * DN1 * PS1 * DN2 * PS2,  # XII
        R * DT1 * DN1 * PS1 * DT2 * DN2 * PS2  # XIII
    ]
end

function rate_normed_sequence_classes(raster, boundary, lag_extents)
    # 0 means same as noise
    raster_size = get_raster_size(raster, boundary)
    raw_sequence_classes = sequence_class_tricorr(raster, boundary, lag_extents)
    (raw_sequence_classes ./ expectation_conditioned_on_spike_count(count(raster), raster_size, lag_extents) ./ calculate_scaling_factor(raster, boundary)) .- 1
end