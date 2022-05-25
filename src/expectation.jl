function slice_meat(raster, boundary::PeriodicExtended)
    bd = boundary.boundary
    fin = size(raster)[end]
    view_slice_last(raster, (bd+1):(fin-bd))
end

function size_meat(raster, boundary::PeriodicExtended)
    (size(raster[1:end-1])..., raster[end]-2 * boundary.boundary)
end

function triplet_count_per_motif_base_node(boundary::Union{Periodic,PeriodicExtended}, lag_extents)
    t_pm = lag_extents[end]
    n_extents = lag_extents[1:end-1]
    n_pm = prod(n_extents .+ 1) - 1

    t_m = floor(t_pm / 2)
    t_p = ceil(t_pm / 2)

    [
        1  # 0
        3*t_pm  # I
        t_pm*(t_pm-1)  # II
        3*n_pm  # III
        n_pm*(n_pm-1)  # IV
        3*n_pm*t_pm  # V
        4*n_pm*t_p + 2*n_pm*t_m  # VI
        4*n_pm*t_m + 2*n_pm*t_p  # VII
        n_pm*t_p*(t_p-1) + 2*n_pm*t_m*t_p + n_pm*t_m*(t_m-1) # VIII
        n_pm*(t_p)*(t_p-1) + n_pm*(t_m)*(t_m-1) + 2n_pm*t_m*t_p  # IX
        n_pm*t_m*(t_m-1) + 2*n_pm*t_p*t_m + n_pm*t_p*(t_p-1)  # X
        n_pm*(n_pm-1)*t_p + 2n_pm*(n_pm-1)*t_m
        n_pm*(n_pm-1)*t_m + 2n_pm*(n_pm-1)*t_p
        n_pm*(n_pm-1)*t_pm*(t_pm-1)
    ]
end

function motif_order()
    [
        1 # 0
        2 # 1
        3 # 2
        2 # 3
        3 # IV
        2 # V
        3
        3
        3
        3
        3
        3
        3
        3
    ]
end

function expectation_of_independent_spiking_conditioned_on_rate(raster::Matrix{Bool}, boundary::Periodic, lag_extents)
    p = mean(raster)
    return prod(size(raster)) .* (p .^ motif_order()) .* 
        triplet_count_per_motif_base_node(boundary, lag_extents) ./ 
        calculate_scaling_factor(raster, boundary)
end

function expectation_of_independent_spiking_conditioned_on_rate(raster::Matrix{Bool}, boundary::PeriodicExtended, lag_extents)
    p = mean(raster)
    return prod(size_meat(raster,boundary)) .* (p .^ motif_order()) .* 
        triplet_count_per_motif_base_node(boundary, lag_extents) ./ 
        calculate_scaling_factor(raster, boundary)
end


function sequence_classes_divide_E_given_rate(raster, boundary, lag_extents)
    # 0 means same as noise
    raw_sequence_classes = sequence_class_tricorr(raster, boundary, lag_extents)
    raw_sequence_classes ./= expectation_of_independent_spiking_conditioned_on_rate(raster, boundary, lag_extents)
end
function sequence_classes_divide_E_given_constituents(raster, boundary, lag_extents)
    # 0 means same as noise
    raw_sequence_classes = sequence_class_tricorr(raster, boundary, lag_extents)
    raw_sequence_classes ./ expectation_conditioned_on_constituent_parts(raw_sequence_classes, raster, boundary, lag_extents)
end

function variance_of_standard_normals(boundary::Periodic, lag_extents)
    t_pm = lag_extents[end]
    n_extents = lag_extents[1:end-1]
    n_pm = prod(n_extents .+ 1) - 1

    t_m = floor(t_pm / 2)
    t_p = ceil(t_pm / 2)

    counts = triplet_count_per_motif_base_node(boundary, lag_extents)
    counts_coefficients = [
        15 # 0
        3 # I
        2 # II
        3 # III
        2 # IV
        3 # V
        2
        2
        2
        2
        2
        2
        2
        2
    ]
    bias = [
        0
        6t_pm
        0
        6n_pm
        0
        3n_pm*t_pm + n_pm*t_pm*(n_pm-1)*(t_pm-1)
        0
        0
        0
        0
        0
        0
        0
        0
    ]
    (counts .* counts_coefficients) .+ bias
end
