

abstract type AbstractDistributionAssumption end
struct IndBernoulli <: AbstractDistributionAssumption end
struct IndStdNormal <: AbstractDistributionAssumption end

abstract type AbstractConditional end
struct None <: AbstractConditional end
struct Rate <: AbstractConditional end
struct Constituents <: AbstractConditional end

export IndBernoulli, IndStdNormal
export None, Rate, Constituents


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

function estimate_μ(assumption::IndBernoulli, condition::Rate, spike_count::Int, raster_size::Tuple, boundary::Periodic, lag_extents)
    p = spike_count / prod(raster_size)
    return prod(raster_size) .* (p .^ motif_order()) .* 
        triplet_count_per_motif_base_node(boundary, lag_extents) ./ 
        calculate_scaling_factor(raster_size, boundary)
end

function estimate_μ(assumption::IndBernoulli, condition::Rate, raster::Matrix{Bool}, boundary::Periodic, lag_extents)
    estimate_μ(assumption, condition, count(raster), size_meat(raster), boundary, lag_extents)
end

function estimate_μ(assumption::IndBernoulli, condition::Rate, raster::Matrix{Bool}, boundary::PeriodicExtended, lag_extents)
    p = mean(raster)
    return prod(size_meat(raster,boundary)) .* (p .^ motif_order()) .* 
        triplet_count_per_motif_base_node(boundary, lag_extents) ./ 
        calculate_scaling_factor(raster, boundary)
end

function estimate_σ(assumption::IndBernoulli, condition::Rate, raster::Matrix, boundary::AbstractBoundaryCondition, lag_extents, n_bootstraps=100)
    estimate_σ(assumption, condition, n_spikes, size(raster), boundary, lag_extents, 100)
end

function estimate_σ(assumption::IndBernoulli, condition::Rate, n_spikes::Int, raster_size::Tuple, boundary::AbstractBoundaryCondition, lag_extents, n_bootstraps=100)
    raster = zeros(Bool, raster_size...)
    raster[1:n_spikes] .= true
    l_contributions = [sequence_class_tricorr(shuffle!(raster), boundary, lag_extents) for _ ∈ 1:(n_bootstraps-1)]
    std(l_contributions)
end

function estimate_σ(::IndStdNormal, condition::None, raster::Matrix, boundary::AbstractBoundaryCondition, lag_extents, n_bootstraps=100)
    l_contributions = [sequence_class_tricorr(randn(size(raster)...), boundary, lag_extents) for _ ∈ 1:(n_bootstraps-1)]
    std(l_contributions)
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
