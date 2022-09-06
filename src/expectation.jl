

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
        TripleCorrelations.calculate_scaling_factor(raster_size, boundary)
end

function estimate_μ(assumption::IndBernoulli, condition::Rate, raster::AbstractMatrix{Bool}, boundary::Periodic, lag_extents)
    estimate_μ(assumption, condition, count(raster), size_meat(raster), boundary, lag_extents)
end

function estimate_μ(assumption::IndBernoulli, condition::Rate, raster::AbstractMatrix{Bool}, boundary::PeriodicExtended, lag_extents)
    p = mean(raster)
    return prod(size_meat(raster,boundary)) .* (p .^ motif_order()) .* 
        triplet_count_per_motif_base_node(boundary, lag_extents) ./ 
        TripleCorrelations.calculate_scaling_factor(raster, boundary)
end

function estimate_σ(assumption::IndBernoulli, condition::Rate, raster::AbstractMatrix, boundary::AbstractBoundaryCondition, lag_extents, n_bootstraps=100)
    estimate_σ(assumption, condition, n_spikes, size(raster), boundary, lag_extents, 100)
end

function estimate_σ(assumption::IndBernoulli, condition::Rate, n_spikes::Int, raster_size::Tuple, boundary::AbstractBoundaryCondition, lag_extents, n_bootstraps=100)
    raster = zeros(Bool, raster_size...)
    raster[1:n_spikes] .= true
    l_contributions = [sequence_class_tricorr(shuffle!(raster), boundary, lag_extents) for _ ∈ 1:(n_bootstraps-1)]
    std(l_contributions)
end

function estimate_σ(::IndStdNormal, condition::None, raster::AbstractMatrix, boundary::AbstractBoundaryCondition, lag_extents, n_bootstraps=100)
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


#################### Constituent


function estimate_μ(assumption::IndBernoulli, condition::Constituents, count::Int, raster_size, boundary, lag_extents, measured)
    rate_expectation = estimate_μ(assumption, Rate(), count, raster_size, boundary, lag_extents)
    [
        rate_expectation[1],
        rate_expectation[2],  # I
        (rate_expectation[3] / (rate_expectation[2] * rate_expectation[1])) * (measured[2] * measured[1]),   # II
        rate_expectation[4],  # III
        (rate_expectation[5] / (rate_expectation[4] * rate_expectation[1])) * (measured[4] * measured[1]),  # IV
        rate_expectation[6],  # V
        (rate_expectation[7] / sqrt(rate_expectation[2] * rate_expectation[4] * rate_expectation[6])) * sqrt(measured[2] * measured[4] * measured[6]),  # VI
        (rate_expectation[8] / sqrt(rate_expectation[2] * rate_expectation[4] * rate_expectation[6])) * sqrt(measured[2] * measured[4] * measured[6]),  # VII
        (rate_expectation[9] / sqrt(rate_expectation[2] * rate_expectation[6]^2)) * sqrt(measured[2] * measured[6]^2), # VIII
        (rate_expectation[10] / sqrt(rate_expectation[2] * rate_expectation[6]^2)) * sqrt(measured[2] * measured[6]^2),  # IX; FIXME what
        (rate_expectation[11] / sqrt(rate_expectation[2] * rate_expectation[6]^2)) * sqrt(measured[2] * measured[6]^2), # X; local dynamics precede
        (rate_expectation[12] / sqrt(rate_expectation[4] * rate_expectation[6]^2)) * sqrt(measured[4] * measured[6]^2),  # XI
        (rate_expectation[13] / sqrt(rate_expectation[4] * rate_expectation[6]^2)) * sqrt(measured[4] * measured[6]^2),  # XII
        (rate_expectation[14] / (rate_expectation[6]^(3/2))) * (measured[6]^(3/2))  # XIII
    ]

end


function expectation_conditioned_on_constituent_parts(actual, raster, boundary::AbstractBoundaryCondition, lag_extents::NTuple{2})
    expected = expectation_conditioned_on_spike_count(raster, boundary, lag_extents)
    
end
