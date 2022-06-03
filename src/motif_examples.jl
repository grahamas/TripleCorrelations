motif_examples = Dict(
    "0" => [0 0 0;
            0 1 0;
            0 0 0],
    "I" => [0 0 0;
             1 1 0;
             0 0 0],
    "II" => [0 0 0;
              1 1 1;
              0 0 0],
    "III" => [0 1 0;
             0 1 0;
             0 0 0],
    "IV" => [0 1 0;
            0 1 0;
            0 1 0],
    "V" => [1 0 0;
             0 1 0;
             0 0 0],
    "VI" => [1 0 0;
              1 1 0;
              0 0 0],
    "VII" => [1 1 0;
             0 1 0;
             0 0 0],
    "VIII" => [0 1 1;
             1 0 0;
             0 0 0],
    "X" => [1 1 0;
            0 0 1;
            0 0 0],
    "IX" => [1 0 1;
             0 1 0;
             0 0 0],
    "XI" => [1 0 0;
              0 1 0;
              0 1 0],
    "XII" => [1 0 0;
               0 1 0;
               1 0 0],
    "XIII" => [1 0 0;
               0 1 0;
               0 0 1]
)

function rand_0(n_range, t_range, n_jitter, t_jitter)
    point = (rand(n_range), rand(t_range))
    (point, point, point)
end

function rand_I(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; t1_lag_range=[-t_jitter:-1..., 1:t_jitter...])
end

function rand_II(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        t1_lag_range=[-t_jitter:-1..., 1:t_jitter...], 
        t2_lag_range=[-t_jitter:-1..., 1:t_jitter...])
end

function rand_III(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...])
end

function rand_IV(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        n2_lag_range=[-n_jitter:-1..., 1:n_jitter...])
end

function rand_V(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t1_lag_range=[-t_jitter:-1..., 1:t_jitter...])
end

function rand_VI(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t2_lag_range=[1:t_jitter...])
end

function rand_VII(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t2_lag_range=[-t_jitter:-1...])
end

function rand_VIII(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t1_lag_range=[-t_jitter:-1...],
        t2_lag_range=[1:t_jitter...])
end

function rand_X(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t1_lag_range=[1:t_jitter...],
        t2_lag_range=[-t_jitter:-1...])
end

function rand_IX(n_range, t_range, n_jitter, t_jitter)
    (n0, t0) = (rand(n_range), rand(t_range))
    n1_lag = 0
    t1_lag = rand([-t_jitter:-2...,2:t_jitter...])
    n2_lag = rand([-n_jitter:-1...,1:n_jitter...])
    t1_sgn = sign(t1_lag) * 1
    t2_lag = rand(t1_sgn:t1_sgn:(t1_lag-t1_sgn))
    (n1, t1) = (n0 + n1_lag, t0 + t1_lag)
    (n2, t2) = (n0 + n2_lag, t0 + t2_lag)
    return ((n0, t0), (n1, t1), (n2, t2))
end

function rand_XI(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        n2_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t2_lag_range=[-t_jitter:-1...])
end

function rand_XII(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        n2_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t2_lag_range=[1:t_jitter...])
end

function rand_XIII(n_range, t_range, n_jitter, t_jitter)
    rand_motif(n_range, t_range; 
        n1_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        n2_lag_range=[-n_jitter:-1..., 1:n_jitter...],
        t1_lag_range=[-t_jitter:-1..., 1:t_jitter...],
        t2_lag_range=[-t_jitter:-1..., 1:t_jitter...])
end

function rand_motif(n_range, t_range; n1_lag_range=0:0, t1_lag_range=0:0, n2_lag_range=0:0, t2_lag_range=0:0)
    # any equality will be with node0
    # all jitters are wrt node0
    (n0, t0) = (rand(n_range), rand(t_range))
    n1_lag = rand(n1_lag_range); t1_lag = rand(t1_lag_range)
    n2_lag = if n2_lag_range != 0:0 
        rand(n2_lag_range[n2_lag_range .!= n1_lag])
    else
        0
    end
    t2_lag = if t2_lag_range != 0:0
        rand(t2_lag_range[t2_lag_range .!= t1_lag])
    else
        0
    end
    (n1, t1) = (n0 + n1_lag, t0 + t1_lag)
    (n2, t2) = (n0 + n2_lag, t0 + t2_lag)
    return ((n0, t0), (n1, t1), (n2, t2))
end

function offset_motif_numeral(n::Integer)
    roman_encode_zero(n-1)
end

# https://www.rosettacode.org/wiki/Roman_numerals/Encode#Julia
function roman_encode_zero(n::Integer)
    if n == 0 return "0" end
    if n < 0 || n > 4999 throw(DomainError(n)) end
 
    DR = [["I", "X", "C", "M"] ["V", "L", "D", "MMM"]]
    rnum = ""
    for (omag, d) in enumerate(digits(n))
        if d == 0
            omr = ""
        elseif d <  4
            omr = DR[omag, 1] ^ d
        elseif d == 4
            omr = DR[omag, 1] * DR[omag, 2]
        elseif d == 5
            omr = DR[omag, 2]
        elseif d <  9
            omr = DR[omag, 2] * DR[omag, 1] ^ (d - 5)
        else
            omr = DR[omag, 1] * DR[omag + 1, 1]
        end
        rnum = omr * rnum
    end
    return rnum
end

function rand_motif(motif_class::String, base_node_ranges::Tuple, lag_jitters::Tuple)
    n_range, t_range = base_node_ranges
    n_jitter, t_jitter = lag_jitters
    rand_motif(motif_class, n_range, t_range, n_jitter, t_jitter)
end

function rand_motif(motif_class::String, n_range::AbstractArray, t_range::AbstractArray, n_jitter, t_jitter)
    if motif_class == "0"
        motif = rand_0(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "I"
        motif = rand_I(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "II"
        motif = rand_II(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "III"
        motif = rand_III(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "IV"
        motif = rand_IV(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "V"
        motif = rand_V(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "VI"
        motif = rand_VI(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "VII"
        motif = rand_VII(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "VIII"
        motif = rand_VIII(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "X"
        motif = rand_X(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "IX"
        motif = rand_IX(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        ) (motif, motif_class)
        motif
    elseif motif_class == "XI"
        motif = rand_XI(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "XII"
        motif = rand_XII(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    elseif motif_class == "XIII"
        motif = rand_XIII(n_range, t_range, n_jitter, t_jitter)
        @assert motif_class == offset_motif_numeral(
            lag_motif_sequence_class(motif...)
        )
        motif
    else
        error("Invalid motif class: $motif_class")
    end
end

function load_raster(fn)
    arr = CSV.read(fn, Array{Bool} ∘ Tables.matrix, delim=",", header=false)
    return OffsetArray(arr, 1:size(arr,1), 0:(size(arr,2)-1))
end

function generate_random_raster(raster_size::Tuple, spike_prob=0.1)
    putative_inputs = rand(Float64, raster_size)
    Array{Bool}(putative_inputs .<= spike_prob)
end

function pad_motif_example(motif, pad_n, pad_t)
    motif_n, motif_t = size(motif)
    padded = zeros(Bool, motif_n + 2pad_n, motif_t + 2pad_t)
    padded[pad_n+1:pad_n+motif_n, pad_t+1:pad_t+motif_t] .= motif
    return padded
end

function pad_top_motif_example(motif, pad_n, pad_t)
    motif_n, motif_t = size(motif)
    padded = zeros(Bool, motif_n + 2pad_n, motif_t + 2pad_t)
    padded[1:motif_n, pad_t+1:pad_t+motif_t] .= motif
    return padded
end

function repeat_padded_motif(motif_name, pad_n, pad_t, repeat_n, repeat_t)
    padded_motif = pad_motif_example(motif_examples[motif_name], pad_n, pad_t)
    repeat(padded_motif, outer=(repeat_n, repeat_t))
end

function repeat_padded_top_motif(motif_name, pad_n, pad_t, repeat_n, repeat_t)
    padded_motif = pad_top_motif_example(motif_examples[motif_name], pad_n, pad_t)
    repeat(padded_motif, outer=(repeat_n, repeat_t))
end