motif_examples = Dict(
    "I" => [0 0 0;
            0 1 0;
            0 0 0],
    "II" => [0 0 0;
             1 1 0;
             0 0 0],
    "III" => [0 0 0;
              1 1 1;
              0 0 0],
    "IV" => [0 1 0;
             0 1 0;
             0 0 0],
    "V" => [0 1 0;
            0 1 0;
            0 1 0],
    "VI" => [1 0 0;
             0 1 0;
             0 0 0],
    "VII" => [1 0 0;
              1 1 0;
              0 0 0],
    "VIII" => [1 1 0;
             0 1 0;
             0 0 0],
    "IX" => [0 1 1;
             1 0 0;
             0 0 0],
    "X" => [1 1 0;
            0 0 1;
            0 0 0],
    "XI" => [1 0 1;
             0 1 0;
             0 0 0],
    "XII" => [1 0 0;
              0 1 0;
              0 1 0],
    "XIII" => [1 0 0;
               0 1 0;
               1 0 0],
    "XIV" => [1 0 0;
               0 1 0;
               0 0 1]
)


function load_raster(fn)
    arr = CSV.read(fn, BitArray ∘ Tables.matrix, delim=",", header=false)
    return OffsetArray(arr, 1:size(arr,1), 0:(size(arr,2)-1))
end

function generate_random_raster(raster_size::Tuple, spike_prob=0.1)
    putative_inputs = rand(Float64, raster_size)
    BitArray(putative_inputs .<= spike_prob)
end

function pad_motif_example(motif, pad_n, pad_t)
    motif_n, motif_t = size(motif)
    padded = zeros(Int, motif_n + 2pad_n, motif_t + 2pad_t)
    padded[pad_n+1:pad_n+motif_n, pad_t+1:pad_t+motif_t] .= motif
    return padded
end

function pad_top_motif_example(motif, pad_n, pad_t)
    motif_n, motif_t = size(motif)
    padded = zeros(Int, motif_n + 2pad_n, motif_t + 2pad_t)
    padded[1+pad_n:motif_n+pad_n, pad_t+1:pad_t+motif_t] .= motif
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