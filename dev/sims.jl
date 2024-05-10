#using Pkg
#Pkg.add(path="/home/dellint5/repos/SymbolicInference.jl")

using SymbolicInference
using Random , Distributions
using RecurrenceAnalysis
using CairoMakie, GLMakie
using Colors

cycles = 70

rec_rate = 0.3
sd_norm = 0.2
p_window_size = 20
linear_basis = (0:0.6:2*pi*cycles)
n_sample = length(linear_basis)

time_series_sin = map(sin,linear_basis)

# White noise
rng = Random.MersenneTwister(42)
dist = Distributions.Normal(0,sd_norm)
white_n_time_series = rand(dist,n_sample)

# Noisy sine wave
time_series_sin_nois = white_n_time_series .+ time_series_sin


res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,rec_rate;fixedrate=true),
    [time_series_sin,white_n_time_series,time_series_sin_nois])
#rec_matrix_motifs(rec_matrix::RecurrenceMatrix;seqs="recurrences",max_window=6,n_motifs=2)
all_probs = map(x-> SymbolicInference.rec_matrix_motifsv3(x;
        max_window=p_window_size,seqs="recurrences"),res_recs)


motifs_dict = all_probs[3]

coordinates = SymbolicInference.extract_recurrences(time_series_sin_nois,motifs_dict;num_windows = 12)
GLMakie.activate!(inline=false)
p = plot_motifs(time_series_sin_nois,coordinates;n_motifs=8)
p

function plot_motifs(time_series::Vector{Float64}, 
    coordinates::Vector{Any}; plot_size=(2000, 1000), n_motifs=2)

    pl = CairoMakie.Figure(size = plot_size)
    p = Axis(pl[1, 1])
    lines!(p,1:length(time_series), time_series, color=:black)
    if n_motifs == 1
        first_coordinate = coordinates[1]
        lines!(p, first_coordinate["x1"], first_coordinate["y1"], linewidth=5, color=:yellow)
        lines!(p, first_coordinate["x2"], first_coordinate["y2"], linewidth=5, color=:purple)
    else
        for motif in coordinates[1:n_motifs]
            R1, G1, B1 = rand(3)
            rcolor = (Colors.RGB(R1, G1, B1),0.3)
            lines!(p, motif["x1"], motif["y1"], linewidth=5, color=rcolor)
            lines!(p, motif["x2"], motif["y2"], linewidth=5, color=rcolor)
            CairoMakie.scatter!(p, motif["x1"], motif["y1"], color=rcolor)
            CairoMakie.scatter!(p, motif["x2"], motif["y2"], color=rcolor)
            CairoMakie.bracket!(p,
                text = "p-value: "*string(round(motif["prob"],digits=3)),
                motif["x1"][1], motif["y1"][1],
                motif["x2"][1], motif["y2"][1], color=rcolor)
        end
    end
    return pl
end
