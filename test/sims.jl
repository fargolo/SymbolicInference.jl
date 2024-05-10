using Plots

#using Pkg
#Pkg.add(path="/home/dellint5/repos/SymbolicInference.jl")

using SymbolicInference
using Random , Distributions
using RecurrenceAnalysis
using CairoMakie, GLMakie

cycles = 8

rec_rate = 0.2
sd_norm = 0.5
p_window_size = 20
linear_basis = (0:0.5:2*pi*cycles)
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
all_probs = map(x-> rec_matrix_motifs(x;
        max_window=p_window_size,seqs="recurrences"),res_recs)


motifs_dict = all_probs[3]

coordinates = extract_recurrences(time_series_sin_nois,motifs_dict)
GLMakie.activate!(inline=true)
p = plot_motifs(time_series_sin_nois,coordinates)
p
CairoMakie.display(p)

function plot_motifs(time_series::Vector{Float64}, 
    coordinates::Vector{Any}; plot_size=(2000, 1000), n_motifs=2)

    p = Axis(Figure(size = plot_size)[1, 1])
    lines!(p,1:length(time_series), time_series, color=:black)
    if n_motifs == 1
        first_coordinate = coordinates[1]
        lines!(p, first_coordinate["x1"], first_coordinate["y1"], linewidth=5, color=:yellow)
        lines!(p, first_coordinate["x2"], first_coordinate["y2"], linewidth=5, color=:purple)
    else
        for motif in coordinates[1:n_motifs]
            R1, G1, B1 = rand(3)
            color = Colors.RGB(R1, G1, B1)
            lines!(p, motif["x1"], motif["y1"], linewidth=5, color=color)
            lines!(p, motif["x2"], motif["y2"], linewidth=5, color=color)
            CairoMakie.scatter!(p, motif["x1"], motif["y1"], color=color)
            CairoMakie.scatter!(p, motif["x2"], motif["y2"], color=color)
            CairoMakie.bracket!(p,motif["x1"], motif["y1"],motif["x2"], motif["y2"])
        end
    end
    display(p)
end
