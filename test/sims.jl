using Plots

#using Pkg
#Pkg.add(path="/home/dellint5/repos/SymbolicInference.jl")

using SymbolicInference
using Random , Distributions
using RecurrenceAnalysis
using CairoMakie, GLMakie

cycles = 12

rec_rate = 0.2
sd_norm = 0.2
p_window_size = 15
linear_basis = (0:0.4:2*pi*cycles)
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

coordinates = extract_recurrences(time_series_sin_nois,motifs_dict;num_windows = 12)
GLMakie.activate!(inline=false)
p = plot_motifs(time_series_sin_nois,coordinates;n_motifs=12)
p

