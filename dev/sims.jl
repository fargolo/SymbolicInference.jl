#using Pkg
#Pkg.add(path="/home/dellint5/repos/SymbolicInference.jl")

using SymbolicInference
using Random , Distributions
using RecurrenceAnalysis
using CairoMakie
using Colors

cycles = 30

rec_rate = 0.6
sd_norm = 0.5
p_window_size = 30
linear_basis = (0:0.6:2*pi*cycles)
n_sample = length(linear_basis)

time_series_sin = map(sin,linear_basis)

# White noise
rng = Random.MersenneTwister(42)
dist = Distributions.Normal(0,sd_norm)
white_n_time_series = rand(dist,n_sample)

# Noisy sine wave
time_series_sin_nois = white_n_time_series .+ time_series_sin

pl_persist_p = persistence_motifs(time_series_sin;n_windows=30)
pl_persist_bar = persistence_barcode(time_series_sin_nois;n_windows=30)

res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,rec_rate;fixedrate=true),
    [time_series_sin,white_n_time_series,time_series_sin_nois])
#rec_matrix_motifs(rec_matrix::RecurrenceMatrix;seqs="recurrences",max_window=6,n_motifs=2)
all_probs = map(x-> SymbolicInference.rec_matrix_motifs(x;
        max_window=p_window_size,seqs="recurrences",n_motifs=2),res_recs)


motifs_dict = all_probs[3]

coordinates = SymbolicInference.extract_recurrences(time_series_sin,motifs_dict;num_windows = 10)
GLMakie.activate!(inline=false)
p = plot_motifs(time_series_sin_nois,coordinates;n_motifs=3)
p
