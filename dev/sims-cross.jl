#using Pkg
#Pkg.add(path="/home/dellint5/repos/SymbolicInference.jl")

using SymbolicInference
using Random , Distributions
using RecurrenceAnalysis
using CairoMakie, GLMakie
using Colors

cycles = 10

rec_rate = 0.3
sd_norm = 0.5
p_window_size = 12
linear_basis = (0:0.6:2*pi*cycles)
n_sample = length(linear_basis)

time_series_sin = map(sin,linear_basis)

# White noise
rng = Random.MersenneTwister(42)
dist = Distributions.Normal(0,sd_norm)
white_n_time_series = rand(dist,n_sample)

# Noisy sine wave
time_series_sin_nois = white_n_time_series .+ time_series_sin


res_recs = map(x -> RecurrenceAnalysis.CrossRecurrenceMatrix(x,time_series_sin,rec_rate;fixedrate=true),
    [time_series_sin,white_n_time_series,time_series_sin_nois])
res_recs_joint = map(x -> RecurrenceAnalysis.JointRecurrenceMatrix(x,time_series_sin,rec_rate;fixedrate=true),
    [time_series_sin,white_n_time_series,time_series_sin_nois])
    #rec_matrix_motifs(rec_matrix::RecurrenceMatrix;seqs="recurrences",max_window=6,n_motifs=2)
all_probs = map(x-> rec_matrix_motifs(x;
        max_window=p_window_size,seqs="recurrences"),res_recs)
all_probs_joint = map(x-> rec_matrix_motifs(x;
        max_window=p_window_size,seqs="recurrences"),res_recs_joint)

motifs_dict = all_probs[3]
motifs_dict_joint = all_probs_joint[3]

coordinates = SymbolicInference.extract_recurrences_cross(time_series_sin_nois,time_series_sin,
            motifs_dict;num_windows = 12)

coordinates_joint = SymbolicInference.extract_recurrences_cross(time_series_sin_nois,time_series_sin,
    motifs_dict_joint;num_windows = 12)

GLMakie.activate!(inline=false)
p = plot_motifs_cross(time_series_sin_nois,time_series_sin,
    coordinates;n_motifs=2)
p2 = plot_motifs_cross(time_series_sin_nois,time_series_sin,
    coordinates_joint;n_motifs=3)


