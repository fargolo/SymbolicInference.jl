module SymbolicInference

using StatsBase
using RecurrenceAnalysis
using AnalyticComb
using Distributions
using LinearAlgebra

using Colors
using CairoMakie, GLMakie


include("recurrence_plots.jl")
include("explore.jl")


export rec_matrix_probs, rec_matrix_motifs,
        extract_recurrences, plot_motifs 

end
