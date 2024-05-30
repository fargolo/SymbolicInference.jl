module SymbolicInference

using StatsBase
using RecurrenceAnalysis
using AnalyticComb
using Distributions
using LinearAlgebra

using Colors
using CairoMakie


include("rqa_interface.jl")
include("explore.jl")
include("tda.jl")


export rec_matrix_probs, rec_matrix_motifs,
        extract_recurrences, plot_motifs, 
        extract_recurrences_cross , plot_motifs_cross ,
        extract_recurrences_joint , plot_motifs_joint ,
        persistence_motifs
        

end
