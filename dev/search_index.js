var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SymbolicInference","category":"page"},{"location":"#SymbolicInference","page":"Home","title":"SymbolicInference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SymbolicInference.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SymbolicInference]","category":"page"},{"location":"#SymbolicInference.extract_recurrences-Tuple{Vector{Float64}, Dict{String, Vector{Any}}}","page":"Home","title":"SymbolicInference.extract_recurrences","text":"extractrecurrences(datasource::Vector{Float64}, motifsdict::Dict{String, Vector}; numwindows::Int64 = 3)\n\nThis function returns x and y coordinates for a given window considering start and size of each motif. The y coordinates are the values from data provided by the user. \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.extract_recurrences_cross-Tuple{Vector{Float64}, Vector{Float64}, Dict{String, Vector{Any}}}","page":"Home","title":"SymbolicInference.extract_recurrences_cross","text":"extract_recurrences_cross(data_source::Vector{Float64}, data_source2::Vector{Float64},\n    motifs_dict::Dict{String, Vector}; num_windows::Int64 = 3)\n\nThis function returns x and y coordinates for a given window      considering start and size of each motif detected from      two time-series in cross-recurrence matrices . The y coordinates are the values from data provided by the user. \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.extract_recurrences_joint-Tuple{Vector{Float64}, Vector{Float64}, Dict{String, Vector{Any}}}","page":"Home","title":"SymbolicInference.extract_recurrences_joint","text":"extract_recurrences_joint(data_source::Vector{Float64}, data_source2::Vector{Float64},\n    motifs_dict::Dict{String, Vector}; num_windows::Int64 = 3)\n\n    This function returns x and y coordinates for a given window \n        considering start and size of each motif detected from \n        two time-series in joint-recurrence matrices .\n\nThe y coordinates are the values from data provided by the user. \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.persistence_barcode-Tuple{Any, Any}","page":"Home","title":"SymbolicInference.persistence_barcode","text":"persistencebarcode(timeseries1,timeseries2; range = collect(0.1:0.1:0.9), nwindows=10,alpha_thresh=0.05)\n\nReturn barcode plot for each window.      Points are plotted whenever the p-value is smaller than alpha_thresh.  \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.persistence_barcode-Tuple{Any}","page":"Home","title":"SymbolicInference.persistence_barcode","text":"persistencebarcode(timeseries; range = collect(0.1:0.1:0.9), nwindows=10,alphathresh=0.05)\n\nReturn barcode plot for each window.      Points are plotted whenever the p-value is smaller than alpha_thresh.  \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.plot_motifs-Tuple{Vector{Float64}, Vector{Any}}","page":"Home","title":"SymbolicInference.plot_motifs","text":"plot_motifs(time_series::Vector{Float64},coordinates::Vector{Any}; plot_size=(2000, 1000), n_motifs=2)\n\nPlot motifs from coordinates extracted with extract_recurrences.  \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.plot_motifs_cross-Tuple{Vector{Float64}, Vector{Float64}, Vector{Any}}","page":"Home","title":"SymbolicInference.plot_motifs_cross","text":"plot_motifs_cross(time_series::Vector{Float64},time_series2::Vector{Float64},\n    coordinates::Vector{Any}; plot_size=(2000, 1000), n_motifs=2)\n\nPlot motifs from coordinates extracted with extract_recurrences_cross.  \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.plot_motifs_joint-Tuple{Vector{Float64}, Vector{Float64}, Vector{Any}}","page":"Home","title":"SymbolicInference.plot_motifs_joint","text":"plot_motifs_joint(time_series::Vector{Float64},time_series2::Vector{Float64},\n    coordinates::Vector{Any}; plot_size=(2000, 1000), n_motifs=2)\n\nPlot motifs from coordinates extracted with extract_recurrences_joint.  \n\n\n\n\n\n","category":"method"},{"location":"#SymbolicInference.rec_matrix_motifs-Tuple{Union{RecurrenceAnalysis.CrossRecurrenceMatrix, RecurrenceAnalysis.JointRecurrenceMatrix, RecurrenceAnalysis.RecurrenceMatrix}}","page":"Home","title":"SymbolicInference.rec_matrix_motifs","text":"rec_matrix_motifs(rec_matrix::RecurrenceMatrix;seqs=\"recurrences\",window_range=collect(1:6),n_motifs=2)\n\nReturns set of probabilities associated with consecutive runs in off-diagonals.\n\nArgument seqs sets the type of consecutive sequences: either 'double' (recurrences and non-recurrences),  'recurrences' or 'poincare' (non-recurrences).  The diagonals given by window_range argument are considered, along with nmotifs for each diagonal.  See `AnalyticComb.weightedbinrunsprob` for definition of symbolic construction.  \n\n\n\n\n\n","category":"method"}]
}
