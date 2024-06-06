function persistence_motifs(time_series; range = collect(0.1:0.1:0.9), n_windows=10)
    all_motifs = []

    for cur_range in range
        time_series_mat = RecurrenceAnalysis.RecurrenceMatrix(time_series, cur_range; fixedrate=true)
        motifs_result = SymbolicInference.rec_matrix_motifs(time_series_mat; max_window = n_windows, n_motifs=1)
        push!(all_motifs, motifs_result)
    end

    num_windows = length(all_motifs[1]["Probs"]) # Assuming all have the same number of motifs

    fig = Figure()
    ax = Axis(fig[1, 1], limits = (0, 1, 0, 1),
    xlabel="Recurrence Rates", 
        ylabel="Probabilities")

    for j in 1:num_windows
        x_vals = range
        y_vals = [all_motifs[i]["Probs"][j] for i in 1:length(range)]
        lines!(ax, x_vals, vcat(y_vals...), label="Window $j")
    end

    Legend(fig[1, 2], ax, "Motifs")
    fig
end

function persistence_barcode(time_series; range = collect(0.1:0.1:0.9), n_windows=10,alpha_thresh=0.05)
    all_motifs = []

    for cur_range in range
        time_series_mat = RecurrenceMatrix(time_series, cur_range; fixedrate=true)
        motifs_result = rec_matrix_motifs(time_series_mat; max_window=n_windows, n_motifs=1)
        push!(all_motifs, motifs_result)
    end

    fig = Figure()
    ax = Axis(fig[1, 1], limits=(0, 1, 0, n_windows + 1),
              xlabel="Recurrence Rates", ylabel="Windows")

    for cur_window in 1:n_windows
        x_vals = range
        y_vals = [all_motifs[j]["Probs"][cur_window] for j in 1:length(range)]
        inds = findall(x -> x < alpha_thresh, vcat(y_vals...))
        scatter!(ax, x_vals[inds], fill(cur_window,length(inds)) , label="Window $cur_window")

    end
    fig
end

# Example usage
#time_series = rand(100) # Replace with your actual time series data
#fig = persistence_barcode(time_series)
#display(fig)