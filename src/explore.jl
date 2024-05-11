
"""
extract_recurrences(data_source::Vector{Float64}, motifs_dict::Dict{String, Vector}; num_windows::Int64 = 3)

This function returns x and y coordinates for a given window considering start and size of each motif.
The y coordinates are the values from data provided by the user. 

"""
function extract_recurrences(data_source::Vector{Float64}, 
    motifs_dict::Dict{String, Vector{Any}}; num_windows::Int64 = 3)
    
    probs_n = first(size(motifs_dict["Probs"][1]))
    probs_flat = vcat(motifs_dict["Probs"]...)
    sorted_probs_inds = sortperm(probs_flat)

    motifs_flat = vcat(motifs_dict["Motifs starts and duration"]...)
    motifs_by_window = motifs_flat[sorted_probs_inds[1:num_windows]]

    full_data_dict = []
    for (index, motif) in enumerate(motifs_by_window)
        motif_start, motif_size = motif
        motif_range = motif_start:(motif_start+motif_size)-1
        motif_y1 = data_source[motif_range]
        
        motif_rec_start = motif_start + div(sorted_probs_inds[index],probs_n)
        motif_rec_range = motif_rec_start:(motif_rec_start+motif_size)-1 
        motif_y2 = data_source[motif_rec_range]
        
        push!(full_data_dict, Dict(
        "x1" => collect(motif_range), "y1" => motif_y1,
        "x2" => collect(motif_rec_range), "y2" => motif_y2,
        "prob" => probs_flat[sorted_probs_inds][index]))
    end
    return full_data_dict
end

"""
    plot_motifs(time_series::Vector{Float64},coordinates::Vector{Any}; plot_size=(2000, 1000), n_motifs=2)

Plot motifs
"""
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
            rcolor = Colors.RGB(R1, G1, B1)
            lines!(p, motif["x1"], motif["y1"], linewidth=5, color=(rcolor,0.3))
            lines!(p, motif["x2"], motif["y2"], linewidth=5, color=(rcolor,0.3))
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
