"""
    rec_matrix_probs(rec_matrix::RecurrenceMatrix;seqs="double",max_window=6)

Returns set of probabilities associated with consecutive runs in off-diagonals.

Argument `seqs` sets the type of consecutive sequences: either 'double' (recurrences and non-recurrences), 
'recurrences' or 'poincare' (non-recurrences).  The first n diagonals (given by `max_window` argument) are considered.
See also `AnalyticComb.p_val_weighted` and  `AnalyticComb.weighted_bin_runs_coeff`.  
"""
function rec_matrix_probs(rec_matrix::RecurrenceMatrix;seqs="recurrences",max_window=6)

    if seqs ∉ ["double","recurrences","poincare"]
        println("'seqs' must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 

    mat_len = size(rec_matrix,1)
    
    if mat_len < max_window
        println("'max_window' must be smaller than matrix length")
        return(NaN)
    end 
    
    
    probs = Float64[]

    # p and q values
    p = RecurrenceAnalysis.recurrencerate(rec_matrix)
    q = 1-p
    println("RR is: ",p)
    println("P and Q are: ",p," and ",q)

    for i in 1:max_window
        cur_len = mat_len - i
        col_counts = StatsBase.rle(LinearAlgebra.diag(Matrix(rec_matrix),i))
        zipped_tups = collect(zip(col_counts[2],col_counts[1]))
        
        if seqs == "recurrences"
            zipped_tups = filter(i -> i[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 1st value is 1 
        elseif seqs == "poincare"
            zipped_tups = filter(i -> i[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 1st value is 1 
        end

        println("\n Zipped tuples: ")
        print(zipped_tups)
        println("\n Current segment length: ")
        print(cur_len)
        
        try
            max_val = maximum(first.(zipped_tups))
            println("\n Diagonal largest sequence size and total length are")
            print(max_val," and ", cur_len)
            cur_p = AnalyticComb.weighted_bin_runs_prob(p,q,max_val,cur_len)
            println("\n Probability:")
            print(cur_p)
            push!(probs,cur_p)
        catch e
            println("\n Sequence has zero occurences.")
            push!(probs,NaN)
        end         
        
    end

    return(probs)

end



"""
    rec_matrix_motifs(rec_matrix::RecurrenceMatrix;seqs="recurrences",max_window=6,n_motifs=2)

Returns set of probabilities associated with consecutive runs in off-diagonals.

Argument `seqs` sets the type of consecutive sequences: either 'double' (recurrences and non-recurrences), 
'recurrences' or 'poincare' (non-recurrences).  The first n diagonals (given by `max_window` argument) are considered,
along with n_motifs for each diagonal.  See `AnalyticComb.weighted_bin_runs_prob` for definition of symbolic construction.  
    
"""
function rec_matrix_motifs(rec_matrix::RecurrenceMatrix;seqs="recurrences",max_window=6,n_motifs=2)

    if seqs ∉ ["double","recurrences","poincare"]
        println("'seqs' must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 

    mat_len = size(rec_matrix,1)
    
    if mat_len < max_window
        println("'max_window' must be smaller than matrix length")
        return(NaN)
    end 
    
    
    probs = Vector{Float64}[]
    motifs_inds_duration = Vector{Tuple{Int64, Int64}}[]

    # p and q values
    p = RecurrenceAnalysis.recurrencerate(rec_matrix)
    q = 1-p
    println("RR is: ",p)
    println("P and Q are: ",p," and ",q)

    for i in 1:max_window
        cur_len = mat_len - i
        col_counts = StatsBase.rle(LinearAlgebra.diag(Matrix(rec_matrix),i))
        zipped_tups = collect(zip(col_counts[2],col_counts[1]))
        
        if seqs == "recurrences"
            zipped_tups = filter(i -> i[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 1st value is 1 
        elseif seqs == "poincare"
            zipped_tups = filter(i -> i[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 1st value is 1 
        end

        println("\n Zipped tuples: ")
        print(zipped_tups)
        println("\n Current segment length: ")
        print(cur_len)
        
        try
            cur_counts = first.(zipped_tups)
            max_vals_inds = partialsortperm(cur_counts, 1:n_motifs, rev=true)
            max_vals = cur_counts[max_vals_inds]

            println("\n Diagonal largest sequences sizes and current total length are")
            print(max_vals," and ", cur_len)
            cur_probs = map(x -> AnalyticComb.weighted_bin_runs_prob(p,q,x,cur_len),max_vals)
            println("\n Probabilities:")
            print(cur_probs)
            push!(probs,cur_probs)


            # Sum counts (stored in col_counts[2]) of sequences appearing before current sequence col_counts[2][max_val_inds[x] - 1] 
            # to obtain position of current sequence  in original diagonal
        
            cur_inds = map(x -> sum(col_counts[2][1:(max_vals_inds[x]-1)]) + 1,1:n_motifs)             
            zipped_motifs_inds_duration = collect(zip(cur_inds, max_vals))
            push!(motifs_inds_duration, zipped_motifs_inds_duration)

        catch e
            println("\n Sequence has zero occurences or n_motifs is larger than the number of motifs in diagonal.")
            push!(probs,NaN)
            push!(motifs_inds_duration,NaN) 
        end         
        
    end
    
    dict_keys = ["Window","Probs","Motifs starts and duration"]

    Dict(zip(dict_keys,
        [collect(1:max_window), # Window
        probs,  # Probs
        motifs_inds_duration])) #Motif starts   

end

"""
    get_recurrence_points_from_data(data_source::Vector{Int64}, motifs::Vector{Tuple{Int64, Int64}}, window=2)

    This function returns x and y coordinates for a given window considering start and size of each motif.
    The y coordinates are the values from data provided by the user. 

"""
function get_recurrence_points_from_data(data_source::Vector{Float64}, 
    motifs_dict::Dict{String, Vector}; num_windows::Int64 = 3)
    
    sorted_probs_inds = sortperm(map(x -> minimum(x), motifs_dict["Probs"]))
    motifs_by_window = motifs_dict["Motifs starts and duration"][sorted_probs_inds[1:num_windows]]

    full_data_dict = []

    for (index, window) in enumerate(motifs_by_window)

        for motif in window
            first_item, first_size = motif
            first_range = first_item: first_item+first_size-1
            first_current_data_y = data_source[first_range]
            

            second_item = first_item + sorted_probs_inds[index]
            second_size = first_size

            second_range = second_item: second_item+second_size -1 
            second_current_data_y = data_source[second_range]

            push!(full_data_dict, Dict("x1" => collect(first_range),  "y1" => first_current_data_y, 
                                       "x2" => collect(second_range), "y2" => second_current_data_y))
        end
    end
    return full_data_dict
end

function plot_motifs(time_series::Vector{Float64}, coordinates::Vector{Any}; 
    plot_size=(2000, 1000), n_motifs=2)
    p = plot(1:length(time_series), time_series, alpha=0.1, 
    label="Time Series", size=plot_size, linewidth=5, color=:black, legend=false)
    if n_motifs == 1
        first_coordinate = coordinates[1]
        plot!(p, first_coordinate["x1"], first_coordinate["y1"], linewidth=5, color=:yellow, alpha=0.3)
        plot!(p, first_coordinate["x2"], first_coordinate["y2"], linewidth=5, color=:purple, alpha=0.4)
    else
        for motif in coordinates[1:n_motifs]
            R1, G1, B1 = rand(3)
            color = Colors.RGB(R1, G1, B1)
            plot!(p, motif["x1"], motif["y1"], linewidth=5, color=color, alpha=0.6)
            plot!(p, motif["x2"], motif["y2"], linewidth=5, color=color, alpha=0.8)
            scatter!(p, motif["x1"], motif["y1"], color=color, alpha=0.6)
            scatter!(p, motif["x2"], motif["y2"], color=color, alpha=0.8)
        end
    end
    display(p)
end