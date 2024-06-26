"""
    rec_matrix_motifs(rec_matrix::RecurrenceMatrix;seqs="recurrences",window_range=collect(1:6),n_motifs=2)

Returns set of probabilities associated with consecutive runs in off-diagonals.

Argument `seqs` sets the type of consecutive sequences: either 'double' (recurrences and non-recurrences), 
'recurrences' or 'poincare' (non-recurrences).  The diagonals given by `window_range` argument are considered,
along with n_motifs for each diagonal.  See `AnalyticComb.weighted_bin_runs_prob` for definition of symbolic construction.  
    
"""
function rec_matrix_motifs(
    rec_matrix::Union{RecurrenceMatrix,CrossRecurrenceMatrix,JointRecurrenceMatrix};
    seqs="recurrences", window_range=collect(1:6), n_motifs=2)

    if seqs ∉ ["double","recurrences","poincare"]
        println("'seqs' must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 

    mat_len = size(rec_matrix,1)
    
    if mat_len < window_range[end]
        println("Largest range in 'window_range' must be smaller than matrix length")
        return(NaN)
    end 
    
    
    probs = []#Vector{Union{Missing,Float64}}[]
    motifs_inds_duration = []#Vector{Union{Missing,{Tuple{Int64, Int64}}}}[]

    # p and q values
    p = RecurrenceAnalysis.recurrencerate(rec_matrix)
    q = 1-p
    println("\n RR is: ",p)
    println("P and Q are: ",p," and ",q)

    for i in window_range
        # add Try catch to entire loop over window
        try
            cur_len = mat_len - i
            col_counts = StatsBase.rle(LinearAlgebra.diag(Matrix(rec_matrix),i))
            zipped_tups = collect(zip(col_counts[2],col_counts[1]))
            
            if seqs == "recurrences"
                zipped_tups = filter(x -> x[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 1st value is 1 
            elseif seqs == "poincare"
                zipped_tups = filter(x -> x[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 1st value is 1 
            end

            println("\n Zipped tuples: ")
            print(zipped_tups)
            println("\n Current segment length: ")
            print(cur_len)
            
            
            cur_counts = first.(zipped_tups)
            max_vals_inds = partialsortperm(cur_counts, 1:n_motifs, rev=true)
            max_vals = cur_counts[max_vals_inds]

            println("\n Diagonal largest sequences sizes and current total length are")
            print(max_vals," and ", cur_len)
            cur_probs = map(x -> AnalyticComb.weighted_bin_runs_pval(p,q,x,cur_len),max_vals)
            println("\n Probabilities:")
            print(cur_probs)
            push!(probs,cur_probs)


            # Sum counts (stored in col_counts[2]) of sequences appearing before current sequence col_counts[2][max_val_inds[x] - 1] 
            # to obtain position of current sequence  in original diagonal
            cur_inds = map(x -> sum(col_counts[2][1:(max_vals_inds[x]-1)]) + 1,1:n_motifs)             
            zipped_motifs_inds_duration = collect(zip(cur_inds, max_vals))
            push!(motifs_inds_duration, zipped_motifs_inds_duration)
                    
        catch e
            println("\n Check if sequence has zero occurences or n_motifs is larger than the number of motifs in diagonal")
            push!(probs,missing)
            push!(motifs_inds_duration,missing)  
        end
    end
    
    dict_keys = ["Window","Probs","Motifs starts and duration"]

    Dict(zip(dict_keys,
        [window_range, # Window
        probs,  # Probs
        motifs_inds_duration])) #Motif starts   

end