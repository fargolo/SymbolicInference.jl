"""
double_inference_weighted(rec_matrix::RecurrenceMatrix;seqs="double",max_window=6)

Returns set of probabilities associated with consecutive runs in off-diagonals.

Argument `seqs` sets the type of consecutive sequences.  
Could be either 'double' (recurrences and non-recurrences), 'recurrences' or 'poincare' (non-recurrences).  
The first n diagonals (given by `max_window` argument) are considered.
See `AnalyticComb.p_val_weighted` and  `AnalyticComb.weighted_bin_runs_coeff`.  
"""
function double_inference_weighted(rec_matrix::RecurrenceMatrix;seqs="double",max_window=6)

    if seqs ∉ ["double","recurrences","poincare"]
        println("seqs must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 

    mat_len = dim(rec_matrix)
    
    if mat_len < max_window
        println("max_window must be smaller than matrix length")
        return(NaN)
    end 
    
    
    probs = Float64[]

    p = RecurrenceAnalysis.recurrencerate(rec_matrix)
    q = 1-p
    println("RR is ",p)
    println("P and Q are ",p," and ",q)

    

    # To do: implement map instead of for loop 
    # sequences = map(x-> diag(Matrix(rec_matrix),x), 1:size(rec_matrix)[1])

    for i in 1:max_window
        cur_len = mat_len - i
        col_counts = StatsBase.rle(LinearAlgebra.diag(Matrix(rec_matrix),i))
        zipped_tups = collect(zip(col_counts[2],col_counts[1]))
        
        if seqs == "recurrences"
            zipped_tups = filter(i -> i[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 1st value is 1 
        elseif seqs == "poincare"
            zipped_tups = filter(i -> i[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 1st value is 1 
        end

        println("Zipped tuples: ")
        print(zipped_tups)
        println("\n Current segment length: ")
        print(cur_len)
        
        try
            max_val = maximum(first.(zipped_tups))
            println("\nDiagonal largest sequence size and total length are")
            print(max_val,cur_len)
            cur_p = AnalyticComb.p_val_weighted(p,q,max_val,cur_len)
            println("\nProbability:")
            print(cur_p)
            push!(probs,cur_p)
        catch e
            println("\nSequence has zero occurences.")
            push!(probs,NaN)
        end         
        
    end

    return(probs)

end


