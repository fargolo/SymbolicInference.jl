"""
    rec_matrix_motifs(rec_matrix::Union{RecurrenceMatrix,CrossRecurrenceMatrix,JointRecurrenceMatrix};
                      seqs="recurrences", window_range=collect(1:6), n_motifs=2, verbose=false)

Returns set of probabilities associated with consecutive runs in off-diagonals.

Argument `seqs` sets the type of consecutive sequences: either 'double' (recurrences and non-recurrences), 
'recurrences' or 'poincare' (non-recurrences).  The diagonals given by `window_range` argument are considered,
along with n_motifs for each diagonal.  See `AnalyticComb.weighted_bin_runs_prob` for definition of symbolic construction.  
    
"""
function rec_matrix_motifs(
    rec_matrix::Union{RecurrenceMatrix,CrossRecurrenceMatrix,JointRecurrenceMatrix};
    seqs="recurrences", window_range=collect(1:6), n_motifs=2, verbose::Bool=false)

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
    if verbose
        println("\n RR is: ",p)
        println("P and Q are: ",p," and ",q)
    end

    for i in window_range
        # add Try catch to entire loop over window
        try
            cur_len = mat_len - i
            col_counts = StatsBase.rle(LinearAlgebra.diag(Matrix(rec_matrix),i))
            zipped_tups = collect(zip(col_counts[2],col_counts[1]))
            
            if seqs == "recurrences"
                zipped_tups = filter(x -> x[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 2nd value is 1 
            elseif seqs == "poincare"
                zipped_tups = filter(x -> x[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 2nd value is 0 
            end

            if verbose
                println("\n Zipped tuples: ")
                print(zipped_tups)
                println("\n Current segment length: ")
                print(cur_len)
            end
            
            cur_counts = first.(zipped_tups)
            max_vals_inds = partialsortperm(cur_counts, 1:n_motifs, rev=true)
            max_vals = cur_counts[max_vals_inds]

            if verbose
                println("\n Diagonal largest sequences sizes and current total length are")
                print(max_vals," and ", cur_len)
            end
            p_eff, q_eff = seqs == "poincare" ? (q, p) : (p, q)
            cur_probs = map(x -> AnalyticComb.weighted_bin_runs_pval(p_eff, q_eff, x, cur_len), max_vals)
            if verbose
                println("\n Probabilities:")
                print(cur_probs)
            end
            push!(probs,cur_probs)


            # Sum counts (stored in col_counts[2]) of sequences appearing before current sequence col_counts[2][max_val_inds[x] - 1] 
            # to obtain position of current sequence  in original diagonal
            cur_inds = map(x -> sum(col_counts[2][1:(max_vals_inds[x]-1)]) + 1,1:n_motifs)             
            zipped_motifs_inds_duration = collect(zip(cur_inds, max_vals))
            push!(motifs_inds_duration, zipped_motifs_inds_duration)
                    
        catch e
            if verbose
                println("\n Check if sequence has zero occurences or n_motifs is larger than the number of motifs in diagonal")
            end
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


"""
    rqa_moments_test(rec_matrix::Union{RecurrenceMatrix,CrossRecurrenceMatrix,JointRecurrenceMatrix};
                     lmin::Int=2, window_range=collect(1:6), seqs="recurrences", metric="determinism")

Perform analytical hypothesis testing on Recurrence Quantification Analysis (RQA) diagonal statistics
using exact moments (`mean`, `var`, `std`) derived from Bivariate Generating Functions (BGFs) via `AnalyticComb.weighted_bin_runs_moments`.

# Arguments
- `rec_matrix`: A `RecurrenceMatrix`, `CrossRecurrenceMatrix`, or `JointRecurrenceMatrix`.
- `lmin::Int=2`: Minimum run length threshold ``l_{\\min}`` (standard RQA default is 2).
- `window_range=collect(1:6)`: Off-diagonal window indices to evaluate.
- `seqs="recurrences"`: Sequence type (`"recurrences"` for 1-runs, `"poincare"` for 0-runs, or `"double"` for both).
- `metric="determinism"`:
  - `"determinism"` (or `"length_ge"`): Tests the total number of points in runs of length ``\\ge l_{\\min}`` (the numerator of RQA Determinism ``DET``).
  - `"motif_count"` (or `"count_ge"`): Tests the number of distinct diagonal lines (runs) of length ``\\ge l_{\\min}``.
  - `"exact_length"` (or `"count_eq"`): Tests the count of diagonal lines of exact length ``l_{\\min}`` (line-length histogram ``P(l_{\\min})``).

# Returns
A `Dict{String, Any}` containing:
- `"Window"`: Evaluated `window_range`
- `"Observed"`: Vector of observed statistics per diagonal
- `"Expected"`: Vector of exact analytical means ``\\mu_i`` per diagonal
- `"Std"`: Vector of exact analytical standard deviations ``\\sigma_i`` per diagonal
- `"ZScores"`: Vector of standardized ``Z``-scores per diagonal
- `"PValues"`: Vector of one-sided upper-tail normal ``p``-values per diagonal
- `"Global"`: A `NamedTuple` `(observed, expected, var, std, z_score, p_value, p_value_twosided, chebyshev_bound, det_observed, det_expected)`
  aggregating across all diagonals in `window_range`.
"""
function rqa_moments_test(
    rec_matrix::Union{RecurrenceMatrix,CrossRecurrenceMatrix,JointRecurrenceMatrix};
    lmin::Int=2, window_range=collect(1:6), seqs="recurrences", metric="determinism")

    if seqs ∉ ["double", "recurrences", "poincare"]
        error("'seqs' must be either 'double', 'recurrences' or 'poincare'")
    end

    mat_len = size(rec_matrix, 1)
    if mat_len - maximum(window_range) < lmin
        error("Diagonal length (mat_len - maximum(window_range)) must be at least lmin ($lmin)")
    end

    bgf_mode = if metric in ("determinism", "length_ge")
        seqs == "double" ? error("metric='determinism' requires seqs='recurrences' or 'poincare'") : "length_ge"
    elseif metric in ("motif_count", "count_ge")
        seqs == "double" ? "double_count_ge" : "count_ge"
    elseif metric in ("exact_length", "count_eq")
        seqs == "double" ? error("metric='exact_length' requires seqs='recurrences' or 'poincare'") : "count_eq"
    else
        error("Unsupported metric '$metric'. Choose 'determinism', 'motif_count', or 'exact_length'.")
    end

    p = Float64(RecurrenceAnalysis.recurrencerate(rec_matrix))
    q = 1.0 - p
    p_eff, q_eff = seqs == "poincare" ? (q, p) : (p, q)

    mat = Matrix(rec_matrix)
    obs_vec = Float64[]
    exp_vec = Float64[]
    var_vec = Float64[]
    std_vec = Float64[]
    z_vec = Float64[]
    pval_vec = Float64[]

    total_active_points = 0
    total_expected_active = 0.0

    for i in window_range
        cur_len = mat_len - i
        diag_vec = LinearAlgebra.diag(mat, i)
        vals, lens = StatsBase.rle(diag_vec)

        if seqs == "recurrences"
            run_lens = [lens[k] for k in eachindex(vals) if vals[k] == 1]
            total_active_points += count(==(1), diag_vec)
        elseif seqs == "poincare"
            run_lens = [lens[k] for k in eachindex(vals) if vals[k] == 0]
            total_active_points += count(==(0), diag_vec)
        else # "double"
            run_lens = lens
            total_active_points += cur_len
        end
        total_expected_active += cur_len * (seqs == "double" ? 1.0 : p_eff)

        obs_val = if bgf_mode == "length_ge"
            Float64(sum(r for r in run_lens if r >= lmin; init=0))
        elseif bgf_mode in ("count_ge", "double_count_ge")
            Float64(count(r -> r >= lmin, run_lens))
        else # "count_eq"
            Float64(count(r -> r == lmin, run_lens))
        end

        mom = AnalyticComb.weighted_bin_runs_moments(p_eff, q_eff, lmin, cur_len; mode=bgf_mode)
        z_i = mom.std > 0 ? (obs_val - mom.mean) / mom.std : 0.0
        p_i = Distributions.ccdf(Distributions.Normal(0.0, 1.0), z_i)

        push!(obs_vec, obs_val)
        push!(exp_vec, mom.mean)
        push!(var_vec, mom.var)
        push!(std_vec, mom.std)
        push!(z_vec, z_i)
        push!(pval_vec, p_i)
    end

    tot_obs = sum(obs_vec)
    tot_exp = sum(exp_vec)
    tot_var = sum(var_vec)
    tot_std = sqrt(max(0.0, tot_var))
    tot_z = tot_std > 0 ? (tot_obs - tot_exp) / tot_std : 0.0
    tot_pval = Distributions.ccdf(Distributions.Normal(0.0, 1.0), tot_z)
    tot_pval_2s = 2.0 * Distributions.ccdf(Distributions.Normal(0.0, 1.0), abs(tot_z))
    cheb_bound = tot_z > 0 ? 1.0 / (1.0 + tot_z^2) : 1.0

    det_obs = total_active_points > 0 ? tot_obs / total_active_points : NaN
    det_exp = total_expected_active > 0 ? tot_exp / total_expected_active : NaN

    global_summary = (
        observed = tot_obs,
        expected = tot_exp,
        var = tot_var,
        std = tot_std,
        z_score = tot_z,
        p_value = tot_pval,
        p_value_twosided = tot_pval_2s,
        chebyshev_bound = cheb_bound,
        det_observed = bgf_mode == "length_ge" ? det_obs : NaN,
        det_expected = bgf_mode == "length_ge" ? det_exp : NaN
    )

    return Dict{String, Any}(
        "Window" => collect(window_range),
        "Observed" => obs_vec,
        "Expected" => exp_vec,
        "Std" => std_vec,
        "ZScores" => z_vec,
        "PValues" => pval_vec,
        "Global" => global_summary
    )
end
