using SymbolicInference
using RecurrenceAnalysis
using Random
using Test

@testset "SymbolicInference.jl" begin
    rng = Random.MersenneTwister(42)
    t = collect(0:0.4:(2*pi*15))
    sin_ts = sin.(t)
    noise_ts = randn(rng, length(t))

    rm_sin = RecurrenceMatrix(sin_ts, 0.3; fixedrate=true)
    rm_noise = RecurrenceMatrix(noise_ts, 0.3; fixedrate=true)

    # Test rec_matrix_motifs
    motifs_sin = rec_matrix_motifs(rm_sin; window_range=collect(1:6), n_motifs=2)
    @test haskey(motifs_sin, "Window")
    @test haskey(motifs_sin, "Probs")
    @test haskey(motifs_sin, "Motifs starts and duration")
    @test length(motifs_sin["Probs"]) == 6

    # Test rqa_moments_test: Determinism (length_ge), Motif count (count_ge), Exact length (count_eq)
    det_sin = rqa_moments_test(rm_sin; lmin=3, window_range=collect(1:10), metric="determinism")
    det_noise = rqa_moments_test(rm_noise; lmin=3, window_range=collect(1:10), metric="determinism")

    # Periodic sine wave must exhibit significantly higher determinism than analytical null and white noise
    @test det_sin["Global"].z_score > 3.0
    @test det_sin["Global"].p_value < 0.005
    @test det_sin["Global"].det_observed > det_sin["Global"].det_expected
    @test det_sin["Global"].z_score > det_noise["Global"].z_score

    # Test motif_count and exact_length modes
    cnt_sin = rqa_moments_test(rm_sin; lmin=3, window_range=collect(1:6), metric="motif_count")
    @test length(cnt_sin["Observed"]) == 6
    @test isfinite(cnt_sin["Global"].z_score)

    eq_noise = rqa_moments_test(rm_noise; lmin=2, window_range=collect(1:6), metric="exact_length")
    @test length(eq_noise["Expected"]) == 6
    @test all(>(0.0), eq_noise["Std"])

    # Test poincare and double sequence modes
    poinc_sin = rqa_moments_test(rm_sin; lmin=3, window_range=collect(1:6), seqs="poincare", metric="determinism")
    @test poinc_sin["Global"].z_score > 3.0

    dbl_sin = rqa_moments_test(rm_sin; lmin=3, window_range=collect(1:6), seqs="double", metric="motif_count")
    @test isfinite(dbl_sin["Global"].expected)
end
