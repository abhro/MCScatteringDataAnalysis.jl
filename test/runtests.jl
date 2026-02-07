using Test, MCScatteringDataAnalysis
using Aqua
using Random: randn, rand
using BiNormalDistributions: BiNormal
using Distributions: Normal

# Run test suite
println("Starting tests")
ti = time()

@testset "MCScatteringDataAnalysis" begin
    # check that all the functions at least run, i.e., there are no
    # syntax/import issues or anything
    @testset "API tests" begin

        dist, score = fit_dist_to_histogram(BiNormal, rand(5); params = rand(5))
        @test dist isa BiNormal
        @test score isa Float64
        dist, score = fit_dist_to_histogram(BiNormal, []; params = rand(5))
        @test ismissing(dist)
        @test ismissing(score)
        @test fit_dist_to_histogram(Normal, rand(5)) isa Normal
        @test fit_dist_to_histogram(Normal, []) |> ismissing
        @test fitnormal(rand(5)) isa Normal
        @test fitnormal([]) |> ismissing

        x, y = get_hist_curve(randn(3000); nbins = 150)
        @test x isa AbstractVector{Float64}
        @test y isa Vector{Float64}
    end

    @testset "Aqua.jl (Code quality)" begin
        Aqua.test_all(MCScatteringDataAnalysis)
    end
end

ti = time() - ti
println()
println("Test took total time of:")
println(round(ti / 60, digits = 3), " minutes")
