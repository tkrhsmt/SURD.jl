using SURD
using Test

@testset "SURD.jl" begin
    @testset "Histogram" begin
        x1 = repeat([0, 0, 1, 1], 4)
        x2 = repeat([0, 1, 0, 1], 4)
        agents = (x1, x2)

        # Balanced factorial data make the histogram answers exact.
        unique = surd(x1, agents, SURD.Histogram())
        @test unique["u"] == Dict("1" => 1.0, "2" => 0.0)
        @test unique["r"]["1,2"] == unique["s"]["1,2"] == 0.0
        @test unique["leak"] ≈ 0.0 atol=1e-12

        xor_target = mod.(x1 .+ x2, 2)
        xor_result = surd(xor_target, agents, SURD.Histogram(); log=true)
        @test all(iszero, values(xor_result["u"]))
        @test xor_result["r"]["1,2"] == 0.0
        @test xor_result["s"]["1,2"] ≈ 1.0
        @test xor_result["leak"] ≈ 0.0 atol=1e-12

        independent = surd(vcat(zeros(16), ones(16)), (repeat(x1, 2), repeat(x2, 2)), SURD.Histogram())
        @test all(iszero, values(independent["u"]))
        @test independent["r"]["1,2"] == independent["s"]["1,2"] == 0.0
        @test independent["leak"] ≈ 1.0
        @test ismissing(surd(zeros(16), agents, SURD.Histogram())["leak"])
        @test_throws MethodError surd(x1, agents)

        @testset "Explicit bin count" begin
            one_bin = surd(x1, agents, SURD.Histogram(1))
            @test all(iszero, values(one_bin["u"]))
            @test ismissing(one_bin["leak"])

            result = surd([0.0, 0.34, 1.0], ([0.0, 0.0, 1.0],), SURD.Histogram(; nbins=3))
            @test result["leak"] ≈ (2 * log(2)) / (3 * log(3))

            narrow = [1.0, nextfloat(1.0)]
            @test surd(narrow, (narrow,), SURD.Histogram(4))["u"]["1"] ≈ 1.0
            @test ismissing(surd(zeros(16), agents, SURD.Histogram(4))["leak"])
            @test_throws ArgumentError SURD.Histogram(0)
        end
    end

    @testset "kNN" begin
        @testset "Three agents and shared radius" begin
            y = [3, 7, 7, 8, 6, 2, 2, 0]
            x1 = [5, 4, 3, 0, 2, 6, 7, 2]
            x2 = [1, 3, 8, 7, 5, 9, 4, 7]
            x3 = [4, 0, 1, 1, 8, 3, 8, 1]
            result = surd(y, (x1, x2, x3), SURD.kNN(; k=3); log=true)

            # Reference values from exact harmonic-number scores and pairwise
            # distances on the eight rows, including points on the radius.
            I_total = 2333//13440
            expected = (
                ("u", Dict("1" => 17//630, "2" => 0, "3" => 53//3360)),
                ("r", Dict("1,2" => 0, "1,3" => 1261//40320, "2,3" => 0, "1,2,3" => 0)),
                ("s", Dict("1,2" => 0, "1,3" => 2203//40320, "2,3" => 3//160, "1,2,3" => 211//8064)),
            )
            for (label, components) in expected
                @test keys(result[label]) == keys(components)
                for (key, value) in components
                    @test result[label][key] ≈ value / I_total atol=1e-12
                end
            end
            @test Set(keys(result)) == Set(("u", "r", "s", "leak"))
            @test sum(sum(values(result[label])) for label in ("u", "r", "s")) ≈ 1.0
            @test result["leak"] ≈ exp(-2 * I_total)
        end

        @testset "Empty Y neighborhood and zero total" begin
            # At the first observation both other Y values lie on the radius.
            # Its fallback must use the lower row index. The third observation
            # uses both scores, making all three corrected scores 1/2.
            result = surd([0, 1, -1], ([0, 0, 3],), SURD.kNN(1))
            @test result["u"]["1"] ≈ 1.0
            @test result["leak"] ≈ exp(-1)
            @test isempty(result["r"]) && isempty(result["s"])

            result = surd([0.0, 1.0], ([0.0, 0.5],), SURD.kNN(1))
            @test result["u"]["1"] ≈ 0.0 atol=1e-12
            @test result["leak"] == 1.0

            # The third observation averages a local score of -1/2. The
            # other two averages are zero, so the total after correction is zero.
            result = surd([0, 0, 1], ([0, 2, 0],), SURD.kNN(1))
            @test result["u"]["1"] ≈ 0.0 atol=1e-12
            @test result["leak"] ≈ 1.0
        end

        @testset "Vector-valued observations and input checks" begin
            y = [7 1 1 7 1 6; 1 7 7 4 5 7]
            x = [0 4 4 2 7 0; 6 2 2 5 6 5]
            result = surd(y, (x,), SURD.kNN(2))
            @test result["u"]["1"] ≈ 1.0
            @test result["leak"] ≈ exp(-187/360)

            @test_throws ArgumentError surd([0, 1], ([0, 1],), SURD.kNN(0))
            @test_throws ArgumentError surd([0, 1], ([0, 1],), SURD.kNN(2))
            @test_throws ArgumentError surd([0, 1], ([0],), SURD.kNN(1))
            @test_throws ArgumentError surd([0, 0], ([0, 0],), SURD.kNN(1))
        end
    end
end
