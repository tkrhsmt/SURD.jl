using SURD
using Test

@testset "SURD.jl" begin

    a = rand(100000)
    b = rand(100000)

    surd_value = surd(a, (a, b); log=true)

    @test round(surd_value["leak"], digits=3) ≈ 0.0
    @test round(surd_value["u"]["1"], digits=3) ≈ 1.0
    @test round(surd_value["u"]["2"], digits=3) ≈ 0.0
    @test round(surd_value["r"]["1,2"], digits=3) ≈ 0.0
    @test round(surd_value["s"]["1,2"], digits=3) ≈ 0.0

end
