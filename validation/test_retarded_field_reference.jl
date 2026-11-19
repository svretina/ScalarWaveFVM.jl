using Test

include("RetardedFieldReference.jl")
using .RetardedFieldReference

@testset "independent retarded-field reference" begin
    t = collect(range(0.0, 1.0; length=21))
    z = fill(0.2, length(t))
    v = zeros(length(t))
    τ = copy(t)
    history = TrajectoryHistory(t, z, v, τ)
    q = 0.7

    for time in (0.0, 0.23, 0.8), x in (-0.9, 0.0, 0.35, 1.1)
        remainder = exact_remainder_1d(history, time, x, q)
        @test remainder.Φ ≈ 0 atol=2e-14
        @test remainder.Π ≈ 0 atol=2e-14
        @test remainder.Ψ ≈ 0 atol=2e-14
    end

    for center in (-0.31, 0.19, 0.72)
        average = cell_average_exact_remainder_1d(history, 0.8, center, 0.1, q)
        @test average.Φ ≈ 0 atol=2e-14
        @test average.Π ≈ 0 atol=2e-14
        @test average.Ψ ≈ 0 atol=2e-14
    end

    field = retarded_scalar_1d(history, 0.8, 0.5, q)
    @test field.tret !== nothing
    @test 0.8 - field.tret ≈ abs(0.5 - trajectory_state(history, field.tret).z)
end
