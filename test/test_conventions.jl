@testset "convention and distributional jumps" begin
    q, z = 1.3, 0.17
    for v in (-0.9, -0.7, -0.3, 0.0, 0.3, 0.7, 0.9)
        γ = ref_gamma(v)
        ΠL, ΠR, ΨL, ΨR = SF.singular_particle_traces(z, v, q)

        @test ΨR - ΨL ≈ q * γ
        @test ΠR - ΠL ≈ -q * γ * v
        delta_coefficient = -v * (ΠR - ΠL) - (ΨR - ΨL)
        @test delta_coefficient ≈ -q / γ
        @test -delta_coefficient ≈ q / γ # coefficient in □Φ

        KL = ΨL + v * ΠL
        KR = ΨR + v * ΠR
        @test KL ≈ -q / (2γ)
        @test KR ≈ q / (2γ)
        @test (KL + KR) / 2 ≈ 0 atol=2e-14
    end
end
