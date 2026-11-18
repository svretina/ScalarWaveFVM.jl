@testset "exact uniformly moving charge" begin
    q, m0, t = 0.9, 2.7, 1.1
    for v in (0.0, -0.9, -0.7, -0.3, 0.3, 0.7, 0.9)
        γ = ref_gamma(v)
        τ = t / γ
        for x in (-0.8t, -0.2t, 0.35t, 0.8t)
            x == v * t && continue
            r = x - v * t
            n = sign(r)
            reference = ref_puncture(x, v * t, v, q, τ)
            @test SF.Φs(x, v * t, v, q, τ) ≈ reference.Φ
            @test SF.Πs(x, v * t, v, q) ≈ reference.Π
            @test SF.Ψs(x, v * t, v, q) ≈ reference.Ψ
            @test reference.Φ ≈ q * γ * (abs(r) + v * r) / 2 - q * t / (2γ)
            @test reference.Π ≈ -q * γ * (1 + v * n) / 2
            @test reference.Ψ ≈ q * γ * (n + v) / 2
        end

        uL = q^2 * γ^2 * (1 - v)^2 / 4
        uR = q^2 * γ^2 * (1 + v)^2 / 4
        pL = -q^2 * γ^2 * (1 - v)^2 / 4
        pR = q^2 * γ^2 * (1 + v)^2 / 4
        Efield = uL * (v + 1) * t + uR * (1 - v) * t
        Pfield = pL * (v + 1) * t + pR * (1 - v) * t
        m = m0 - q^2 * t / (2γ)

        @test Efield ≈ q^2 * t / 2
        @test Pfield ≈ q^2 * v * t / 2
        @test γ * m + Efield ≈ γ * m0
        @test γ * m * v + Pfield ≈ γ * m0 * v

        Q, dQ, params = constant_field_problem(m=m0, v=v, q=q)
        ODE.particle_rhs_relativistic!(dQ, Q, params, 0.0)
        @test dQ.x[2][1] ≈ -q^2 / (2γ)
        @test dQ.x[2][3] ≈ 0 atol=2e-14
    end
end
