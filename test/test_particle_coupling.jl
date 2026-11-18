@testset "homogeneous regular scalar field" begin
    m, q, z, v = 2.2, 0.6, 0.13, -0.4
    A, B, C = 1.1, -0.3, 0.25
    x = [-1.5, -0.5, 0.5, 1.5]
    Φ = A .+ C .* x
    Qfield = hcat(fill(B, length(x)), fill(C, length(x)), Φ)
    Qparticle = [m, z, v, 0.2]
    Q = (x=(Qfield, Qparticle),)
    dQ = (x=(zeros(size(Qfield)), zeros(4)),)
    params = (q=q, dxPotential=z -> 0.0, sf=true, x=x, h=1.0, N=length(x),
              equation=EQ.LinearScalarWaveEquation1D(1.0),
              slope_limiter=(a, b) -> (a + b) / 2)

    for sample_z in (-0.49, -1e-10, 1e-10, 0.49)
        Φp, Πp, Ψp = ODE.regular_fields_at_particle(Qfield, sample_z, params)
        @test Φp ≈ A + C * sample_z
        @test Πp ≈ B
        @test Ψp ≈ C
    end

    ODE.particle_rhs_relativistic!(dQ, Q, params, 0.0)
    invγ = sqrt(1 - v^2)
    @test dQ.x[2][3] ≈ -q * invγ^2 * (C + v * B) / m
    @test dQ.x[2][1] ≈ q * (B + v * C) - q^2 * invγ / 2
    @test dQ.x[2][4] ≈ invγ
end
