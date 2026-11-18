@testset "relativistic particle in a constant mechanical force" begin
    m, z0, v0, F, tf = 1.8, -0.2, 0.35, 0.4, 1.3
    x = [-1.5, -0.5, 0.5, 1.5]
    Qfield = zeros(length(x), 3)
    u0 = ArrayPartition(Qfield, [m, z0, v0, 0.0])
    params = (q=0.0, dxPotential=z -> -F, sf=true, x=x, h=1.0, N=length(x),
              equation=EQ.LinearScalarWaveEquation1D(1.0), slope_limiter=minmod)

    function rhs!(du, u, p, t)
        fill!(du.x[1], 0)
        ODE.particle_rhs_relativistic!(du, u, p, t)
    end
    sol = solve(ODEProblem(rhs!, u0, (0.0, tf), params), Tsit5();
                abstol=1e-12, reltol=1e-12)
    numerical = sol.u[end].x[2]

    p0 = ref_gamma(v0) * m * v0
    pt = p0 + F * tf
    expected_v = pt / sqrt(m^2 + pt^2)
    expected_z = z0 + (sqrt(m^2 + pt^2) - sqrt(m^2 + p0^2)) / F
    @test numerical[1] ≈ m atol=2e-11
    @test numerical[2] ≈ expected_z atol=2e-10
    @test numerical[3] ≈ expected_v atol=2e-10
end
