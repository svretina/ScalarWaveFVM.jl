function isolated_coupled_run(N; tf=6.0)
    x, h = centered_grid(4.0, N)
    m0, z0, v0, q = 3.0, 0.17, 0.35, 0.4
    u0 = ArrayPartition(zeros(N, 3), [m0, z0, v0, 0.0])
    params = (q, m0, dxPotential=z -> 0.0, particle_feedback=true,
              field_source=true, x, h, N,
              equation=EQ.LinearScalarWaveEquation1D(1.0),
              numerical_flux=NF.flux_godunov, slope_limiter=LIM.MC)
    sol = solve(ODEProblem(ODE.radiation_reaction_rhs_relativistic!, u0,
                           (0.0, tf), params), SSPRK54();
                adaptive=false, dt=0.2h, saveat=range(0.0, tf; length=31))
    return sol, params
end

@testset "long-time isolated coupled health" begin
    sol1, params1 = isolated_coupled_run(80)
    sol2, params2 = isolated_coupled_run(160)
    @test sol1.retcode == ReturnCode.Success
    @test sol2.retcode == ReturnCode.Success
    health = DIAG.runtime_health(sol2.u)
    @test health.valid
    @test health.minimum_mass > 2.5
    @test health.maximum_speed ≈ 0.35 atol=2e-13

    p1, p2 = sol1.u[end].x[2], sol2.u[end].x[2]
    @test p2[2] ≈ 0.17 + 0.35 * 6 atol=2e-12
    @test p2[3] ≈ 0.35 atol=2e-13
    @test p2[1] ≈ 3.0 - 0.4^2 * p2[4] / 2 atol=2e-12
    @test p1[2:4] ≈ p2[2:4] atol=3e-12

    obs = DIAG.particle_observables(sol2.u[end].x[1], p2, params2)
    @test obs.mass_constraint ≈ 0 atol=3e-12
    balance = DIAG.remainder_energy_balance((; sol=sol1))
    @test maximum(abs, balance.residual) ≈ 0 atol=1e-14
end
