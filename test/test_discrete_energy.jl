@testset "exact semidiscrete energy identity" begin
    N, L, q = 48, 2.0, 0.7
    x, h = centered_grid(L, N)
    z, v, a = 0.13, -0.27, 0.31
    Π = @. 0.2 + 0.1sin(2π * x / (2L))
    Ψ = @. -0.15 + 0.08cos(4π * x / (2L))
    Q = hcat(Π, Ψ, zeros(N))
    source = zeros(N, 3)
    for i in 1:N
        ODE.add_singular_sources!(source, i, x[i], h, z, v, q, a)
    end

    for numerical_flux in (NF.flux_godunov, NF.flux_rusanov),
        slope_limiter in (LIM.upwind, LIM.MC)
        params = (equation=EQ.LinearScalarWaveEquation1D(1.0), numerical_flux,
                  slope_limiter, h, N, x)
        dQ = similar(Q)
        ODE.muscl_optimized!(dQ, Q, params, 0.0)
        dQ .+= source
        direct_rate = h * sum(@view(Q[:, 1]) .* @view(dQ[:, 1]) .+
                              @view(Q[:, 2]) .* @view(dQ[:, 2]))
        terms = DIAG.semidiscrete_energy_terms(Q, source, params)
        @test direct_rate ≈ terms.rate atol=2e-13
        @test terms.dissipation >= 0

        cells = 9:37
        direct_subdomain_rate = h * sum(
            @view(Q[cells, 1]) .* @view(dQ[cells, 1]) .+
            @view(Q[cells, 2]) .* @view(dQ[cells, 2]))
        subdomain_terms = DIAG.semidiscrete_energy_terms(Q, source, params,
                                                         cells)
        @test direct_subdomain_rate ≈ subdomain_terms.rate atol=2e-13
        @test subdomain_terms.boundary_flux ≈
              subdomain_terms.boundary_flux_left +
              subdomain_terms.boundary_flux_right
        if slope_limiter === LIM.upwind
            @test terms.reconstruction_residual ≈ 0 atol=2e-14
            @test subdomain_terms.reconstruction_residual ≈ 0 atol=2e-14
        end
    end
end
