@testset "explicit remainder initialization and coupling modes" begin
    x = [-0.75, -0.25, 0.25, 0.75]
    z, v, q, τ = 0.1, -0.2, 0.6, 0.4
    regular(x) = (0.3 + 0.1x, -0.25, 0.2)
    full(x) = begin
        r = regular(x)
        (SF.Φs(x, z, v, q, τ) + r[1],
         SF.Πs(x, z, v, q) + r[2],
         SF.Ψs(x, z, v, q) + r[3])
    end

    zero_data = SF.initialize_remainder(x; mode=:zero_remainder, q=q)
    @test zero_data == zeros(4, 3)
    retarded = SF.initialize_remainder(x; mode=:retarded, z, v, q, τ,
                                       full_field=full)
    @test retarded[:, 1] ≈ fill(-0.25, 4)
    @test retarded[:, 2] ≈ fill(0.2, 4)
    @test retarded[:, 3] ≈ 0.3 .+ 0.1 .* x
    custom = reshape(collect(1.0:12.0), 4, 3)
    @test SF.initialize_remainder(x; mode=:custom, custom) == custom
    @test_throws ArgumentError SF.initialize_remainder(x; mode=:retarded, q=q)

    m0, Φr = 2.0, 0.3
    @test ODE.dynamical_mass(m0, q, Φr, τ) ≈ m0 + q * Φr - q^2 * τ / 2
    rates = ODE.particle_rates(m0, z, v, q, -0.25, 0.2, z -> 0.0;
                               particle_feedback=false)
    @test rates.dm == 0
    @test rates.dv == 0

    grid, h = centered_grid(1.0, 20)
    base = (equation=EQ.LinearScalarWaveEquation1D(1.0),
            numerical_flux=NF.flux_godunov, slope_limiter=LIM.MC,
            x=grid, h, N=length(grid), q, z0=0.0, A=0.1, ω=1.2)
    d_off = zeros(20, 3)
    d_on = zeros(20, 3)
    Q0 = zeros(20, 3)
    ODE.field_rhs_forced!(d_off, Q0, merge(base, (field_source=false,)), 0.5)
    ODE.field_rhs_forced!(d_on, Q0, merge(base, (field_source=true,)), 0.5)
    @test d_off == zeros(20, 3)
    @test any(x -> !iszero(x), d_on)
end
