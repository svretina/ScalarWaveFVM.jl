using Test
using ForwardDiff: derivative

module ScalarWaveFVMTestHarness
include(joinpath(@__DIR__, "..", "src", "Reconstructions.jl"))
include(joinpath(@__DIR__, "..", "src", "ParticleMotion.jl"))
include(joinpath(@__DIR__, "..", "src", "ScalarField.jl"))
include(joinpath(@__DIR__, "..", "src", "RetardedPotentials.jl"))
include(joinpath(@__DIR__, "..", "src", "Equations.jl"))
include(joinpath(@__DIR__, "..", "src", "ODE.jl"))
include(joinpath(@__DIR__, "..", "src", "Energies.jl"))
end

const SF = ScalarWaveFVMTestHarness.ScalarField
const RP = ScalarWaveFVMTestHarness.RetardedPotentials
const ODE = ScalarWaveFVMTestHarness.ODE
const EQ = ScalarWaveFVMTestHarness.Equations
const EN = ScalarWaveFVMTestHarness.Energies

minmod(a, b) = signbit(a) == signbit(b) ? sign(a) * min(abs(a), abs(b)) : zero(a + b)

function constant_field_problem(; m=2.0, z=0.0, v=0.0, q=0.7,
                                Π=0.0, Ψ=0.0, Φ=0.0,
                                dV=z -> 0.0, sf=true)
    x = [-1.0, 0.0, 1.0]
    Qfield = hcat(fill(Π, 3), fill(Ψ, 3), fill(Φ, 3))
    Qparticle = [m, z, v, 0.0]
    Q = (x=(Qfield, Qparticle),)
    dQ = (x=(zeros(size(Qfield)), zeros(4)),)
    params = (q=q, dxPotential=dV, sf=sf, x=x, h=1.0, N=3,
              equation=EQ.LinearScalarWaveEquation1D(1.0),
              slope_limiter=minmod)
    return Q, dQ, params
end

@testset "ScalarWaveFVM.jl" begin
    @testset "retarded-matched 1+1 puncture" begin
        q = 1.7
        v = 0.3
        z = -0.2
        τ = 0.8
        γ = inv(sqrt(1 - v^2))

        @test SF.Φs(z, z, v, q, τ) == -q * τ / 2
        for x in (-1.1, 0.6)
            r = x - z
            n = sign(r)
            @test SF.Φs(x, z, v, q, τ) ≈ q * γ * (abs(r) + v * r) / 2 - q * τ / 2
            @test SF.Πs(x, z, v, q) ≈ -q * γ * (1 + v * n) / 2
            @test SF.Ψs(x, z, v, q) ≈ q * γ * (n + v) / 2
        end

        ΠL, ΠR, ΨL, ΨR = SF.singular_particle_traces(z, v, q)
        @test ΨL + v * ΠL ≈ -q / (2γ)
        @test ΨR + v * ΠR ≈ q / (2γ)
        @test ((ΨL + v * ΠL) + (ΨR + v * ΠR)) / 2 ≈ 0.0 atol=1e-14

        x = 0.6
        @test SF.dvΦs(x, z, v, q) ≈ derivative(vv -> SF.Φs(x, z, vv, q, τ), v)
        @test SF.dvΠs(x, z, v, q) ≈ derivative(vv -> SF.Πs(x, z, vv, q), v)
        @test SF.dvΨs(x, z, v, q) ≈ derivative(vv -> SF.Ψs(x, z, vv, q), v)

        # For uniform motion, the instantaneous puncture equals the exact
        # switched-on retarded potential inside the causal future.
        t, x, z0 = 1.4, 0.8, 0.1
        zt = z0 + v * t
        τt = t / γ
        tr = (t - x + z0) / (1 - v) # right side of the particle
        @test x > zt
        @test SF.Φs(x, zt, v, q, τt) ≈ -q * tr / (2γ)
    end


    @testset "exact retarded fields" begin
        q = 1.3
        stationary = (z0=0.0, A=0.0, ω=1.0, q=q)
        x, t = 0.2, 1.0
        @test RP.Φ_ret(x, stationary, t) ≈ -q * (t - abs(x)) / 2
        @test RP.Π_ret(x, stationary, t) ≈ -q / 2
        @test RP.Ψ_ret(x, stationary, t) ≈ q / 2
        @test RP.Φ_ret(2.0, stationary, t) == 0.0

        accelerated = (z0=0.0, A=0.1, ω=1.7, q=q)
        zt = RP.z(t, accelerated)
        vt = RP.v(t, accelerated)
        τt, _ = ScalarWaveFVMTestHarness.RetardedPotentials.quadgk(
            s -> sqrt(1 - RP.v(s, accelerated)^2), 0.0, t)
        ϵ = 1e-7
        for side in (-1.0, 1.0)
            xside = zt + side * ϵ
            @test RP.Φ_ret(xside, accelerated, t) - SF.Φs(xside, zt, vt, q, τt) ≈ 0.0 atol=1e-6
            @test RP.Π_ret(xside, accelerated, t) - SF.Πs(xside, zt, vt, q) ≈ 0.0 atol=1e-6
            @test RP.Ψ_ret(xside, accelerated, t) - SF.Ψs(xside, zt, vt, q) ≈ 0.0 atol=1e-6
        end
    end

    @testset "exact finite-volume puncture averages" begin
        xi, h, z, v, q, τ = 0.0, 1.0, 0.19, 0.2, 1.7, 0.4
        left, right = xi - h / 2, xi + h / 2
        lenL, lenR = z - left, right - z

        expected_Φ = (lenL * (SF.Φs(left, z, v, q, τ) + SF.Φs(z, z, v, q, τ)) +
                      lenR * (SF.Φs(z, z, v, q, τ) + SF.Φs(right, z, v, q, τ))) / (2h)
        expected_Π = (lenL * SF.Πs(left, z, v, q) + lenR * SF.Πs(right, z, v, q)) / h
        expected_Ψ = (lenL * SF.Ψs(left, z, v, q) + lenR * SF.Ψs(right, z, v, q)) / h

        @test SF.cell_average_Φs(xi, h, z, v, q, τ) ≈ expected_Φ
        @test SF.cell_average_Πs(xi, h, z, v, q) ≈ expected_Π
        @test SF.cell_average_Ψs(xi, h, z, v, q) ≈ expected_Ψ
        @test SF.cell_average_dvΦs(xi, h, z, v, q) ≈
              derivative(vv -> SF.cell_average_Φs(xi, h, z, vv, q, τ), v)
        @test SF.cell_average_dvΠs(xi, h, z, v, q) ≈
              derivative(vv -> SF.cell_average_Πs(xi, h, z, vv, q), v)
        @test SF.cell_average_dvΨs(xi, h, z, v, q) ≈
              derivative(vv -> SF.cell_average_Ψs(xi, h, z, vv, q), v)

        xs = [-1.0, 0.0, 1.0]
        @test SF.cell_averages_Φs(xs, h, z, v, q, τ) ≈
              [SF.cell_average_Φs(x, h, z, v, q, τ) for x in xs]
    end

    @testset "singular source insertion" begin
        xi, h, z, v, q, a = 0.0, 1.0, 0.19, 0.2, 1.7, -0.4
        dQ = zeros(1, 3)
        ODE.add_singular_sources!(dQ, 1, xi, h, z, v, q, a)
        @test dQ[1, 1] ≈ -a * SF.cell_average_dvΠs(xi, h, z, v, q)
        @test dQ[1, 2] ≈ -a * SF.cell_average_dvΨs(xi, h, z, v, q)
        @test dQ[1, 3] ≈ -a * SF.cell_average_dvΦs(xi, h, z, v, q)
    end

    @testset "regular-field extraction API" begin
        Q, _, params = constant_field_problem(Π=0.3, Ψ=-0.2, Φ=0.4)
        Φ, Π, Ψ = ODE.regular_fields_at_particle(Q.x[1], Q.x[2][2], params)
        @test Φ ≈ 0.4
        @test Π ≈ 0.3
        @test Ψ ≈ -0.2
    end

    @testset "relativistic particle equations" begin
        m, q = 2.5, 0.7

        Q, dQ, params = constant_field_problem(m=m, q=q)
        ODE.particle_rhs_relativistic!(dQ, Q, params, 0.0)
        @test dQ.x[2] ≈ [-q^2 / 2, 0.0, 0.0, 1.0]

        v = 0.4
        invγ = sqrt(1 - v^2)
        Q, dQ, params = constant_field_problem(m=m, v=v, q=q)
        ODE.particle_rhs_relativistic!(dQ, Q, params, 0.0)
        @test dQ.x[2][1] ≈ -q^2 * invγ / 2
        @test dQ.x[2][2] == v
        @test dQ.x[2][3] ≈ 0.0 atol=1e-14
        @test dQ.x[2][4] ≈ invγ

        Π, Ψ, dV = 0.3, -0.2, 0.6
        Q, dQ, params = constant_field_problem(m=m, v=v, q=q, Π=Π, Ψ=Ψ,
                                                dV=z -> 0.6)
        ODE.particle_rhs_relativistic!(dQ, Q, params, 0.0)
        @test dQ.x[2][1] ≈ q * (Π + v * Ψ) - q^2 * invγ / 2
        @test dQ.x[2][3] ≈
              (-q * (1 - v^2) * (Ψ + v * Π) - dV * invγ^3) / m
    end

    @testset "scalar feedback switch" begin
        v, q = 0.4, 0.7
        Q, dQ, params = constant_field_problem(v=v, q=q, Π=0.3, Ψ=-0.2, sf=false)
        ODE.particle_rhs_relativistic!(dQ, Q, params, 0.0)
        @test dQ.x[2][1] == 0.0
        @test dQ.x[2][3] == 0.0
    end

    @testset "Newton-Cotes weights" begin
        x = collect(0.0:0.25:1.0)
        @test EN.newton_cotes(x, 0.25; order=1) ≈ 0.5
        @test EN.newton_cotes(x .^ 2, 0.25; order=2) ≈ 1 / 3
    end
end
