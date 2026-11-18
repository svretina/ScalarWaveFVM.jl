using Test

module ScalarWaveFVMTestHarness
include(joinpath(@__DIR__, "..", "src", "Reconstructions.jl"))
include(joinpath(@__DIR__, "..", "src", "ParticleMotion.jl"))
include(joinpath(@__DIR__, "..", "src", "ScalarField.jl"))
include(joinpath(@__DIR__, "..", "src", "Equations.jl"))
include(joinpath(@__DIR__, "..", "src", "ODE.jl"))
end

function expected_singular_sources(x_i, half_h, _h, z, v, q, a)
    left_face = x_i - half_h
    right_face = x_i + half_h
    h = right_face - left_face

    one_v = one(v)
    inv_L = inv((one_v + v)^2)
    inv_R = inv((one_v - v)^2)

    len_L = clamp(z - left_face, zero(h), h)
    len_R = h - len_L

    src_pi_L = -0.5 * q * a * inv_L
    src_pi_R = 0.5 * q * a * inv_R
    src_psi_L = -0.5 * q * a * inv_L
    src_psi_R = -0.5 * q * a * inv_R

    upper_L = min(z, right_face)
    phi_int_L = upper_L > left_face ? z * (upper_L - left_face) -
                                      0.5 * (upper_L^2 - left_face^2) : zero(h)

    lower_R = max(z, left_face)
    phi_int_R = right_face > lower_R ? 0.5 * (right_face^2 - lower_R^2) -
                                       z * (right_face - lower_R) : zero(h)

    src_phi = _h * (0.5 * q * a * inv_L * phi_int_L - 0.5 * q * a * inv_R * phi_int_R)
    src_pi = _h * (len_L * src_pi_L + len_R * src_pi_R)
    src_psi = _h * (len_L * src_psi_L + len_R * src_psi_R)

    return src_pi, src_psi, src_phi
end

linear_preserving_limiter(a, b) = a == b ? a : zero(a)

function minmod_limiter(a, b)
    if signbit(a) == signbit(b) && !iszero(a) && !iszero(b)
        return sign(a) * min(abs(a), abs(b))
    elseif iszero(a) || iszero(b)
        return zero(a + b)
    else
        return zero(a + b)
    end
end

function singular_average_constants(v, q)
    one_v = one(v)
    denL = one_v + v
    denR = one_v - v
    invL = inv(denL^2)
    invR = inv(denR^2)
    return (; ΠL=0.5 * q * v / denL,
            ΠR=-0.5 * q * v / denR,
            ΨL=-0.5 * q / denL,
            ΨR=0.5 * q / denR,
            dΠL=0.5 * q * invL,
            dΠR=-0.5 * q * invR,
            dΨL=0.5 * q * invL,
            dΨR=0.5 * q * invR,
            invL,
            invR)
end

@testset "ScalarWaveFVM.jl" begin
    @testset "add_singular_sources! adds all analytic source terms" begin
        x_i = 0.0
        half_h = 0.5
        h = 2 * half_h
        _h = inv(h)
        v = 0.2
        q = 1.7
        a = -0.4
        idx = 2

        for z in (1.1, -1.1, 0.19, 0.5, -0.5)
            dQ = fill(1.25, 3, 3)
            expected_pi, expected_psi, expected_phi =
                expected_singular_sources(x_i, half_h, _h, z, v, q, a)

            ScalarWaveFVMTestHarness.ODE.add_singular_sources!(dQ,
                                                               idx,
                                                               x_i,
                                                               half_h,
                                                               _h,
                                                               z,
                                                               v,
                                                               q,
                                                               a)

            @test dQ[idx, 1] ≈ 1.25 + expected_pi atol=1e-12 rtol=1e-12
            @test dQ[idx, 2] ≈ 1.25 + expected_psi atol=1e-12 rtol=1e-12
            @test dQ[idx, 3] ≈ 1.25 + expected_phi atol=1e-12 rtol=1e-12
            @test dQ[1, :] == fill(1.25, 3)
            @test dQ[3, :] == fill(1.25, 3)
        end
    end

    @testset "non_radiative_particle_rhs_relativistic! matches m gamma^3 vdot = -k z" begin
        du = zeros(3)
        u = [2.0, 0.25, 0.4]
        params = (k=3.0,)

        ScalarWaveFVMTestHarness.ODE.non_radiative_particle_rhs_relativistic!(du,
                                                                              u,
                                                                              params,
                                                                              0.0)

        @test du[1] == 0.0
        @test du[2] == u[3]
        @test du[3] ≈ -params.k * u[2] * (1 - u[3]^2)^1.5 atol=1e-12 rtol=1e-12
    end

    @testset "particle_rhs_relativistic! evolves mass with gathered singular Φs term" begin
        x = [-1.0, 0.0, 1.0]
        h = 1.0
        N = length(x)
        equation = ScalarWaveFVMTestHarness.Equations.LinearScalarWaveEquation1D(1.0)
        Πr = 0.3
        Ψr = -0.2
        Φr = 0.4
        q = 1.7
        m0 = 2.5
        z = 0.0
        v = 0.25

        Q_field = [Πr Ψr Φr
                   Πr Ψr Φr
                   Πr Ψr Φr]
        Q_particle = [m0, z, v]
        Q = (x=(Q_field, Q_particle),)
        dQ = (x=(zeros(size(Q_field)), zeros(3)),)
        params = (q=q, m0=m0, dxPotential=z -> 0.0, x=x, h=h, N=N,
                  equation=equation, slope_limiter=minmod_limiter)

        ScalarWaveFVMTestHarness.ODE.particle_rhs_relativistic!(dQ, Q, params, 0.0)

        expected_accel = -((-1 + v^2) * (q * (v * Πr + Ψr))) / Q_particle[1]
        expected_dmdt = -q * (Πr + v * Ψr)

        @test dQ.x[2][1] ≈ expected_dmdt atol=1e-12 rtol=1e-12
        @test dQ.x[2][2] == v
        @test dQ.x[2][3] ≈ expected_accel atol=1e-12 rtol=1e-12
    end

    @testset "gather_muscl reconstructs analytic scalar data" begin
        x = [-1.5, -0.5, 0.5, 1.5]
        h = 1.0
        N = length(x)
        equation = ScalarWaveFVMTestHarness.Equations.LinearScalarWaveEquation1D(1.0)

        constant_params = (x=x, h=h, N=N, equation=equation,
                           slope_limiter=minmod_limiter)
        q_constant = fill(7.25, N)
        z_boundary = -1.9
        @test ScalarWaveFVMTestHarness.ODE.gather_muscl(q_constant,
                                                        z_boundary,
                                                        constant_params) ≈ 7.25 atol=1e-12 rtol=1e-12

        linear_params = (x=x, h=h, N=N, equation=equation,
                         slope_limiter=linear_preserving_limiter)
        f(x) = 2.0 + 3.0 * x
        q_linear = f.(x)
        z_interior = 0.25
        @test ScalarWaveFVMTestHarness.ODE.gather_muscl(q_linear,
                                                        z_interior,
                                                        linear_params) ≈ f(z_interior) atol=1e-12 rtol=1e-12
        @test ScalarWaveFVMTestHarness.ODE.gather_muscl(q_linear,
                                                        x[3],
                                                        linear_params) ≈ q_linear[3] atol=1e-12 rtol=1e-12

        q_extremum = [0.0, 1.0, 0.0, 0.0]
        z_peak_cell = -0.25
        @test ScalarWaveFVMTestHarness.ODE.gather_muscl(q_extremum,
                                                        z_peak_cell,
                                                        constant_params) ≈ 1.0 atol=1e-12 rtol=1e-12
    end

    @testset "gather_muscl on local cell-average callbacks matches array gather" begin
        x = [-1.5, -0.5, 0.5, 1.5]
        h = 1.0
        N = length(x)
        equation = ScalarWaveFVMTestHarness.Equations.LinearScalarWaveEquation1D(1.0)
        params = (x=x, h=h, N=N, equation=equation, slope_limiter=minmod_limiter)
        z = 0.19
        v = 0.2
        q = 1.7

        dvΦ_avg = ScalarWaveFVMTestHarness.ScalarField.cell_averages_dvΦs(x, h, z, v, q)
        gathered_array = ScalarWaveFVMTestHarness.ODE.gather_muscl(dvΦ_avg, z, params)
        gathered_local = ScalarWaveFVMTestHarness.ODE.gather_muscl(ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΦs,
                                                                   z,
                                                                   params,
                                                                   v,
                                                                   q)

        @test gathered_local ≈ gathered_array atol=1e-12 rtol=1e-12
    end

    @testset "gather_muscl reconstructs analytic characteristic-linear states" begin
        x = [-1.5, -0.5, 0.5, 1.5]
        h = 1.0
        N = length(x)
        equation = ScalarWaveFVMTestHarness.Equations.LinearScalarWaveEquation1D(1.0)
        params = (x=x, h=h, N=N, equation=equation,
                  slope_limiter=linear_preserving_limiter)

        w1(x) = 1.0 + 2.0 * x
        w2(x) = 2.0 - x
        qvec(x) = [w1(x) - w2(x) w1(x) + w2(x)]
        Q = reduce(vcat, [qvec(xi) for xi in x])

        z = 0.25
        gathered = ScalarWaveFVMTestHarness.ODE.gather_muscl(Q, z, params)
        expected = qvec(z)

        @test gathered[1] ≈ expected[1] atol=1e-12 rtol=1e-12
        @test gathered[2] ≈ expected[2] atol=1e-12 rtol=1e-12
    end

    @testset "exact singular cell averages" begin
        xi = 0.0
        h = 1.0
        left_face = xi - 0.5 * h
        right_face = xi + 0.5 * h
        v = 0.2
        q = 1.7
        coeffs = singular_average_constants(v, q)

        @testset "one-sided singular traces match the current convention" begin
            left_trace = ScalarWaveFVMTestHarness.ScalarField.singular_face_trace(zero(xi),
                                                                                  zero(xi),
                                                                                  v,
                                                                                  q,
                                                                                  -one(v))
            right_trace = ScalarWaveFVMTestHarness.ScalarField.singular_face_trace(zero(xi),
                                                                                   zero(xi),
                                                                                   v,
                                                                                   q,
                                                                                   one(v))

            @test left_trace[1] ≈ coeffs.ΠL atol=1e-12 rtol=1e-12
            @test left_trace[2] ≈ coeffs.ΨL atol=1e-12 rtol=1e-12
            @test right_trace[1] ≈ coeffs.ΠR atol=1e-12 rtol=1e-12
            @test right_trace[2] ≈ coeffs.ΨR atol=1e-12 rtol=1e-12

            particle_traces = ScalarWaveFVMTestHarness.ScalarField.singular_particle_traces(zero(xi),
                                                                                             v,
                                                                                             q)
            @test particle_traces[1] ≈ coeffs.ΠL atol=1e-12 rtol=1e-12
            @test particle_traces[2] ≈ coeffs.ΠR atol=1e-12 rtol=1e-12
            @test particle_traces[3] ≈ coeffs.ΨL atol=1e-12 rtol=1e-12
            @test particle_traces[4] ≈ coeffs.ΨR atol=1e-12 rtol=1e-12
        end

        @testset "field averages: whole-cell left/right constants" begin
            z_left = 1.1
            z_right = -1.1

            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Πs(xi,
                                                                       h,
                                                                       z_left,
                                                                       v,
                                                                       q) ≈ coeffs.ΠL atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Πs(xi,
                                                                       h,
                                                                       z_right,
                                                                       v,
                                                                       q) ≈ coeffs.ΠR atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Ψs(xi,
                                                                       h,
                                                                       z_left,
                                                                       v,
                                                                       q) ≈ coeffs.ΨL atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Ψs(xi,
                                                                       h,
                                                                       z_right,
                                                                       v,
                                                                       q) ≈ coeffs.ΨR atol=1e-12 rtol=1e-12
        end

        @testset "field averages: split-cell weighted and exact formulas" begin
            z_split = 0.19
            lenL = z_split - left_face
            lenR = right_face - z_split
            expected_Π = (lenL * coeffs.ΠL + lenR * coeffs.ΠR) / h
            expected_Ψ = (lenL * coeffs.ΨL + lenR * coeffs.ΨR) / h
            expected_Φ_split = 0.25 * q *
                               (lenL^2 / (one(v) + v) + lenR^2 / (one(v) - v)) / h

            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Πs(xi,
                                                                       h,
                                                                       z_split,
                                                                       v,
                                                                       q) ≈ expected_Π atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Ψs(xi,
                                                                       h,
                                                                       z_split,
                                                                       v,
                                                                       q) ≈ expected_Ψ atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Φs(xi,
                                                                       h,
                                                                       z_split,
                                                                       v,
                                                                       q) ≈ expected_Φ_split atol=1e-12 rtol=1e-12
        end

        @testset "field averages: whole-cell linear Φ formulas" begin
            z_left = 1.1
            z_right = -1.1

            expected_Φ_left = 0.5 * q * (z_left - xi) / (one(v) + v)
            expected_Φ_right = 0.5 * q * (xi - z_right) / (one(v) - v)

            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Φs(xi,
                                                                       h,
                                                                       z_left,
                                                                       v,
                                                                       q) ≈ expected_Φ_left atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_Φs(xi,
                                                                       h,
                                                                       z_right,
                                                                       v,
                                                                       q) ≈ expected_Φ_right atol=1e-12 rtol=1e-12
        end

        @testset "derivative averages: whole-cell left/right constants" begin
            z_left = 1.1
            z_right = -1.1

            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΠs(xi,
                                                                         h,
                                                                         z_left,
                                                                         v,
                                                                         q) ≈ coeffs.dΠL atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΠs(xi,
                                                                         h,
                                                                         z_right,
                                                                         v,
                                                                         q) ≈ coeffs.dΠR atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΨs(xi,
                                                                         h,
                                                                         z_left,
                                                                         v,
                                                                         q) ≈ coeffs.dΨL atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΨs(xi,
                                                                         h,
                                                                         z_right,
                                                                         v,
                                                                         q) ≈ coeffs.dΨR atol=1e-12 rtol=1e-12
        end

        @testset "derivative averages: split-cell weighted and exact formulas" begin
            z_split = 0.19
            lenL = z_split - left_face
            lenR = right_face - z_split
            expected_dΠ = (lenL * coeffs.dΠL + lenR * coeffs.dΠR) / h
            expected_dΨ = (lenL * coeffs.dΨL + lenR * coeffs.dΨR) / h
            expected_dΦ = 0.25 * q * (-coeffs.invL * lenL^2 + coeffs.invR * lenR^2) / h

            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΠs(xi,
                                                                         h,
                                                                         z_split,
                                                                         v,
                                                                         q) ≈ expected_dΠ atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΨs(xi,
                                                                         h,
                                                                         z_split,
                                                                         v,
                                                                         q) ≈ expected_dΨ atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΦs(xi,
                                                                         h,
                                                                         z_split,
                                                                         v,
                                                                         q) ≈ expected_dΦ atol=1e-12 rtol=1e-12
        end

        @testset "derivative averages: whole-cell linear dvΦ formulas" begin
            z_left = 1.1
            z_right = -1.1
            expected_dΦ_left = -0.5 * q * coeffs.invL * (z_left - xi)
            expected_dΦ_right = 0.5 * q * coeffs.invR * (xi - z_right)

            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΦs(xi,
                                                                         h,
                                                                         z_left,
                                                                         v,
                                                                         q) ≈ expected_dΦ_left atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΦs(xi,
                                                                         h,
                                                                         z_right,
                                                                         v,
                                                                         q) ≈ expected_dΦ_right atol=1e-12 rtol=1e-12
        end

        @testset "wrapper smoke tests" begin
            x = [-1.5, -0.5, 0.5]
            z = 0.1

            @test ScalarWaveFVMTestHarness.ScalarField.cell_averages_Φs(x, h, z, v, q) ≈
                  [ScalarWaveFVMTestHarness.ScalarField.cell_average_Φs(xi_, h, z, v, q)
                   for xi_ in x] atol=1e-12 rtol=1e-12
            @test ScalarWaveFVMTestHarness.ScalarField.cell_averages_dvΨs(x, h, z, v, q) ≈
                  [ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΨs(xi_, h, z, v, q)
                   for xi_ in x] atol=1e-12 rtol=1e-12
        end

        @testset "source consistency with add_singular_sources!" begin
            a = -0.4
            half_h = 0.5 * h
            _h = inv(h)

            for z in (1.1, -1.1, 0.19)
                dQ = zeros(1, 3)
                ScalarWaveFVMTestHarness.ODE.add_singular_sources!(dQ,
                                                                   1,
                                                                   xi,
                                                                   half_h,
                                                                   _h,
                                                                   z,
                                                                   v,
                                                                   q,
                                                                   a)

                @test dQ[1, 1] ≈ -a *
                      ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΠs(xi,
                                                                             h,
                                                                             z,
                                                                             v,
                                                                             q) atol=1e-12 rtol=1e-12
                @test dQ[1, 2] ≈ -a *
                      ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΨs(xi,
                                                                             h,
                                                                             z,
                                                                             v,
                                                                             q) atol=1e-12 rtol=1e-12
                @test dQ[1, 3] ≈ -a *
                      ScalarWaveFVMTestHarness.ScalarField.cell_average_dvΦs(xi,
                                                                             h,
                                                                             z,
                                                                             v,
                                                                             q) atol=1e-12 rtol=1e-12
            end
        end
    end
end
