module ODE

using ..Reconstructions
using ..ParticleMotion
using ..ScalarField
using StaticArrays
using Polyester
using DataInterpolations

function trapezoidal(y, h)
    n = length(y) - 1
    integral = zero(Float64)
    @inbounds for i in 1:n
        integral += y[i]
    end
    @inbounds integral += 0.5 * (y[1] + y[n + 1])  # Adjust endpoints
    integral *= h
    return integral
end

# LaxWendroff schemes
function LxW!(dQ, Q, params, t)
    equation = params.equation
    # numerical_flux = params.numerical_flux
    flux_limiter = params.flux_limiter
    h = params.h
    N = params.N
    x = params.x
    cfl = params.cfl
    c = equation.velocity
    dt = params.dt

    half_h = h * 0.5
    _h = 1.0 / h
    r1 = @SVector [c, one(c)]
    r2 = @SVector [-c, one(c)]
    λp = c
    λm = -c
    _2c = 1.0 / (2c)
    @fastmath @inbounds begin
        @batch for i in 1:(N + 1)
            ## Left Face
            ΔΠl = (Q[i, 1] - Q[i - 1, 1]) * _2c
            ΔΨl = (Q[i, 2] - Q[i - 1, 2]) * 0.5
            a1l = ΔΠl + ΔΨl
            a2l = -ΔΠl + ΔΨl
            # W1l = a1l * r1 # left moving at left face
            W2l = a2l * r2 # right moving at left face

            ## Right Face
            ΔΠr = (Q[i + 1, 1] - Q[i, 1]) * _2c
            ΔΨr = (Q[i + 1, 2] - Q[i, 2]) * 0.5
            a1r = ΔΠr + ΔΨr
            a2r = -ΔΠr + ΔΨr
            W1r = a1r * r1 # left moving at right face
            # W2r = a2r * r2 # right moving at right face

            ### Higher resolution stuff - Limiters
            ## Left Face
            # for λ<0
            I = i + 1
            ΔΠlp = (Q[I, 1] - Q[I - 1, 1]) * _2c
            ΔΨlp = (Q[I, 2] - Q[I - 1, 2]) * 0.5
            a1lp = ΔΠlp + ΔΨlp
            a1lp == a1l ? θ1l = 1.0 : θ1l = a1lp / a1l

            a1l_tilde = flux_limiter(θ1l) * a1l
            W1l_tilde = a1l_tilde * r1
            # for λ>0
            I = i - 1
            ΔΠlm = (Q[I, 1] - Q[I - 1, 1]) * _2c
            ΔΨlm = (Q[I, 2] - Q[I - 1, 2]) * 0.5
            a2lm = -ΔΠlm + ΔΨlm
            a2lm == a2l ? θ2l = 1.0 : θ2l = a2lm / a2l

            a2l_tilde = flux_limiter(θ2l) * a2l
            W2l_tilde = a2l_tilde * r2

            ## Right Face
            # for λ<0
            I = i + 1
            ΔΠrp = (Q[I + 1, 1] - Q[I, 1]) * _2c
            ΔΨrp = (Q[I + 1, 2] - Q[I, 2]) * 0.5
            a1rp = ΔΠrp + ΔΨrp

            a1rp == a1r ? θ1r = 1.0 : θ1r = a1rp / a1r

            a1r_tilde = flux_limiter(θ1r) * a1r
            W1r_tilde = a1r_tilde * r1
            # for λ>0
            I = i - 1
            ΔΠrm = (Q[I + 1, 1] - Q[I, 1]) * _2c
            ΔΨrm = (Q[I + 1, 2] - Q[I, 2]) * 0.5
            a2rm = -ΔΠrm + ΔΨrm
            a2rm == a2r ? θ2r = 1.0 : θ2r = a2rm / a2r
            a2r_tilde = flux_limiter(θ2r) * a2r
            W2r_tilde = a2r_tilde * r2

            tmp = λp * (1.0 - cfl)

            Fl_tilde = tmp * (W1l_tilde + W2l_tilde)
            Fr_tilde = tmp * (W1r_tilde + W2r_tilde)

            dQ[i, :] .= -(λm * W1r + λp * W2l) * _h - 0.5(Fr_tilde - Fl_tilde) * _h
        end
    end
    return nothing
end

@inline function get_stencil(Q, i, N)
    @inbounds begin
        if i == 0
            qi = Q[N]
            ql = Q[N - 1]
            qr = Q[1]
        elseif i == 1
            qi = Q[1]
            ql = Q[N]
            qr = Q[2]
        elseif i == N
            qi = Q[N]
            ql = Q[N - 1]
            qr = Q[1]
        elseif i == N + 1
            qi = Q[1]
            ql = Q[N]
            qr = Q[2]
        else
            qi = Q[i]
            ql = Q[i - 1]
            qr = Q[i + 1]
        end
    end
    return qi, ql, qr
end

function muscl!(dQ, Q, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
    slope_limiter = params.slope_limiter

    h = params.h
    N = params.N

    half_h = 0.5h
    _h = 1.0 / h

    @fastmath @inbounds begin
        Π = @view Q[:, 1]
        Ψ = @view Q[:, 2]
        @batch for i in 1:N
            # Calculating slopes for Π
            I = i - 1
            qi, ql, qr = get_stencil(Π, I, N)
            σ_im1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i + 1
            qi, ql, qr = get_stencil(Π, I, N)
            σ_ip1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i
            qi, ql, qr = get_stencil(Π, I, N)
            σ_i = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            # reconstruct Π values at faces
            πlr = qi - σ_i * half_h
            πrl = qi + σ_i * half_h
            πll = ql + σ_im1 * half_h
            πrr = qr - σ_ip1 * half_h

            ## Calculating slopes for Ψ
            I = i - 1
            qi, ql, qr = get_stencil(Ψ, I, N)
            σ_im1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i + 1
            qi, ql, qr = get_stencil(Ψ, I, N)
            σ_ip1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i
            qi, ql, qr = get_stencil(Ψ, I, N)
            σ_i = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            # reconstruct Ψ values at faces
            ψlr = qi - σ_i * half_h
            ψrl = qi + σ_i * half_h
            ψll = ql + σ_im1 * half_h
            ψrr = qr - σ_ip1 * half_h

            qll = @SVector [πll, ψll]
            qlr = @SVector [πlr, ψlr]

            qrl = @SVector [πrl, ψrl]
            qrr = @SVector [πrr, ψrr]

            ## Fl: flux on the left face, accounts for incoming right going flux
            ## Fr: flux on the right face, accounts for incoming left going flux

            Fl = numerical_flux(qll, qlr, equation)
            Fr = numerical_flux(qrl, qrr, equation)

            dQ[i, 1:2] .= -(Fr - Fl) * _h
            dQ[i, 3] = Π[i]
        end
    end
    return nothing
end

function field_rhs_forced!(dQ, Q, params, t)
    x = params.x
    N = params.N
    q = params.q
    h = params.h
    L = params.L
    x10 = params.x1
    A = params.A
    ω = params.ω

    muscl!(dQ, Q, params, t)
    z, v, a = ParticleMotion.oscillator(t, x10, A, ω)
    # x2, v2, a2 = ParticleMotion.oscillator2(t, x20, vmax,
    #                                         A, direction2)
    # @inline s1(x) = a1 * ScalarField.∂vΠ(x, q1, x1, v1)
    # @inline s2(x) = a2 * ScalarField.∂vΠ(x, q2, x2, v2)
    dtΠ = @view dQ[:, 1]
    @inline s1(x) = a * ScalarField.Πs001(x, z, v, q)
    @inbounds @batch for i in 1:N
        # left_face = x[i] - h2
        # right_face = x[i] + h2
        # dtΠ[i] -= a1 * ScalarField.cell_average(left_face, right_face, q, x1, v1, h)

        dtΠ[i] -= s1(x[i])
    end
    return nothing
end

function source_rhs_forced!(dQ, Q, params, t)
    x = params.x
    N = params.N
    q = params.q
    h = params.h
    L = params.L
    x10 = params.x1
    A = params.A
    ω = params.ω

    z, v, a = ParticleMotion.oscillator(t, x10, A, ω)
    # x2, v2, a2 = ParticleMotion.oscillator2(t, x20, vmax,
    #                                         A, direction2)
    # @inline s1(x) = a1 * ScalarField.∂vΠ(x, q1, x1, v1)
    # @inline s2(x) = a2 * ScalarField.∂vΠ(x, q2, x2, v2)
    dtΠ = @view dQ[:, 1]
    @inline s1(xx) = a * ScalarField.Πs001(xx, z, v, q)
    dtΠ[i] .= -s1.(x)
    return nothing
end

function field_rhs!(dQ, Q, params, t)
    x = params.x
    N = params.N
    q = params.q

    # q2 = params.q2
    sf = params.sf
    Q_field = Q.x[1]
    Q_particle = Q.x[2]
    dQ_field = dQ.x[1]
    dQ_particle = dQ.x[2]

    a = dQ_particle[3]
    z = Q_particle[2]
    v = Q_particle[3] ### ASK ERIK dQ or Q ?
    # v1 = dQ_particle[2]
    ###
    # a2 = dQ_particle[6]
    # x2 = Q_particle[5]
    # v2 = dQ_particle[5] ##  or dQ_particle[5]??
    muscl!(dQ_field, Q_field, params, t)
    # xp1, vp1, ap1 = ParticleMotion.oscillator(t, x1, L, direction1)
    # xp2, vp2, ap2 = ParticleMotion.oscillator(t, x2, L, direction2)
    if sf
        h = params.h
        h2 = 0.5h
        @inline s1(x) = a * ScalarField.Πs001(x, z, v, q)
        @inline s2(x) = a * ScalarField.Ψs001(x, z, v, q)
        @inline s3(x) = a * ScalarField.Φs001(x, z, v, q)
        dtΠ = @view dQ_field[:, 1]
        dtΨ = @view dQ_field[:, 2]
        dtΦ = @view dQ_field[:, 3]
        @inbounds @batch for i in 1:N
            # left_face = x[i] - h2
            # right_face = x[i] + h2
            # dtΠ[i] -= a1 * ScalarField.cell_average(left_face, right_face, q, x1, v1, h)
            # old = dtΨ[i]
            ###
            dtΠ[i] -= s1(x[i])
            # dtΠ[i] -= q * sign(x[i] - z) * (1.0 + 2v^2) * a / (2(1 - v^2)^2.5)

            # dtΨ[i] -= s2(x[i])
            # @show s2(x[i])
            # @show dtΨ[i] - old
            # @show "======"
            dtΦ[i] -= s3(x[i])
        end
    end
    return nothing
end

function field_rhs2!(dQ, Q, params, t)
    x = params.x
    N = params.N
    q1 = params.q1
    q2 = params.q2
    sf = params.sf
    h = params.h
    h2 = 0.5h
    # _h = 1 / h
    Q_field = Q.x[1]
    Q_particle = Q.x[2]
    dQ_field = dQ.x[1]
    dQ_particle = dQ.x[2]

    a1 = dQ_particle[3]
    x1 = Q_particle[2]
    v1 = Q_particle[3]
    ###
    a2 = dQ_particle[6]
    x2 = Q_particle[5]
    v2 = Q_particle[6]
    muscl!(dQ_field, Q_field, params, t)
    if sf
        # @inline s1(x) = a1 * ScalarField.∂vΠ(x, q1, x1, v1)
        # @inline s2(x) = a2 * ScalarField.∂vΠ(x, q2, x2, v2)
        dtΠ = @view dQ_field[:, 1]

        @inbounds @batch for i in 1:N
            left_face = x[i] - h2
            right_face = x[i] + h2
            # subcell integration
            singular_term1_averaged = -a1 *
                                      ScalarField.cell_average(left_face, right_face,
                                                               q1, x1, v1, h)
            singular_term2_averaged = -a2 *
                                      ScalarField.cell_average(left_face, right_face,
                                                               q2, x2, v2, h)

            # dtΠ[i] -= s1(x[i]) + s2(x[i])
            dtΠ[i] += singular_term1_averaged + singular_term2_averaged
        end
    end
    return nothing
end

function particle_rhs!(dQ, Q, params, t)
    Q_field = Q.x[1]
    Q_particle = Q.x[2]
    dQ_particle = dQ.x[2]

    Π = @view Q_field[:, 1]
    Ψ = @view Q_field[:, 2]
    x = params.x
    interpolation_method = params.interpolation_method
    q = params.q

    m1 = Q_particle[1]
    x1 = Q_particle[2]
    v1 = Q_particle[3]
    abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ_particle[3])")

    interpolator_Ψ = interpolation_method(Ψ, x; extrapolation=ExtrapolationType.Extension)
    interpolator_Π = interpolation_method(Π, x; extrapolation=ExtrapolationType.Extension)
    Π1 = interpolator_Π(x1)
    Ψ1 = interpolator_Ψ(x1)

    v12 = v1 * v1
    a12 = one(v1) - v12

    dQ_particle[1] = -q * (Π1 + v1 * Ψ1)
    dQ_particle[2] = v1
    dQ_particle[3] = q * a12 * (v1 * Π1 + Ψ1) / m1
    return nothing
end

function interacting_particle_rhs!(dQ, Q, params, t)
    @inbounds begin
        Q_field = Q.x[1]
        Q_particle = Q.x[2]
        dQ_particle = dQ.x[2]

        Π = @view Q_field[:, 1]
        Ψ = @view Q_field[:, 2]
        x = params.x
        interpolation_method = params.interpolation_method
        q1 = params.q1
        q2 = params.q2

        m1 = Q_particle[1]
        x1 = Q_particle[2]
        v1 = Q_particle[3]

        m2 = Q_particle[4]
        x2 = Q_particle[5]
        v2 = Q_particle[6]

        vel = 0.5
        abs(v1) <= vel || throw("v1=$v1 > $(vel), at t=$(t), a1 = $(dQ[3])")
        abs(v2) <= vel || throw("v2=$v2 > $(vel), at t=$(t), a1 = $(dQ[6])")

        interpolator_Ψ = interpolation_method(Ψ, x;
                                              extrapolation=ExtrapolationType.Extension)
        interpolator_Π = interpolation_method(Π, x;
                                              extrapolation=ExtrapolationType.Extension)
        Πr_at_1 = interpolator_Π(x1)
        Ψr_at_1 = interpolator_Ψ(x1)
        Πr_at_2 = interpolator_Π(x2)
        Ψr_at_2 = interpolator_Ψ(x2)

        # value of singular field of first particle
        # evaluated at position of 2nd particle
        Πs1_at_2 = ScalarField.Πs(x2, q1, x1, v1)
        Ψs1_at_2 = ScalarField.Ψs(x2, q1, x1, v1)
        # value of singular field of second particle
        # evaluated at position of 1st particle
        Πs2_at_1 = ScalarField.Πs(x1, q2, x2, v2)
        Ψs2_at_1 = ScalarField.Ψs(x1, q2, x2, v2)

        v12 = v1 * v1
        a12 = one(v1) - v12
        v22 = v2 * v2
        a22 = one(v2) - v22

        Π1 = Πr_at_1 + Πs2_at_1
        Ψ1 = Ψr_at_1 + Ψs2_at_1
        Π2 = Πr_at_2 + Πs1_at_2
        Ψ2 = Ψr_at_2 + Ψs1_at_2
        # Particle 1
        dQ_particle[1] = -q1 * (Π1 + v1 * Ψ1)
        dQ_particle[2] = v1
        dQ_particle[3] = q1 * a12 * (v1 * Π1 + Ψ1) / m1
        # Particle 2
        dQ_particle[4] = -q2 * (Π2 + v2 * Ψ2)
        dQ_particle[5] = v2
        dQ_particle[6] = q2 * a22 * (v2 * Π2 + Ψ2) / m2
        # params.acc[1] = max(dQ_particle[3], dQ_particle[6])
    end
    return nothing
end

function particle_rhs_functions!(dQ, Q, params, t)
    Π = params.pifield
    Ψ = params.psifield
    q1 = params.q1

    @inbounds m1 = Q[1]
    @inbounds x1 = Q[2]
    @inbounds v1 = Q[3]
    abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ[3])")

    Π1 = Π(t, x1)
    Ψ1 = Ψ(t, x1)

    v12 = v1 * v1
    a12 = one(v1) - v12

    @inbounds dQ[1] = -q1 * (Π1 + v1 * Ψ1)
    @inbounds dQ[2] = v1
    @inbounds dQ[3] = q1 * a12 * (v1 * Π1 + Ψ1) / m1
    return nothing
end

function particle_rhs_interpolation!(dQ, Q, params, t)
    Π = params.pifield
    Ψ = params.psifield
    q1 = params.q1
    x = params.x
    interpolation_method = params.interpolation_method
    @inbounds m1 = Q[1]
    @inbounds x1 = Q[2]
    @inbounds v1 = Q[3]
    abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ[3])")
    Πvals = Π.(t, x)
    Ψvals = Ψ.(t, x)

    interpolator_Ψ = interpolation_method(Ψvals, x;
                                          extrapolation=ExtrapolationType.Extension)
    interpolator_Π = interpolation_method(Πvals, x;
                                          extrapolation=ExtrapolationType.Extension)
    Π1 = interpolator_Π(x1)
    Ψ1 = interpolator_Ψ(x1)

    v12 = v1 * v1
    a12 = one(v1) - v12

    v∇Φ1 = Π1 + v1 * Ψ1

    @inbounds dQ[1] = -q1 * v∇Φ1
    @inbounds dQ[2] = v1
    @inbounds dQ[3] = q1 * a12 * (v1 * Π1 + Ψ1) / m1
    return nothing
end

function coupled_rhs!(dQ, Q, params, t)
    particle_rhs!(dQ, Q, params, t)
    field_rhs!(dQ, Q, params, t)
    return nothing
end

function interacting_coupled_rhs!(dQ, Q, params, t)
    interacting_particle_rhs!(dQ, Q, params, t)
    field_rhs2!(dQ, Q, params, t)
    return nothing
end

function particle_in_potential_rhs!(dQ, Q, params, t)
    particle_with_potential_rhs!(dQ, Q, params, t)
    sf = params.sf
    if sf
        field_rhs!(dQ, Q, params, t)
    end
    return nothing
end

function deriv(y, h)
    n = length(y)
    dy = zeros(eltype(y), n)
    for i in 2:(n - 1)
        dy[i] = (y[i + 1] - y[i - 1]) / (2h)
    end

    dy[1] = (y[1] - y[2]) / h
    dy[n] = (y[n - 1] - y[n]) / h
    return dy
end

function particle_with_potential_rhs!(dQ, Q, params, t)
    Q_field = Q.x[1]
    Q_particle = Q.x[2]
    dQ_particle = dQ.x[2]
    dQ_field = dQ.x[1]

    q = params.q
    xs = params.x
    x0 = params.x1
    dxV = params.dxPotential
    V = params.Potential
    dt = params.dt
    sf = params.sf
    dx = params.dx
    L0 = params.L0
    N = params.N
    interpolation_method = params.interpolation_method
    slope_limiter = params.slope_limiter

    m0 = Q_particle[1]
    z = Q_particle[2]
    v = Q_particle[3]

    # @show dQ_particle[3]
    # abs(z) > 3x0 && throw("Unbound orbit! z=$(z)")

    # vel = 0.999999
    # abs(v) <= vel || throw("v=$v > $(vel), at t=$(t), a = $(dQ_particle[3])")

    vel = 1.0
    abs(tanh(v)) <= vel || throw("v=$(tanh(v)) > $(vel), at t=$(t), a = $(dQ_particle[3])")

    # a2 = one(v) - v * v
    γ = cosh(v) #1.0 / sqrt(a2)

    x_range = (-L0 / 2, L0 / 2)
    mask = (xs .>= x_range[1]) .& (xs .<= x_range[2])
    xsm = xs[mask]

    function delta(x, z, h)
        x_left = x - h / 2
        x_right = x + h / 2
        return (z >= x_left && z < x_right) ? 1.0 : 0.0
    end
    δ = delta.(xsm, z, dx)
    # dxΠr2 = zeros(N)
    # dxΨr2 = zeros(N)
    if sf
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]
        interpolator_Ψ = interpolation_method(Ψr, xsm;
                                              extrapolation=ExtrapolationType.Extension)
        interpolator_Π = interpolation_method(Πr, xsm;
                                              extrapolation=ExtrapolationType.Extension)
        interpolator_Φ = interpolation_method(Φr, xsm;
                                              extrapolation=ExtrapolationType.Extension)
        Φrz = interpolator_Φ(z)
        Πrz = interpolator_Π(z)
        Ψrz = interpolator_Ψ(z)

        # Φrz = sum(δ .* Φr) #interpolator_Φ(z)
        # Πrz = sum(δ .* Πr) #interpolator_Π(z)
        # Ψrz = sum(δ .* Ψr) #interpolator_Ψ(z)

        dtΠr = @view dQ_field[mask, 1]
        dtΨr = @view dQ_field[mask, 2]
        # dxΠr = @view dQ_field[mask, 2]
        # dxΨr = dtΠr .+ dQ_particle[3] * ScalarField.Πs001.(xsm, z, v, q)
        # dxΨr = deriv(Ψr, dx)

        @inbounds dQ_particle[1] = 0.0
        # @inbounds dQ_particle[2] = v

        @inbounds dQ_particle[2] = tanh(v) # rapidity

        # term1 = @. 0.25q^2 * sign.(xsm - z)^2 * γ^3

        # term2 = @. 0.5q * γ^5 * sign.(xsm - z) * (3Πr * v + Ψr * (1 + 2 * v^2))
        # m = -(m0 + q * Φrz) * γ^3 + trapezoidal(term2, dx)

        # temp = @. -0.5q * γ^3 * sign.(xsm - z) * (v * dtΨr + dtΠr)
        # numerator = q * v * (Ψrz * v + Πrz) * γ^3 - q * dxV(z) +
        #             trapezoidal(temp, dx)

        term1 = @. 0.25q^2 * sign.(xsm - z)^2 * γ^3

        term2 = @. 0.5q * γ^5 * sign.(xsm - z) * (3Πr * tanh(v) + Ψr * (1 + 2 * tanh(v)^2))
        m = -(m0 + q * Φrz) * γ^3 + trapezoidal(term2, dx)

        temp = @. -0.5q * γ^3 * sign.(xsm - z) * (tanh(v) * dtΨr + dtΠr)
        numerator = q * v * (Ψrz * tanh(v) + Πrz) * γ^3 - q * dxV(z) +
                    trapezoidal(temp, dx)

        @inbounds dQ_particle[3] = cosh(v)^2 * numerator / m
        params.acc[1] = dQ_particle[3]
    else
        # @inbounds dQ_particle[1] = 0.0
        # @inbounds dQ_particle[2] = v
        # @inbounds dQ_particle[3] = q * a2 * dxV(z) / (m - q * V(z))

        ## in this case the conserved energy is simply
        ## E = γ*m(t)
        # @inbounds dQ_particle[1] = 0.0
        # @inbounds dQ_particle[2] = v
        # @inbounds dQ_particle[3] = q * a2 * dxV(z) / m

        @inbounds dQ_particle[1] = 0.0
        @inbounds dQ_particle[2] = v
        @inbounds dQ_particle[3] = q * (1.0 - v * v)^1.5 * dxV(z) / m
    end
    return nothing
end

function muscl_forced!(dQ, Q, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
    slope_limiter = params.slope_limiter
    x = params.x
    N = params.N
    q = params.q
    x10 = params.x1
    A = params.A
    ω = params.ω

    h = params.h
    half_h = 0.5h
    _h = 1.0 / h

    # Get particle motion
    x1, v1, a1 = ParticleMotion.oscillator(t, x10, A, ω)

    @fastmath @inbounds begin
        Π = @view Q[:, 1]
        Ψ = @view Q[:, 2]
        dtΠ = @view dQ[:, 1]
        h2 = 0.5h

        for i in 1:N
            # Calculating slopes for Π
            I = i - 1
            qi, ql, qr = get_stencil(Π, I, N)
            σ_im1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i + 1
            qi, ql, qr = get_stencil(Π, I, N)
            σ_ip1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i
            qi, ql, qr = get_stencil(Π, I, N)
            σ_i = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            # reconstruct Π values at faces
            πlr = qi - σ_i * half_h
            πrl = qi + σ_i * half_h
            πll = ql + σ_im1 * half_h
            πrr = qr - σ_ip1 * half_h

            ## Calculating slopes for Ψ
            I = i - 1
            qi, ql, qr = get_stencil(Ψ, I, N)
            σ_im1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i + 1
            qi, ql, qr = get_stencil(Ψ, I, N)
            σ_ip1 = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            I = i
            qi, ql, qr = get_stencil(Ψ, I, N)
            σ_i = compute_limited_slope(qi, ql, qr, h, slope_limiter)

            # reconstruct Ψ values at faces
            ψlr = qi - σ_i * half_h
            ψrl = qi + σ_i * half_h
            ψll = ql + σ_im1 * half_h
            ψrr = qr - σ_ip1 * half_h

            qll = @SVector [πll, ψll]
            qlr = @SVector [πlr, ψlr]
            qrl = @SVector [πrl, ψrl]
            qrr = @SVector [πrr, ψrr]

            ## Fl: flux on the left face, accounts for incoming right going flux
            ## Fr: flux on the right face, accounts for incoming left going flux
            Fl = numerical_flux(qll, qlr, equation)
            Fr = numerical_flux(qrl, qrr, equation)

            # Calculate MUSCL flux contribution
            dQ[i, 1:2] .= -(Fr - Fl) * _h
            dQ[i, 3] = Π[i]

            # Add forced term contribution in the same loop
            left_face = x[i] - h2
            right_face = x[i] + h2
            dtΠ[i] -= a1 * ScalarField.cell_average(left_face, right_face, q, x1, v1, h)
        end
    end
    return nothing
end

function interacting_particle_harmonic_rhs!(dQ, Q, params, t)
    @inbounds begin
        Q_field = Q.x[1]
        Q_particle = Q.x[2]
        dQ_particle = dQ.x[2]

        x = params.x
        q1 = params.q1
        q2 = params.q2

        m1 = Q_particle[1]
        x1 = Q_particle[2]
        v1 = Q_particle[3]

        m2 = Q_particle[4]
        x2 = Q_particle[5]
        v2 = Q_particle[6]

        vel = 0.5
        abs(v1) <= vel || throw("v1=$v1 > $(vel), at t=$(t), a1 = $(dQ[3])")
        abs(v2) <= vel || throw("v2=$v2 > $(vel), at t=$(t), a1 = $(dQ[6])")

        Πs1(x) = ScalarField.Πs(x, q1, x1, v1)
        Ψs1(x) = ScalarField.Ψs(x, q1, x1, v1)

        Πs2(x) = ScalarField.Πs(x, q2, x2, v2)
        Ψs2(x) = ScalarField.Ψs(x, q2, x2, v2)

        v12 = v1 * v1
        a12 = one(v1) - v12
        v22 = v2 * v2
        a22 = one(v2) - v22
        γ1 = 1 / sqrt(a12)
        γ2 = 1 / sqrt(a22)
        k = 1.8e-5
        b = 1e-5
        # Particle 1
        dQ_particle[1] = 0#-q1 * (Π1 + v1 * Ψ1)
        dQ_particle[2] = v1
        dQ_particle[3] = -(k * (x1 - x2) + b * (v1 + v2)^2) /
                         (m1 * (γ1 + v1^2 / a12^(3 / 2)))
        # Particle 2
        dQ_particle[4] = 0#-q2 * (Π2 + v2 * Ψ2)
        dQ_particle[5] = v2
        dQ_particle[6] = (k * (x1 - x2) + b * (v1 + v2)^2) /
                         (m2 * (γ2 + v2^2 / a22^(3 / 2)))
    end
    return nothing
end

function interacting_particle_rhs2!(dQ, Q, params, t)
    @inbounds begin
        Q_field = Q.x[1]
        Q_particle = Q.x[2]
        dQ_particle = dQ.x[2]

        x = params.x
        q1 = params.q1
        q2 = params.q2

        m1 = Q_particle[1]
        x1 = Q_particle[2]
        v1 = Q_particle[3]

        m2 = Q_particle[4]
        x2 = Q_particle[5]
        v2 = Q_particle[6]

        vel = 0.5
        abs(v1) <= vel || throw("v1=$v1 > $(vel), at t=$(t), a1 = $(dQ[3])")
        abs(v2) <= vel || throw("v2=$v2 > $(vel), at t=$(t), a1 = $(dQ[6])")

        Πs1(x) = 0#ScalarField.Πs(x, q1, x1, 0)
        Ψs1(x) = ScalarField.Ψs(x, q1, x1, 0)

        Πs2(x) = 0#ScalarField.Πs(x, q2, x2, 0)
        Ψs2(x) = ScalarField.Ψs(x, q2, x2, 0)

        v12 = v1 * v1
        a12 = one(v1) - v12
        v22 = v2 * v2
        a22 = one(v2) - v22
        γ1 = 1 / sqrt(a12)
        γ2 = 1 / sqrt(a22)

        # Particle 1
        dQ_particle[1] = 0# -q1 * (Πs2(x1) + v1 * Ψs2(x1))
        dQ_particle[2] = v1
        dQ_particle[3] = -(q1 + q2) * sign(x1 - x2) * a12^(3 / 2) / m1
        # Particle 2
        dQ_particle[4] = 0#-q2 * (Πs1(x2) + v2 * Ψs1(x2))
        dQ_particle[5] = v2
        dQ_particle[6] = (q1 + q2) * sign(x1 - x2) * a22^(3 / 2) / m2
    end
    return nothing
end

function interacting_particle_rhs22!(dQ, Q, params, t)
    @inbounds begin
        Q_field = Q.x[1]
        Q_particle = Q.x[2]
        dQ_particle = dQ.x[2]

        x = params.x
        q1 = params.q1
        q2 = params.q2

        m1 = Q_particle[1]
        x1 = Q_particle[2]
        v1 = Q_particle[3]

        m2 = Q_particle[4]
        x2 = Q_particle[5]
        v2 = Q_particle[6]

        vel = 0.5
        abs(v1) <= vel || throw("v1=$v1 > $(vel), at t=$(t), a1 = $(dQ[3])")
        abs(v2) <= vel || throw("v2=$v2 > $(vel), at t=$(t), a1 = $(dQ[6])")

        # Πs1(x) = 0#ScalarField.Πs(x, q1, x1, 0)
        # Ψs1(x) = ScalarField.Ψs(x, q1, x1, 0)

        # Πs2(x) = 0#ScalarField.Πs(x, q2, x2, 0)
        # Ψs2(x) = ScalarField.Ψs(x, q2, x2, 0)

        v12 = v1 * v1
        a12 = one(v1) - v12
        v22 = v2 * v2
        a22 = one(v2) - v22
        γ1 = 1 / sqrt(a12)
        γ2 = 1 / sqrt(a22)

        # Particle 1
        dQ_particle[1] = 0
        dQ_particle[2] = v1
        # Particle 2
        dQ_particle[4] = 0
        dQ_particle[5] = v2

        nom1 = (x1 - x2) * (-1 + v1^2) *
               (q1 / γ1 + q2 / γ2 + q2 * v1 * v2 / γ2 - v1^2 * (q1 / γ1 + 2q2 / γ2))
        denom1 = abs(x1 - x2) * (2m1 * (-1 + v1^2) + abs(x1 - x2) * (q2 + 2q2 * v1^2)) / γ2

        nom2 = (x1 - x2) * (-1 + v2^2) *
               (q1 / γ1 + q2 / γ2 + q1 * v1 * v2 / γ1 - v2^2 * (2q1 / γ1 + q2 / γ2))
        denom2 = abs(x1 - x2) * (2m2 * (-1 + v2^2) + abs(x1 - x2) * (q1 + 2q1 * v2^2)) / γ1

        dQ_particle[3] = -nom1 / denom1
        dQ_particle[6] = nom2 / denom2
    end
    return nothing
end

end # end of module
