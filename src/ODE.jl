module ODE

using ..Reconstructions
using ..ParticleMotion
using ..ScalarField
using StaticArrays
using Polyester
using DataInterpolations

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

            dQ[i, :] .= -(Fr - Fl) * _h
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
    x1, v1, a1 = ParticleMotion.oscillator(t, x10, A, ω)
    # x2, v2, a2 = ParticleMotion.oscillator2(t, x20, vmax,
    #                                         A, direction2)
    # @inline s1(x) = a1 * ScalarField.∂vΠ(x, q1, x1, v1)
    # @inline s2(x) = a2 * ScalarField.∂vΠ(x, q2, x2, v2)
    dtΠ = @view dQ[:, 1]
    h2 = 0.5h
    @inbounds for i in 1:N
        left_face = x[i] - h2
        right_face = x[i] + h2
        dtΠ[i] -= a1 * ScalarField.cell_average(left_face, right_face, q, x1, v1, h)
    end
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

    a1 = dQ_particle[3]
    x1 = Q_particle[2]
    v1 = dQ_particle[2] ### ASK ERIK dQ or Q ?
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
        # @inline s1(x) = a1 * ScalarField.∂vΠ(x, q1, x1, v1)
        dtΠ = @view dQ_field[:, 1]
        @inbounds for i in 1:N
            left_face = x[i] - h2
            right_face = x[i] + h2
            dtΠ[i] -= a1 * ScalarField.cell_average(left_face, right_face, q, x1, v1, h)
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
    v1 = dQ_particle[2] ### ASK ERIK dQ or Q ?
    ###
    a2 = dQ_particle[6]
    x2 = Q_particle[5]
    v2 = dQ_particle[5] ##  or dQ_particle[5] or Q_particle[3]??
    muscl!(dQ_field, Q_field, params, t)
    if sf
        @inline s1(x) = a1 * ScalarField.∂vΠ(x, q1, x1, v1)
        @inline s2(x) = a2 * ScalarField.∂vΠ(x, q2, x2, v2)
        dtΠ = @view dQ_field[:, 1]

        @inbounds for i in 1:N
            # left_face = x[i] - h2
            # right_face = x[i] + h2
            # # subcell integration
            # singular_term1_averaged = a1 *
            #                           ScalarField.cell_average(left_face, right_face,
            #                                                    q1, x1, v1, h)
            # singular_term2_averaged = a2 *
            #                           ScalarField.cell_average(left_face, right_face,
            #                                                    q2, x2, v2, h)

            dtΠ[i] -= s1(x[i]) + s2(x[i])
            # dtΠ[i] -= singular_term1_averaged + singular_term2_averaged
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

    # v∇Φ1 = Π1 + v1 * Ψ1

    dQ_particle[1] = -q * (Π1 + v1 * Ψ1)
    dQ_particle[2] = v1
    dQ_particle[3] = q * a12 * (v1 * Π1 + Ψ1) / m1
    return nothing
end

function interacting_particle_rhs!(dQ, Q, params, t)
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

    abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ[3])")
    abs(v2) <= 1 || throw("v2=$v2 > 1, at t=$(t), a1 = $(dQ[6])")
    # vel = 0.9
    # v1 = max(-vel, min(vel, v1))
    # Q_particle[3] = v1
    # v2 = max(-vel, min(vel, v2))
    # Q_particle[6] = v2

    interpolator_Ψ = interpolation_method(Ψ, x; extrapolation=ExtrapolationType.Extension)
    interpolator_Π = interpolation_method(Π, x; extrapolation=ExtrapolationType.Extension)
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
    params.acc[1] = max(dQ_particle[3], dQ_particle[6])
    return nothing
end

function particle_rhs_functions!(dQ, Q, params, t)
    Π = params.pifield
    Ψ = params.psifield
    q1 = params.q1

    @inbounds m1 = Q[1]
    @inbounds x1 = Q[2]
    @inbounds v1 = Q[3]
    # abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ[3])")
    v1 = max(-1.0, min(1.0, v1))
    Q[3] = v1
    Π1 = Π(t, x1)
    Ψ1 = Ψ(t, x1)

    v12 = v1 * v1
    a12 = one(v1) - v12

    @inbounds dQ[1] = -q1 * (Π1 + v1 * Ψ1)
    @inbounds dQ[2] = v1
    @inbounds dQ[3] = q1 * a12 * (v1 * Π1 + Ψ1) #/ m1
    return nothing
end

function particle_rhs_interpolation!(dQ, Q, params, t)
    Π = params.pifield
    Ψ = params.psifield
    q1 = params.q1
    x = params.x
    interpolation_method = params.interpolation_method
    m1 = Q[1]
    x1 = Q[2]
    v1 = Q[3]
    abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ_particle[3])")

    interpolator_Ψ = interpolation_method(Ψ, x; extrapolation=ExtrapolationType.Extension)
    interpolator_Π = interpolation_method(Π, x; extrapolation=ExtrapolationType.Extension)
    Π1 = interpolator_Π(x1)
    Ψ1 = interpolator_Ψ(x1)

    v12 = v1 * v1
    a12 = one(v1) - v12

    v∇Φ1 = Π1 + v1 * Ψ1

    dQ[1] = -q1 * v∇Φ1
    dQ[2] = v1
    dQ[3] = q1 * a12 * (v1 * Π1 + Ψ1) / m1
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

end # end of module
