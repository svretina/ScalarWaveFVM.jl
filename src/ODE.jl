module ODE

using ..Reconstructions
using StaticArrays

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

    @fastmath @inbounds for i in 1:(N + 1)
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
    return nothing
end

function get_stencil(Q, i, N)
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

    Π = @view Q[:, 1]
    Ψ = @view Q[:, 2]

    @fastmath @inbounds for i in 1:N
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
    return nothing
end

function rhs!(dQ, Q, params, t)
    x1 = params.x1
    x2 = params.x2
    x = params.x
    L = params.L
    direction1 = params.direction1
    direction2 = params.direction2
    q1 = params.q1
    q2 = params.q2
    N = params.N
    h = params.h

    muscl!(dQ, Q, params, t)
    xp1, vp2, ap2 = ParticleMotion.oscillator(t, x1, L, direction1)
    xp1, vp2, ap2 = ParticleMotion.oscillator(t, x2, L, direction2)
    s1(x) = ap1 * ScalarField.∂vΠ(x, q1, xp1, vp1)
    s2(x) = ap2 * ScalarField.∂vΠ(x, q2, xp2, vp2)
    dtΠ = @view dQ[:, 1]

    @batch for i in 1:N
        dtΠ[i] -= s1(x[i]) + s2(x[i])
    end
    # the source term should be just averaged over the cell
    # meaning that i need to calculate
    # an integral
    #
    # as 0th order approximation i will just sample
    # the source term at the cell center
    #
    # then as a 1st order approximation
    # I can evaluate it at the faces and take the average
    # of it
    #
    # the exact thing to do is to
    # take the integral of it
    #
    # If this does not work, consider implementing
    # the Fractional Step method of LeVeque Chapter 17.1 pg 377

end

end # end of module
