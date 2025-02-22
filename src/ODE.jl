module ODE

using ..Reconstructions
using ..ParticleMotion
using ..ScalarField
using StaticArrays
using Polyester

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
    xp1, vp1, ap1 = ParticleMotion.oscillator(t, x1, L, direction1)
    xp2, vp2, ap2 = ParticleMotion.oscillator(t, x2, L, direction2)
    @inline s1(x) = ap1 * ScalarField.∂vΠ(x, q1, xp1, vp1)
    @inline s2(x) = ap2 * ScalarField.∂vΠ(x, q2, xp2, vp2)
    dtΠ = @view dQ[:, 1]

    for i in 1:N
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

## potential field is provided as a function
function particle_motion!(du, U, params, t)
    fΠ = params.pifield
    fΨ = params.psifield

    q1 = params.q1
    q2 = params.q2

    m1 = U[1]
    x1 = U[2]
    v1 = U[3]
    @assert v1 <= 1
    m2 = U[4]
    x2 = U[5]
    v2 = U[6]
    @assert v2 <= 1

    Π1 = fΠ(t, x1)
    Ψ1 = fΨ(t, x1)
    Π2 = fΠ(t, x2)
    Ψ2 = fΨ(t, x2)

    v12 = v1 * v1
    a1 = sqrt(one(v1) - v12) # reciprocal
    a12 = a1 * a1
    #a14 = a12 * a12
    @assert a1 > 0
    γ1 = one(a1) / a1
    v22 = v2 * v2
    a2 = sqrt(one(v2) - v22) # reciprocal
    a22 = a2 * a2
    # a24 = a22 * a22
    @assert a2 > 0
    γ2 = one(a2) / a2 # sqrt(β² + 1)

    # temporary variables
    # check mathematica file: RelativisticMotion.nb
    # uᵗΠ1 = uᵗ1 * Π1
    # u∇Φ1 = uᵗΠ1 + u1 * Ψ1

    # uᵗΠ2 = uᵗ2 * Π2
    # u∇Φ2 = uᵗΠ2 + u2 * Ψ2

    # # particle 1
    # du[1] = -q1 * u∇Φ1 ## mass evolution
    # du[2] = u1       ## x-coordinate
    # du[3] = (a14 * q1 / m1) * (u1 * u∇Φ1 / a12 + Ψ1) ## x-velocity
    # # particle 2
    # du[4] = -q2 * u∇Φ2 ## mass evolution
    # du[5] = u2       ## x-coordinate
    # du[6] = (a24 * q2 / m2) * (u2 * u∇Φ2 / a22 + Ψ2) ## x-velocity

    v∇Φ1 = -Π1 + v1 * Ψ1
    v∇Φ2 = -Π2 + v2 * Ψ2
    # mass
    du[1] = -q1 * v∇Φ1
    du[4] = -q2 * v∇Φ2
    # position
    du[2] = v1
    du[5] = v2
    # velocity
    γ12 = γ1 * γ1
    γ22 = γ2 * γ2

    du[3] = (q1 / m1) * (-γ1 * Π1 / (one(γ1) + v12 * γ12) + a12 * Ψ1)
    du[6] = (q2 / m2) * (-γ2 * Π2 / (one(γ2) + v22 * γ22) + a22 * Ψ2)

    return nothing
end

end # end of module
