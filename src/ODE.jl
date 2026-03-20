module ODE

using ..Reconstructions
using ..ParticleMotion
using ..ScalarField
using ..Equations
using StaticArrays
using Polyester
using DataInterpolations

function muscl!(dQ, Q, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
    slope_limiter = params.slope_limiter

    h = params.h
    N = params.N
    # x = params.x
    # L0 = params.L0
    # half_h = 0.5h
    _h = one(h) / h
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)
    @inbounds begin
        q = @view Q[:, 1:2]
        # @batch for i in 1:N
        for i in 1:N
            Qi = SVector{2}(q[i, :])
            Qim1 = SVector{2}(q[i - 1, :])
            Qip1 = SVector{2}(q[i + 1, :])
            Qim2 = SVector{2}(q[i - 2, :])
            Qip2 = SVector{2}(q[i + 2, :])
            ## i cell
            Δminus = Qi - Qim1
            Δplus = Qip1 - Qi

            # calculate characteristic jumps
            aminus = L * Δminus
            aplus = L * Δplus
            # limit values
            si = slope_limiter.(aminus, aplus)
            # project back
            σi = R * si
            # reconstruct at faces
            Q_i_plus = Qi + 0.5σi
            Q_i_minus = Qi - 0.5σi

            ## i-1 cell
            Δminus = Qim1 - Qim2
            Δplus = Qi - Qim1
            # calculate characteristic jumps
            aminus = L * Δminus
            aplus = L * Δplus
            # limit values
            si = slope_limiter.(aminus, aplus)
            # project back
            σi = R * si
            # reconstruct at faces
            Q_im1_plus = Qim1 + 0.5σi
            # Q_im1_minus = q[i-1, :] - 0.5σi

            ## i+1 cell
            Δminus = Qip1 - Qi
            Δplus = Qip2 - Qip1
            # calculate characteristic jumps
            aminus = L * Δminus
            aplus = L * Δplus
            # limit values
            si = slope_limiter.(aminus, aplus)
            # project back
            σi = R * si
            # reconstruct at faces
            # Q_ip1_plus = q[i+1, :] + 0.5σi
            Q_ip1_minus = Qip1 - 0.5σi

            ## Fl: flux on the left face, accounts for incoming right going flux
            ## Fr: flux on the right face, accounts for incoming left going flux
            Fl = numerical_flux(Q_im1_plus, Q_i_minus, equation)
            Fr = numerical_flux(Q_i_plus, Q_ip1_minus, equation)

            dQ[i, 1:2] .= -(Fr - Fl) * _h
            # dQ[i, 3] = Q[i, 1]
        end
    end
    dQ[:, 3] .= @view Q[:, 1]
    return nothing
end

@inline function get_q_periodic(Q, idx, N)
    # mod(idx - 1, N) maps 1..N to 0..N-1
    # + 1 shifts it back to 1..N
    i = mod(idx - 1, N) + 1
    return SVector{2}(Q[i, 1], Q[i, 2])
end

# Helper: Purely stack-allocated reconstruction
@inline function reconstruct_cell(Q, idx, N, L, R, slope_limiter)
    # 1. Load data directly into registers (avoids SubArray views)
    Qi = get_q_periodic(Q, idx, N)
    Qim1 = get_q_periodic(Q, idx - 1, N)
    Qip1 = get_q_periodic(Q, idx + 1, N)

    # 2. Vector Math
    Δminus = Qi - Qim1
    Δplus = Qip1 - Qi

    # 3. Characteristic jumps 
    # specific dispatch ensures L is SMatrix, so this is allocation-free
    aminus = L * Δminus
    aplus = L * Δplus

    # 4. Limiters
    si = slope_limiter.(aminus, aplus)

    # 5. Project back
    σi = R * si

    # 6. Reconstruct faces (Left Face, Right Face)
    q_minus = Qi - 0.5 * σi
    q_plus = Qi + 0.5 * σi

    return q_minus, q_plus
end

function muscl_optimized!(dQ, Q, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
    slope_limiter = params.slope_limiter

    h = params.h
    N = params.N
    _h = one(h) / h

    # These return @SMatrix directly based on our definition.
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)

    @inbounds begin
        # --- Startup: Interface 1/2 ---
        # We need the Right face of ghost cell '0' and Left face of cell '1'
        _, q_0_plus = reconstruct_cell(Q, 0, N, L, R, slope_limiter)
        q_1_minus, q_1_plus = reconstruct_cell(Q, 1, N, L, R, slope_limiter)

        # Compute Flux at first interface
        Fl = numerical_flux(q_0_plus, q_1_minus, equation)

        # Save the "Right face of current cell" state for the loop
        q_i_plus = q_1_plus

        # --- Main Loop: Sliding Window ---
        for i in 1:N
            # We already have Fl (Flux on Left). 
            # We need Fr (Flux on Right).

            # Reconstruct the NEXT cell (i+1) to get its Left Face
            q_ip1_minus, q_ip1_plus = reconstruct_cell(Q, i + 1, N, L, R, slope_limiter)

            # Compute Flux
            Fr = numerical_flux(q_i_plus, q_ip1_minus, equation)

            # Update dQ
            # Manual indexing avoids allocating a view for dQ[i, 1:2]
            diff_flux = -(Fr - Fl) * _h
            dQ[i, 1] = diff_flux[1]
            dQ[i, 2] = diff_flux[2]

            # Slide window: Current Right Flux becomes Next Left Flux
            Fl = Fr
            # Current i+1 reconstruction becomes i for next step
            q_i_plus = q_ip1_plus
        end
        dQ[:, 3] .= @view Q[:, 1]
    end
    return nothing
end

# Pure scalar WENO-5 reconstruction
# Reconstructs value at right face (i+1/2) using stencil (i-2, i-1, i, i+1, i+2)
@inline function weno5_scalar(v1, v2, v3, v4, v5)
    epsilon = 1e-6 # To avoid division by zero

    # 1. Smoothness Indicators (Jiang & Shu)
    # beta0 corresponds to stencil {i-2, i-1, i}
    beta0 = 13 / 12 * (v1 - 2 * v2 + v3)^2 + 1 / 4 * (v1 - 4 * v2 + 3 * v3)^2
    # beta1 corresponds to stencil {i-1, i, i+1}
    beta1 = 13 / 12 * (v2 - 2 * v3 + v4)^2 + 1 / 4 * (v2 - v4)^2
    # beta2 corresponds to stencil {i, i+1, i+2}
    beta2 = 13 / 12 * (v3 - 2 * v4 + v5)^2 + 1 / 4 * (3 * v3 - 4 * v4 + v5)^2

    # 2. Linear Weights
    d0 = 0.1
    d1 = 0.6
    d2 = 0.3

    # 3. Alpha Weights (Non-linear)
    alpha0 = d0 / (epsilon + beta0)^2
    alpha1 = d1 / (epsilon + beta1)^2
    alpha2 = d2 / (epsilon + beta2)^2
    sum_alpha = alpha0 + alpha1 + alpha2

    # 4. Normalized Weights
    w0 = alpha0 / sum_alpha
    w1 = alpha1 / sum_alpha
    w2 = alpha2 / sum_alpha

    # 5. Stencil Reconstructions (q_k)
    q0 = 1 / 3 * v1 - 7 / 6 * v2 + 11 / 6 * v3
    q1 = -1 / 6 * v2 + 5 / 6 * v3 + 1 / 3 * v4
    q2 = 1 / 3 * v3 + 5 / 6 * v4 - 1 / 6 * v5

    # 6. Final Combination
    return w0 * q0 + w1 * q1 + w2 * q2
end

# Characteristic-wise WENO5 Reconstruction
@inline function reconstruct_cell_weno5(Q, idx, N, L, R)
    # 1. Load Stencil (5 points: i-2 to i+2)
    Q_mm = get_q_periodic(Q, idx - 2, N) # i-2
    Q_m = get_q_periodic(Q, idx - 1, N) # i-1
    Q_c = get_q_periodic(Q, idx, N) # i (Center)
    Q_p = get_q_periodic(Q, idx + 1, N) # i+1
    Q_pp = get_q_periodic(Q, idx + 2, N) # i+2

    # 2. Project to Characteristic Variables (W = L * Q)
    #    L is SMatrix, so this is non-allocating.
    W_mm = L * Q_mm
    W_m = L * Q_m
    W_c = L * Q_c
    W_p = L * Q_p
    W_pp = L * Q_pp

    # 3. Perform WENO5 on each component independently
    #    We assume a 2-component system (Wave Equation)

    # --- Component 1 ---
    w1_plus = weno5_scalar(W_mm[1], W_m[1], W_c[1], W_p[1], W_pp[1])
    # To get Left Face (i-1/2), we reverse the stencil inputs
    w1_minus = weno5_scalar(W_pp[1], W_p[1], W_c[1], W_m[1], W_mm[1])

    # --- Component 2 ---
    w2_plus = weno5_scalar(W_mm[2], W_m[2], W_c[2], W_p[2], W_pp[2])
    w2_minus = weno5_scalar(W_pp[2], W_p[2], W_c[2], W_m[2], W_mm[2])

    # 4. Assemble Characteristic Vectors
    W_plus = SVector{2}(w1_plus, w2_plus)
    W_minus = SVector{2}(w1_minus, w2_minus)

    # 5. Project back to Physical Space (Q = R * W)
    q_plus = R * W_plus
    q_minus = R * W_minus

    return q_minus, q_plus
end

function muscl_weno5!(dQ, Q, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
    # slope_limiter is NOT used for WENO (it has built-in limiting)

    h = params.h
    N = params.N
    _h = one(h) / h

    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)

    @inbounds begin
        # --- Startup Phase ---
        # We need fluxes at the first few interfaces. 
        # Since the stencil is wider, we just rely on the periodic helper 
        # to handle the wrapping logic at i=0, i=1.

        # Get q_plus (Right Face) of cell 0 (Left neighbor of cell 1)
        _, q_0_plus = reconstruct_cell_weno5(Q, 0, N, L, R)

        # Get q_minus (Left Face) of cell 1
        q_1_minus, q_1_plus = reconstruct_cell_weno5(Q, 1, N, L, R)

        # Compute first Flux at interface 1/2
        Fl = numerical_flux(q_0_plus, q_1_minus, equation)

        # Save state
        q_i_plus = q_1_plus

        # --- Main Loop ---
        for i in 1:N
            # 1. Reconstruct neighbor (i+1)
            q_ip1_minus, q_ip1_plus = reconstruct_cell_weno5(Q, i + 1, N, L, R)

            # 2. Compute Flux (Flux at i + 1/2)
            Fr = numerical_flux(q_i_plus, q_ip1_minus, equation)

            # 3. Update dQ
            diff_flux = -(Fr - Fl) * _h
            dQ[i, 1] = diff_flux[1]
            dQ[i, 2] = diff_flux[2]

            # 4. Slide Window
            Fl = Fr
            q_i_plus = q_ip1_plus
        end
        dQ[:, 3] .= @view Q[:, 1]
    end
    return nothing
end

# Helper: purely scalar math, updates dQ in place
@inline function add_singular_sources!(dQ, idx, x_i, half_h, _h, z, v, q, a)
    # Pre-calculate faces
    left_face = x_i - half_h
    right_face = x_i + half_h

    one_v = one(v)
    inv_L = inv((one_v + v)^2)
    inv_R = inv((one_v - v)^2)

    # -a * ∂ᵥΠs and -a * ∂ₓ∂ᵥΦs = -a * ∂ᵥΨs are piecewise constant.
    pi_val_L = -0.5 * q * a * inv_L
    pi_val_R = 0.5 * q * a * inv_R
    psi_val_L = -0.5 * q * a * inv_L
    psi_val_R = -0.5 * q * a * inv_R

    # -a * ∂ᵥΦs is piecewise linear, so we keep exact cell averages.
    phi_slope_L = 0.5 * q * a * inv_L
    phi_slope_R = -0.5 * q * a * inv_R
    if z >= right_face
        # CASE A: Particle is to the right. 
        # Entire cell is in the Left Regime (x < z).
        @inbounds begin
            dQ[idx, 1] += pi_val_L
            dQ[idx, 2] += psi_val_L
            dQ[idx, 3] += phi_slope_L * (z - x_i)
        end

    elseif z <= left_face
        # CASE B: Particle is to the left.
        # Entire cell is in the Right Regime (x > z).
        @inbounds begin
            dQ[idx, 1] += pi_val_R
            dQ[idx, 2] += psi_val_R
            dQ[idx, 3] += phi_slope_R * (x_i - z)
        end

    else
        # CASE C: Particle is Inside (Split Domain).
        # We sum the contributions from both sides.

        # Fraction of cell to the left of z: (z - left_face) / h
        frac_L = (z - left_face) * _h

        # Fraction of cell to the right of z
        frac_R = one(frac_L) - frac_L

        dx_L = z - left_face
        dx_R = right_face - z

        @inbounds begin
            dQ[idx, 1] += frac_L * pi_val_L + frac_R * pi_val_R
            dQ[idx, 2] += frac_L * psi_val_L + frac_R * psi_val_R
            dQ[idx, 3] += 0.5 * _h * (phi_slope_L * dx_L^2 + phi_slope_R * dx_R^2)
        end
    end
    return nothing
end

function field_rhs_forced!(dQ, Q, params, t)
    x = params.x
    N = params.N
    q = params.q
    h = params.h

    z0 = params.z0
    A = params.A
    ω = params.ω

    muscl_optimized!(dQ, Q, params, t)

    z, v, a = ParticleMotion.oscillator(t, z0, A, ω)

    dtΠ = @view dQ[:, 1]
    dtΦ = @view dQ[:, 3]
    _h = one(h) / h
    half_h = 0.5 * h
    @inbounds for i in 1:N
        add_singular_sources!(dQ, i, x[i], half_h, _h, z, v, q, a)
    end
    return nothing
end

function field_rhs!(dQ, Q, params, t)
    x = params.x
    N = params.N
    q = params.q
    h = params.h

    Q_field = Q.x[1]
    Q_particle = Q.x[2]
    dQ_field = dQ.x[1]
    dQ_particle = dQ.x[2]

    a = dQ_particle[3]
    z = Q_particle[2]
    v = Q_particle[3] ### ASK ERIK dQ or Q ?

    muscl_optimized!(dQ_field, Q_field, params, t)

    # dtΠ = @view dQ_field[:, 1]
    # dtΦ = @view dQ_field[:, 3]

    _h = one(h) / h
    half_h = 0.5 * h
    @inbounds for i in 1:N
        add_singular_sources!(dQ_field, i, x[i], half_h, _h, z, v, q, a)
    end
    return nothing
end

function radiation_reaction_rhs!(dQ, Q, params, t)
    particle_rhs!(dQ, Q, params, t)
    field_rhs!(dQ, Q, params, t)
    return nothing
end

function radiation_reaction_rhs_relativistic!(dQ, Q, params, t)
    particle_rhs_relativistic!(dQ, Q, params, t)
    field_rhs!(dQ, Q, params, t)
    return nothing
end

@inline function find_cell_index(x, dx, grid_origin)
    return floor(Int64, (x - grid_origin) / dx) + 1
end

function particle_rhs!(dQ, Q, params, t)
    Q_field = Q.x[1]
    Q_particle = Q.x[2]
    dQ_particle = dQ.x[2]

    q = params.q
    h = params.h
    xs = params.x
    dxV = params.dxPotential
    sf = params.sf
    m0 = params.m0

    interpolation_method = params.interpolation_method

    # m = Q_particle[1]
    z = Q_particle[2]
    v = Q_particle[3]

    vel = 0.999999
    abs(v) <= vel || throw("v=$v > $(vel), at t=$(t), a = $(dQ_particle[3])")

    Ψr = @view Q_field[:, 2]
    Φr = @view Q_field[:, 3]
    interpolator_Ψ = interpolation_method(Ψr, xs;
                                          extrapolation=ExtrapolationType.Extension)
    interpolator_Φ = interpolation_method(Φr, xs;
                                          extrapolation=ExtrapolationType.Extension)
    Φrz = interpolator_Φ(z)
    Ψrz = interpolator_Ψ(z)

    @inbounds dQ_particle[2] = v
    @inbounds dQ_particle[3] = -q * (Ψrz + dxV(z)) / (m0 + q * Φrz)
    return nothing
end

function particle_rhs_relativistic!(dQ, Q, params, t)
    Q_field = Q.x[1]
    Q_particle = Q.x[2]
    dQ_particle = dQ.x[2]

    q = params.q
    h = params.h
    xs = params.x
    dxV = params.dxPotential
    # V = params.Potential
    # sf = params.sf
    m0 = params.m0

    interpolation_method = params.interpolation_method

    # m = Q_particle[1]
    z = Q_particle[2]
    v = Q_particle[3]

    Πr = @view Q_field[:, 1]
    Ψr = @view Q_field[:, 2]
    Φr = @view Q_field[:, 3]
    interpolator_Ψ = interpolation_method(Ψr, xs;
                                          extrapolation=ExtrapolationType.Extension)
    interpolator_Π = interpolation_method(Πr, xs;
                                          extrapolation=ExtrapolationType.Extension)
    interpolator_Φ = interpolation_method(Φr, xs;
                                          extrapolation=ExtrapolationType.Extension)
    Φrz = interpolator_Φ(z)
    Πrz = interpolator_Π(z)
    Ψrz = interpolator_Ψ(z)

    vel = 0.999999
    abs(v) <= vel ||
        throw("v=$v > $(vel), at t=$(t), a = $(dQ_particle[3]), den= $(m0 - q * Φrz)")

    @inbounds dQ_particle[1] = q * (v * Ψrz + Πrz)
    @inbounds dQ_particle[2] = v

    num = -q * dxV(z) * (1 - v * v)^1.5 + q * (1 - v * v) * (v * Πrz + Ψrz)
    den = m0 - q * Φrz
    @inbounds dQ_particle[3] = num / den
    return nothing
end

# function particle_rhs!(dQ, Q, params, t)
#     Q_field = Q.x[1]
#     Q_particle = Q.x[2]
#     dQ_particle = dQ.x[2]

#     Π = @view Q_field[:, 1]
#     Ψ = @view Q_field[:, 2]
#     x = params.x
#     interpolation_method = params.interpolation_method
#     q = params.q

#     m1 = Q_particle[1]
#     x1 = Q_particle[2]
#     v1 = Q_particle[3]
#     abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ_particle[3])")

#     interpolator_Ψ = interpolation_method(Ψ, x; extrapolation=ExtrapolationType.Extension)
#     interpolator_Π = interpolation_method(Π, x; extrapolation=ExtrapolationType.Extension)
#     Π1 = interpolator_Π(x1)
#     Ψ1 = interpolator_Ψ(x1)

#     v12 = v1 * v1
#     a12 = one(v1) - v12

#     dQ_particle[1] = -q * (Π1 + v1 * Ψ1)
#     dQ_particle[2] = v1
#     dQ_particle[3] = q * a12 * (v1 * Π1 + Ψ1) / m1
#     return nothing
# end

# function particle_rhs_functions!(dQ, Q, params, t)
#     Π = params.pifield
#     Ψ = params.psifield
#     q1 = params.q1

#     @inbounds m1 = Q[1]
#     @inbounds x1 = Q[2]
#     @inbounds v1 = Q[3]
#     abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ[3])")

#     Π1 = Π(t, x1)
#     Ψ1 = Ψ(t, x1)

#     v12 = v1 * v1
#     a12 = one(v1) - v12

#     @inbounds dQ[1] = -q1 * (Π1 + v1 * Ψ1)
#     @inbounds dQ[2] = v1
#     @inbounds dQ[3] = q1 * a12 * (v1 * Π1 + Ψ1) / m1
#     return nothing
# end

# function particle_rhs_interpolation!(dQ, Q, params, t)
#     Π = params.pifield
#     Ψ = params.psifield
#     q1 = params.q1
#     x = params.x
#     interpolation_method = params.interpolation_method
#     @inbounds m1 = Q[1]
#     @inbounds x1 = Q[2]
#     @inbounds v1 = Q[3]
#     abs(v1) <= 1 || throw("v1=$v1 > 1, at t=$(t), a1 = $(dQ[3])")
#     Πvals = Π.(t, x)
#     Ψvals = Ψ.(t, x)

#     interpolator_Ψ = interpolation_method(Ψvals, x;
#                                           extrapolation=ExtrapolationType.Extension)
#     interpolator_Π = interpolation_method(Πvals, x;
#                                           extrapolation=ExtrapolationType.Extension)
#     Π1 = interpolator_Π(x1)
#     Ψ1 = interpolator_Ψ(x1)

#     v12 = v1 * v1
#     a12 = one(v1) - v12

#     v∇Φ1 = Π1 + v1 * Ψ1

#     @inbounds dQ[1] = -q1 * v∇Φ1
#     @inbounds dQ[2] = v1
#     @inbounds dQ[3] = q1 * a12 * (v1 * Π1 + Ψ1) / m1
#     return nothing
# end

# @inline function get_stencil2(Q, i, N)
#     @inbounds begin
#         if i == 0
#             qi = Q[N, :]
#             ql = Q[N - 1, :]
#             qr = Q[1, :]
#         elseif i == 1
#             qi = Q[1, :]
#             ql = Q[N, :]
#             qr = Q[2, :]
#         elseif i == N
#             qi = Q[N, :]
#             ql = Q[N - 1, :]
#             qr = Q[1, :]
#         elseif i == N + 1
#             qi = Q[1, :]
#             ql = Q[N, :]
#             qr = Q[2, :]
#         else
#             qi = Q[i, :]
#             ql = Q[i - 1, :]
#             qr = Q[i + 1, :]
#         end
#     end
#     return qi, ql, qr
# end

# @inline function get_stencil(Q, i, N)
#     @inbounds begin
#         if i == 0
#             qi = Q[N]
#             ql = Q[N - 1]
#             qr = Q[1]
#         elseif i == 1
#             qi = Q[1]
#             ql = Q[N]
#             qr = Q[2]
#         elseif i == N
#             qi = Q[N]
#             ql = Q[N - 1]
#             qr = Q[1]
#         elseif i == N + 1
#             qi = Q[1]
#             ql = Q[N]
#             qr = Q[2]
#         else
#             qi = Q[i]
#             ql = Q[i - 1]
#             qr = Q[i + 1]
#         end
#     end
#     return qi, ql, qr
# end

# function trapezoidal(y, h)
#     n = length(y) - 1
#     integral = zero(Float64)
#     @inbounds for i in 1:n
#         integral += y[i]
#     end
#     @inbounds integral += 0.5 * (y[1] + y[n + 1])  # Adjust endpoints
#     integral *= h
#     return integral
# end

# # LaxWendroff schemes
# function LxW!(dQ, Q, params, t)
#     equation = params.equation
#     # numerical_flux = params.numerical_flux
#     flux_limiter = params.flux_limiter
#     h = params.h
#     N = params.N
#     x = params.x
#     cfl = params.cfl
#     c = equation.velocity
#     dt = params.dt

#     half_h = h * 0.5
#     _h = 1.0 / h
#     r1 = @SVector [c, one(c)]
#     r2 = @SVector [-c, one(c)]
#     λp = c
#     λm = -c
#     _2c = 1.0 / (2c)
#     @fastmath @inbounds begin
#         @batch for i in 1:(N + 1)
#             ## Left Face
#             ΔΠl = (Q[i, 1] - Q[i - 1, 1]) * _2c
#             ΔΨl = (Q[i, 2] - Q[i - 1, 2]) * 0.5
#             a1l = ΔΠl + ΔΨl
#             a2l = -ΔΠl + ΔΨl
#             # W1l = a1l * r1 # left moving at left face
#             W2l = a2l * r2 # right moving at left face

#             ## Right Face
#             ΔΠr = (Q[i + 1, 1] - Q[i, 1]) * _2c
#             ΔΨr = (Q[i + 1, 2] - Q[i, 2]) * 0.5
#             a1r = ΔΠr + ΔΨr
#             a2r = -ΔΠr + ΔΨr
#             W1r = a1r * r1 # left moving at right face
#             # W2r = a2r * r2 # right moving at right face

#             ### Higher resolution stuff - Limiters
#             ## Left Face
#             # for λ<0
#             I = i + 1
#             ΔΠlp = (Q[I, 1] - Q[I - 1, 1]) * _2c
#             ΔΨlp = (Q[I, 2] - Q[I - 1, 2]) * 0.5
#             a1lp = ΔΠlp + ΔΨlp
#             a1lp == a1l ? θ1l = 1.0 : θ1l = a1lp / a1l

#             a1l_tilde = flux_limiter(θ1l) * a1l
#             W1l_tilde = a1l_tilde * r1
#             # for λ>0
#             I = i - 1
#             ΔΠlm = (Q[I, 1] - Q[I - 1, 1]) * _2c
#             ΔΨlm = (Q[I, 2] - Q[I - 1, 2]) * 0.5
#             a2lm = -ΔΠlm + ΔΨlm
#             a2lm == a2l ? θ2l = 1.0 : θ2l = a2lm / a2l

#             a2l_tilde = flux_limiter(θ2l) * a2l
#             W2l_tilde = a2l_tilde * r2

#             ## Right Face
#             # for λ<0
#             I = i + 1
#             ΔΠrp = (Q[I + 1, 1] - Q[I, 1]) * _2c
#             ΔΨrp = (Q[I + 1, 2] - Q[I, 2]) * 0.5
#             a1rp = ΔΠrp + ΔΨrp

#             a1rp == a1r ? θ1r = 1.0 : θ1r = a1rp / a1r

#             a1r_tilde = flux_limiter(θ1r) * a1r
#             W1r_tilde = a1r_tilde * r1
#             # for λ>0
#             I = i - 1
#             ΔΠrm = (Q[I + 1, 1] - Q[I, 1]) * _2c
#             ΔΨrm = (Q[I + 1, 2] - Q[I, 2]) * 0.5
#             a2rm = -ΔΠrm + ΔΨrm
#             a2rm == a2r ? θ2r = 1.0 : θ2r = a2rm / a2r
#             a2r_tilde = flux_limiter(θ2r) * a2r
#             W2r_tilde = a2r_tilde * r2

#             tmp = λp * (1.0 - cfl)

#             Fl_tilde = tmp * (W1l_tilde + W2l_tilde)
#             Fr_tilde = tmp * (W1r_tilde + W2r_tilde)

#             dQ[i, :] .= -(λm * W1r + λp * W2l) * _h - 0.5(Fr_tilde - Fl_tilde) * _h
#         end
#     end
#     return nothing
# end

end # end of module
