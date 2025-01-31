module ScalarWaveFVM

include("Equations.jl")

include("Limiters.jl")
include("Reconstructions.jl")
include("NumericalFluxes.jl")

using ..Equations
using ..Reconstructions
using ..NumericalFluxes

using GridFunctions
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqLowOrderRK
using StaticArrays
using CircularArrays
using ForwardDiff
using LinearAlgebra

@inline function Gaussian1D(t::Real, x::Real, A::Real, σ::Real)
    L = 100
    xc = t
    xc = mod(xc + L, 2L) - L

    x_shifted1 = x - xc + 2L
    g1 = A * exp(-x_shifted1^2 / (2 * σ^2))
    x_shiftedm1 = x - xc - 2L
    gm1 = A * exp(-x_shiftedm1^2 / (2 * σ^2))
    x_shifted0 = x - xc
    g0 = A * exp(-x_shifted0^2 / (2 * σ^2))
    return g1 + gm1 + g0
end

@inline function dtGaussian1D(t::Real, x::Real, A::Real, σ::Real)
    return ForwardDiff.derivative(t1 -> Gaussian1D(t1, x, A, σ), t)
end

@inline function dxGaussian1D(t::Real, x::Real, A::Real, σ::Real)
    return ForwardDiff.derivative(x1 -> Gaussian1D(t, x1, A, σ), x)
end

function run(L, N, tf, cfl)
    grid = UniformGrid(Float64[-L, L], N)
    x = coords(grid)
    Π = CircularVector(zeros(Float64, N + 1))
    Ψ = CircularVector(zeros(Float64, N + 1))

    # Ψ[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= -1.0 # -1 for right moving, +1 for left moving
    # Π[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= 1.0

    σ = 20.0
    A = 10.0
    Ψ .= dxGaussian1D.(0.0, x, A, σ)
    Π .= dtGaussian1D.(0.0, x, A, σ)

    statevector = hcat(Π, Ψ)

    tspan = (0.0, tf)
    dx = spacing(grid)

    # alg = SSPRK54()
    alg = SSPRK22()
    # alg = Euler()

    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    dt = cfl * dx / c
    @show dx, dt
    numerical_flux = NumericalFluxes.flux_godunov
    # reconstruct_func = Reconstructions.piecewise_linear
    # slope_func = Reconstructions.constant_slope
    flux_limiter = Limiters.Lax_Wendroff
    params = (equation=equation, numerical_flux=numerical_flux,
              #reconstruct_func=reconstruct_func,
              #slope_func=slope_func,
              flux_limiter=flux_limiter,
              h=dx, N=N, x=GridFunctions.coords(grid), cfl=cfl,
              dt=dt)

    ode = ODEProblem{true}(rhs2!, statevector, tspan, params)
    # return statevector, params
    sol = solve(ode, alg; adaptive=false, dt=dt)
    return GridFunctions.coords(grid), sol
    # du = similar(statevector)
    # rhs!(du, statevector, params, 0.1)

    # return du
end

function rhs2!(dQ, Q, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
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
    #loop_over_1:(N
    @inbounds for i in 1:(N + 1)
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

        ## 4 for cfl=1/4, 2 for cfl=1/2
        ## WHY?! else I get phase error!
        ## for those values, tmp = 0 and
        ## i fall back to the upwind method,
        ## thus the 1st order convergence.
        ## although I get only 1st order conv
        tmp = λp * (1.0 - cfl)

        Fl_tilde = tmp * (W1l_tilde + W2l_tilde)
        Fr_tilde = tmp * (W1r_tilde + W2r_tilde)

        dQ[i, :] .= -(λm * W1r + λp * W2l) * _h - 0.5 * (Fr_tilde - Fl_tilde) * _h
    end
    return nothing
end

function rhs!(dQ, Q, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
    flux_limiter = params.flux_limiter
    h = params.h
    N = params.N
    x = params.x
    cfl = params.cfl
    c = equation.velocity

    half_h = h * 0.5
    _h = 1.0 / h
    @inbounds for i in 1:N
        # left_face = x[i] - half_h
        # right_face = x[i] + half_h

        ## no need for the reconstructions with the limiters!
        ## rewrite!
        # construct reconstruction functions based on neighboring cells
        # for Π

        # Riπ = reconstruct_func(x[i], Q[i, 1], Q[i - 1, 1], Q[i + 1, 1], params)
        # Rlπ = reconstruct_func(x[i], Q[i - 1, 1], Q[i - 2, 1], Q[i, 1], params)
        # Rrπ = reconstruct_func(x[i], Q[i + 1, 1], Q[i, 1], Q[i + 2, 1], params)

        # # construct reconstruction functions based on neighboring cells
        # # for Ψ
        # Riψ = reconstruct_func(x[i], Q[i, 2], Q[i - 1, 2], Q[i + 1, 2], params)
        # Rlψ = reconstruct_func(x[i], Q[i - 1, 2], Q[i - 2, 2], Q[i, 2], params)
        # Rrψ = reconstruct_func(x[i], Q[i + 1, 2], Q[i, 2], Q[i + 2, 2], params)

        # # reconstruct Π values at faces
        # πlr = Riπ(left_face)
        # πll = Rlπ(left_face)

        # πrr = Rrπ(right_face)
        # πrl = Riπ(right_face)

        # # reconstruct Ψ values at faces
        # ψlr = Riψ(left_face)
        # ψll = Rlψ(left_face)

        # ψrr = Rrψ(right_face)
        # ψrl = Riψ(right_face)

        # ### Construct differences at faces
        # ###
        # Δql = @SVector [πlr - πll, ψlr - ψll]
        # Δqr = @SVector [πrr - πrl, ψrr - ψrl]
        Δql = @SVector [Q[i, 1] - Q[i - 1, 1], Q[i, 2] - Q[i - 1, 2]]
        Δqr = @SVector [Q[i + 1, 1] - Q[i, 1], Q[i + 1, 2] - Q[i, 2]]
        ## +1  right direction
        ## -1 left direction
        ## Fl: flux on the left face, accounts for incoming right going flux
        ## Fr: flux on the right face, accounts for incoming left going flux
        if flux_limiter == Limiters.upwind
            Fl_tilde = @SVector [0.0, 0.0]
            Fr_tilde = @SVector [0.0, 0.0]
        else
            θl = Limiters.θ(i, Q, equation, +1)
            θr = Limiters.θ(i, Q, equation, -1)

            ϕl = flux_limiter.(θl)
            ϕr = flux_limiter.(θr)

            A = Equations.Aabs(equation)

            Fl_tilde = 0.5 * A * (I - cfl * A) * Δql .* ϕl
            Fr_tilde = 0.5 * A * (I - cfl * A) * Δqr .* ϕr
        end
        # @show ϕl, Fl_tilde

        Fl = numerical_flux(Q[i - 1, :], Q[i, :], equation; cfl=cfl) + Fl_tilde
        Fr = numerical_flux(Q[i, :], Q[i + 1, :], equation; cfl=cfl) + Fr_tilde

        dQ[i, :] .= -(Fr - Fl) * _h
    end
    return nothing
end

end
