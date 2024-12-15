module ScalarWaveFVM

include("Reconstructions.jl")
include("NumericalFluxes.jl")

using ..Reconstructions
using ..NumericalFluxes
using OrdinaryDiffEqSSPRK

struct LinearScalarAdvectionEquation1D{RealT<:Real}
    advection_velocity::RealT
end

# Calculate 1D flux in for a single point
@inline function flux(u, equation::LinearScalarAdvectionEquation1D)
    a = equation.advection_velocity
    return a * u
end

function run(L, N, tf, cfl)
    xs = LinRange(-L, L, N + 1)
    statevector = zeros(N + 1)
    statevector[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= 1.0
    tspan = (0.0, tf)
    @show tf
    dx = (xs[2] - xs[1])
    dt = cfl * dx
    alg = SSPRK54()

    equation = LinearScalarAdvectionEquation1D(1.0)
    numerical_flux = NumericalFluxes.flux_godunov
    reconstruct_func = Reconstructions.piecewise_linear
    slope_func = Reconstructions.constant_slope
    params = (equation=equation, numerical_flux=numerical_flux,
              reconstruct_func=reconstruct_func, slope_func=slope_func,
              h=dx, N=N, x=xs)
    @show params
    ode = ODEProblem{true}(rhs!, statevector, tspan, params)
    sol = solve(ode, alg; adaptive=false, dt=dt)
    return xs, sol
end

function rhs!(du, U, params, t)
    equation = params.equation
    numerical_flux = params.numerical_flux
    reconstruct_func = params.reconstruct_func
    h = params.h
    N = params.N
    x = params.x
    # apply boundary conditions
    du[1] = 0.0
    du[end] = 0.0
    # evolution in the bulk/
    for i in 2:(N - 1)
        left_face = x[i] - 0.5h
        right_face = x[i] + 0.5h

        # Rj = reconstruct_func(x[i], U[i], U[i - 1], U[i + 1], params)
        # Rl = reconstruct_func(x[i], U[i - 1], U[i - 2], U[i], params)
        # Rr = reconstruct_func(x[i], U[i + 1], U[i], U[i + 2], params)
        Rj = reconstruct_func(x[i], U[i], U[i - 1], U[i + 1], params)
        Rl = reconstruct_func(x[i], U[i - 1], U[i - 1], U[i], params)
        Rr = reconstruct_func(x[i], U[i + 1], U[i], U[i + 1], params)

        ulr = Rj(left_face, U[i])
        ull = Rl(left_face, U[i - 1])

        urr = Rr(right_face, U[i + 1])
        url = Rj(right_face, U[i])

        Fl = numerical_flux(ull, ulr, equation)
        Fr = numerical_flux(url, urr, equation)

        # Fl = numerical_flux(ull, ulr, equation)
        # Fr = numerical_flux(url, urr, equation)
        du[i] = -(Fr - Fl) / h
    end

    return nothing
end

end
