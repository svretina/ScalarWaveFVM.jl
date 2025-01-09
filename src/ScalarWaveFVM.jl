module ScalarWaveFVM

include("Equations.jl")
include("Reconstructions.jl")
include("NumericalFluxes.jl")

using ..Reconstructions
using ..NumericalFluxes

using OrdinaryDiffEqSSPRK
using StaticArrays

function run(L, N, tf, cfl)
    xs = LinRange(-L, L, N + 1)

    Π = zeros(Float64, N + 1)
    Ψ = zeros(Float64, N + 1)
    Ψ[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= 1.0

    tspan = (0.0, tf)
    @show tf
    dx = (xs[2] - xs[1])
    dt = cfl * dx
    alg = SSPRK54()

    equation = LinearScalarWavequation1D(1.0)
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

    # evolution in the bulk/
    half_h = h * 0.5
    for i in 2:(N - 1)
        left_face = x[i] - half_h
        right_face = x[i] + half_h

        Ri = reconstruct_func(x[i], U[i], U[i], U[i], params)
        Rl = reconstruct_func(x[i], U[i], U[i], U[i], params)
        Rr = reconstruct_func(x[i], U[i], U[i], U[i], params)

        ulr = Ri(left_face, U[i])
        ull = Rl(left_face, U[i - 1])

        urr = Rr(right_face, U[i + 1])
        url = Ri(right_face, U[i])

        Fl = numerical_flux(ull, ulr, equation, +1)
        Fr = numerical_flux(url, urr, equation, -1)

        # Fl = numerical_flux(ull, ulr, equation)
        # Fr = numerical_flux(url, urr, equation)
        du[i] = -(Fr - Fl) / h
    end
    return nothing
end

end
