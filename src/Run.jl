module Run

using ..Equations
using ..ODE
using ..InitialData
using ..Reconstructions
using ..NumericalFluxes
using ..Limiters

using GridFunctions
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqLowOrderRK

function run(L, N, tf, cfl, muscl=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    grid = UniformStaggeredGrid1D(Float64[0, L], N)
    x = coords(grid) # cell centers

    Π = zeros(Float64, length(grid))
    Ψ = zeros(Float64, length(grid))

    # Ψ[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= -1.0 # -1 for right moving, +1 for left moving
    # Π[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= 1.0

    #  I get a reflection at the boundary with these?!
    # Ψ = InitialData.dx_triangular_pulse.(x, 0.0, L, 1.0, 20, c)
    # Π = InitialData.dt_triangular_pulse.(x, 0.0, L, 1.0, 20, c)

    σ = 30.0
    A = 1.0
    n = 2
    Ψ .= InitialData.dxSineWave.(0.0, x, n, c, A, L)
    Π .= InitialData.dtSineWave.(0.0, x, n, c, A, L)

    # Ψ .= InitialData.dxGaussian1D.(0.0, x, A, σ)
    # Π .= InitialData.dtGaussian1D.(0.0, x, A, σ)
    statevector = hcat(Π, Ψ)

    tspan = (0.0, tf)
    dx = spacing(grid)

    dt = cfl * dx / c
    @show dx, dt
    t = 0.0:dt:tf
    numerical_flux = NumericalFluxes.flux_godunov
    slope_limiter = Limiters.MC
    flux_limiter = Limiters.superbee

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl,
              dt=dt, x1=L / 2, x2=-L / 2,
              direction1=+1, direction2=-1,
              q1=1, q2=-1)
    if muscl
        # alg = SSPRK54()
        alg = SSPRK22()
        ode = ODEProblem{true}(ODE.muscl!, statevector, tspan, params; saveat=t)
    else
        alg = Euler()
        ode = ODEProblem{true}(ODE.LxW!, statevector, tspan, params)
    end

    # return statevector, params
    sol = solve(ode, alg; adaptive=false, dt=dt)
    return GridFunctions.coords(grid), sol
end

end # end of module
