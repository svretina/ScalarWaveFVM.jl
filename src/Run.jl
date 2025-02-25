module Run

using ..Equations
using ..ODE
using ..InitialData
using ..Reconstructions
using ..NumericalFluxes
using ..Limiters

using GridFunctions
using RecursiveArrayTools
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqLowOrderRK

function run(L, N, tf, cfl, sf=false, muscl=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    grid = UniformStaggeredGrid1D(Float64[-L, L], N)
    x = coords(grid) # cell centers

    Π = zeros(Float64, length(grid))
    Ψ = zeros(Float64, length(grid))

    # Ψ[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= -1.0 # -1 for right moving, +1 for left moving
    # Π[(div(N, 4) + 1):(3 * div(N, 4) + 1)] .= 1.0

    #  I get a reflection at the boundary with these?!
    # Ψ = InitialData.dx_triangular_pulse.(x, 0.0, L, 1.0, 20, c)
    # Π = InitialData.dt_triangular_pulse.(x, 0.0, L, 1.0, 20, c)

    # σ = 30.0
    # A = 1.0
    # n = 2
    # Ψ .= InitialData.dxSineWave.(0.0, x, n, c, A, L)
    # Π .= InitialData.dtSineWave.(0.0, x, n, c, A, L)

    # Ψ .= InitialData.dxGaussian1D.(0.0, x, A, σ)
    # Π .= InitialData.dtGaussian1D.(0.0, x, A, σ)
    statevector = hcat(Π, Ψ)

    tspan = (0.0, tf)
    dx = spacing(grid)

    dt = cfl * dx / c
    @show dx, dt
    t = 0.0:dt:tf
    numerical_flux = NumericalFluxes.flux_godunov
    slope_limiter = Limiters.minmod
    flux_limiter = Limiters.minmod

    q1 = +1
    q2 = +1                        
    direction1 = +1
    direction2 = -1

    # check for signs
    # either qs have same sign and directions not
    # or qs have different sign and directions have same sign
    cond1 = (q1 * q2 > 0 && direction1 * direction2 < 0)
    cond2 = (q1 * q2 < 0 && direction1 * direction2 > 0)

    @assert cond1 || cond2

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L,
              dt=dt, x1=L / 2, x2=-L / 2,
              direction1=direction1, direction2=direction2,
              q1=q1, q2=q2)
    if muscl
        alg = SSPRK54()
        # alg = SSPRK22()
        if sf
            ode = ODEProblem{true}(ODE.field_rhs_forced!,
                                   statevector, tspan, params;
                                   saveat=t)
        else
            ode = ODEProblem{true}(ODE.muscl!, statevector, tspan, params; saveat=t)
        end
    else
        alg = Euler()
        ode = ODEProblem{true}(ODE.LxW!, statevector, tspan, params)
    end

    # return statevector, params
    sol = solve(ode, alg; adaptive=false, dt=dt)
    return GridFunctions.coords(grid), sol
end

function particle_given_potential(L, tf, dt)
    tspan = (0.0, tf)
    t = 0.0:dt:tf
    h = 1.0
    x1 = L / 2 + h / 2
    x2 = -L / 2 + h / 2
    x = (-L):h:L
    f(t, x) = 1 / (1 + (x - x1)^2) + 1 / (1 + (x - x2)^2)
    ft(t, x) = zero(x)
    function fx(t, x)
        tmp = (2(x - x1) / (one(x) + 0.5(x - x1)^2)^2) +
              (2(x - x2) / (one(x) + 0.5(x - x2)^2)^2)
        return tmp
    end

    Π = ft.(0.0, x)
    Ψ = fx.(0.0, x)
    q1 = -1.0
    q2 = -1.0

    m10 = 1.0
    m20 = 1.0
    x10 = x1
    x20 = x2
    v10 = 0.1
    v20 = 0.1

    # PotentialEnergy = q1 * f(0.0, x10)
    # KineticEnergyNR = 0.5 * m10 * v10^2
    # KineticEnergyRel = (sqrt(1 / (1 - v10^2)) - 1) * m10 # (γ-1)m
    # Energy = KineticEnergyRel + PotentialEnergy
    # @show Energy, PotentialEnergy, KineticEnergyNR, KineticEnergyRel
    # # @assert Energy < 0

    statevector = [m10, x10, v10, m20, x20, v20]
    params = (q1=q1, q2=q2, pifield=Π, psifield=Ψ, x=x)
    alg = SSPRK54()
    ode = ODEProblem{true}(ODE.particle_motion!, statevector,
                           tspan, params; saveat=t)
    sol = solve(ode, alg; adaptive=false, dt=dt)
    return sol
end

function coupled_system(L, N, tf, cfl)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    grid = UniformStaggeredGrid1D(Float64[-L, L], N)
    x = coords(grid) # cell centers

    Π = zeros(Float64, length(grid))
    Ψ = zeros(Float64, length(grid))

    field_statevector = hcat(Π, Ψ)
    q1 = -1.0
    q2 = -1.0
    m10 = 1.0
    m20 = 1.0
    x10 = -L / 2
    x20 = L / 2
    v10 = 0.1
    v20 = -0.1
    cond1 = (q1 * q2 > 0 && sign(v10) * sign(v20) < 0)
    cond2 = (q1 * q2 < 0 && sign(v10) * sign(v20) > 0)
    @assert cond1 || cond2
    # check for signs
    # either qs have same sign and directions not
    # or qs have different sign and directions have same sign
                                   
    particle_statevector = [m10, x10, v10, m20, x20, v20]
    statevector = ArrayPartition(field_statevector, particle_statevector)

    tspan = (0.0, tf)
    dx = spacing(grid)

    dt = cfl * dx / c
    @show dx, dt
    t = 0.0:dt:tf
    numerical_flux = NumericalFluxes.flux_godunov
    slope_limiter = Limiters.minmod
    flux_limiter = Limiters.minmod

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L,
              dt=dt, x1=L / 2, x2=-L / 2,
              direction1=direction1, direction2=direction2,
              q1=q1, q2=q2, pifield=Π, psifield=Ψ)

    alg = SSPRK54()
    ode = ODEProblem{true}(ODE.coupled_rhs!, statevector, tspan,
                           params; saveat=t)

    # return statevector, params
    sol = solve(ode, alg; adaptive=false, dt=dt)
    # PotentialEnergy = q1 * f(0.0, x10)
    # KineticEnergyNR = 0.5 * m10 * v10^2
    # KineticEnergyRel = (sqrt(1 / (1 - v10^2)) - 1) * m10 # (γ-1)m
    # Energy = KineticEnergyRel + PotentialEnergy
    # @show Energy, PotentialEnergy, KineticEnergyNR, KineticEnergyRel
    # # @assert Energy < 0
    return GridFunctions.coords(grid), sol
end

end # end of module
