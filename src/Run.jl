module Run

using ..Equations
using ..ODE
using ..InitialData
using ..Reconstructions
using ..NumericalFluxes
using ..Limiters

using GridFunctions
using RecursiveArrayTools
using OrdinaryDiffEq
using DataInterpolations
using DiffEqCallbacks

abstract type SimulationParameters{T} end

struct SimulationParametersPotential{T} <: SimulationParameters{T}
    L::Int64
    N::Int64
    dx::T
    x::Vector{T}
    q::T
    m::T
    x0::T
    v0::T
    sf::Bool
    λ::T
    cfl::T
end

struct SimulationParametersForced <: SimulationParameters{Any}
    L::Any
    N::Any
    dx::Any
    x::Any
    q::Any
    x0::Any
    A::Any
    ω::Any
    vmax::Any
    direction::Any
    Nosc::Any
    cfl::Any
end

struct SimulationParametersInteracting{T} <: SimulationParameters{T}
    L::Int64
    N::Int64
    dx::T
    x::Vector{T}
    q1::T
    m1::T
    x1::T
    v1::T
    q2::T
    m2::T
    x2::T
    v2::T
    cfl::T
end

struct Simulation{T}
    params::SimulationParameters{T}
    sol::SciMLBase.ODESolution
end

function forced_motion(Nosc, dx, cfl, vmax, q, sf=true, muscl=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    # vmax = 0.999
    # ω = 2π
    # T = 1.0
    # A = vmax / ω

    A = 1.0
    ω = round(vmax / A; digits=9)
    T = (2π / ω)

    L = ceil(4Nosc * T / 3)
    @show ω
    @show A
    @show L
    @show T
    N = 2Int64(L)
    N = Int64(ceil(N / dx))
    @show N
    grid = UniformStaggeredGrid1D(Float64[-L, L], N)
    x = coords(grid) # cell centers
    @show length(x)
    Π = zeros(Float64, length(grid))
    Ψ = zeros(Float64, length(grid))

    statevector = hcat(Π, Ψ)

    dx = spacing(grid)

    dt = round(cfl * dx / c; digits=8)
    @show dx, dt

    numerical_flux = NumericalFluxes.flux_rusanov
    slope_limiter = Limiters.minmod
    flux_limiter = Limiters.minmod

    # q = +1
    x10 = 0.0
    direction = +1

    tf = Nosc * T
    tspan = (0.0, tf)
    t = 0.0:dt:(tf + dt)
    @show length(t)
    sim = SimulationParametersForced(L, N, dx, x,
                                     q, x10, A,
                                     ω, vmax, direction, Nosc,
                                     cfl)

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L, vmax=vmax, A=A,
              dt=dt, x1=x10, ω=ω, q=q)
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
    res = Simulation(sim, sol)
    return res
end

function particle_given_potential(tf, dt)
    tspan = (0.0, tf)
    t = 0.0:dt:tf

    n = 50
    c = 1
    A = 1 / (2π)
    λ = 10
    Lwave = λ * n

    # Ψ(tt, xx) = InitialData.dxSineWave(tt, xx, n, c, A, Lwave)
    # Π(tt, xx) = InitialData.dtSineWave(tt, xx, n, c, A, Lwave)
    Ψ(tt, xx) = 10.0 #InitialData.dxSineWave(tt, xx, n, c, A, Lwave)
    Π(tt, xx) = 0.0 #InitialData.dtSineWave(tt, xx, n, c, A, Lwave)
    q1 = -1.0
    m10 = 1.0
    x10 = 0.0
    v10 = 0.0

    statevector = [m10, x10, v10] #, m20, x20, v20]
    params = (q1=q1, pifield=Π, psifield=Ψ)
    alg = RK4()
    ode = ODEProblem{true}(ODE.particle_rhs_functions!, statevector,
                           tspan, params; saveat=t,
                           #callback=cb, save_everystep=false
                           )
    sol = solve(ode, alg; adaptive=false, dt=dt)
    return sol
end

function particle_given_interpolation(L, dx, tf, dt)
    tspan = (0.0, tf)
    t = 0.0:dt:tf

    x = (-L):dx:L

    n = 50
    c = 1
    A = 1 / (2π)
    λ = 10
    Lwave = λ * n

    Ψ(tt, xx) = InitialData.dxSineWave(tt, xx, n, c, A, Lwave)
    Π(tt, xx) = InitialData.dtSineWave(tt, xx, n, c, A, Lwave)

    q1 = -1.0
    m10 = 1.0
    x10 = 0.0
    v10 = 0.0
    interpolation_method = QuadraticInterpolation
    statevector = [m10, x10, v10] #, m20, x20, v20]
    params = (q1=q1, pifield=Π.(0.0, x), psifield=Ψ.(0.0, x),
              interpolation_method=interpolation_method, x=x)
    alg = SSPRK54()
    ode = ODEProblem{true}(ODE.particle_rhs_interpolation!, statevector,
                           tspan, params; saveat=t)
    sol = solve(ode, alg; adaptive=false, dt=dt)
    return sol
end

function particle_in_potential(L, N, tf, cfl, q, sf=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    grid = UniformStaggeredGrid1D(Float64[-L, L], N)
    x = GridFunctions.Grids.coords(grid) # cell centers
    tspan = (0.0, tf)
    dx = spacing(grid)

    dt = cfl * dx / c
    @show dx, dt
    t = 0.0:dt:tf

    # Π = zeros(Float64, length(grid))
    # Ψ = zeros(Float64, length(grid))

    #n = 5
    A = 1 / (2π)
    λ = 200.0 # L / n
    @show 2A * π / λ
    @show λ
    Π = zeros(Float64, length(grid))

    Ψ = InitialData.dxSineWave.(0.0, x, A, λ, c)
    # Ψ = InitialData.dxHarmonic(0.0, x, 0.005)
    # Π = InitialData.dtSineWave.(0.0, x, A, λ, c)

    field_statevector = hcat(Π, Ψ)

    ## Particle 1
    # q = 0.05
    m10 = 1.0
    x10 = 0.0
    v10 = 0.0

    println("Particle:")
    println("q = $q")
    println("m = $m10")
    println("x = $x10")
    println("v = $v10")

    particle_statevector = [m10, x10, v10]#, m20, x20, v20]
    statevector = ArrayPartition(field_statevector,
                                 particle_statevector)

    numerical_flux = NumericalFluxes.flux_rusanov
    slope_limiter = Limiters.minmod
    flux_limiter = Limiters.minmod
    interpolation_method = QuadraticInterpolation # LinearInterpolation # QuadraticInterpolation

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L,
              dt=dt, x1=x10,
              q=q, sf=sf,
              pifield=Π, psifield=Ψ,
              interpolation_method=interpolation_method)
    sim = SimulationParametersPotential(L, N, dx, x, q,
                                        m10, x10, v10,
                                        sf, λ, cfl)
    alg = SSPRK54()
    ode = ODEProblem{true}(ODE.coupled_rhs!, statevector, tspan,
                           params; saveat=t)

    sol = solve(ode, alg; dt=dt)#, callback=cb)
    res = Simulation(sim, sol)
    return res
end

function interacting_particles(L, N, tf, cfl, v0, q0, same_q=true, sf=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    grid = UniformStaggeredGrid1D(Float64[-L, L], N)
    x = GridFunctions.Grids.coords(grid) # cell centers
    tspan = (0.0, tf)
    dx = spacing(grid)

    dt = cfl * dx / c
    @show dx, dt
    t = 0.0:dt:tf

    Π = zeros(Float64, length(grid))
    Ψ = zeros(Float64, length(grid))

    field_statevector = hcat(Π, Ψ)

    ## Particle 1
    q1 = Float64(q0)
    m10 = 1.0
    x10 = -L / 4
    v10 = v0
    ## Particle 2
    if same_q
        q2 = Float64(q0)
    else
        q2 = -Float64(q0)
    end
    m20 = 1.0
    x20 = L / 4
    v20 = -v0

    println("Particle 1:         Particle 2:")
    println("q = $q1             q = $q2")
    println("m = $m10            m = $m20")
    println("x = $x10            x = $x20")
    println("v = $v10            v = $v20")

    particle_statevector = [m10, x10, v10, m20, x20, v20]
    statevector = ArrayPartition(field_statevector,
                                 particle_statevector)

    numerical_flux = NumericalFluxes.flux_rusanov
    slope_limiter = Limiters.minmod
    flux_limiter = Limiters.minmod
    interpolation_method = ConstantInterpolation # LinearInterpolation # QuadraticInterpolation

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L,
              dt=dt, x1=x10, x2=x20,
              q1=q1, q2=q2, sf=sf,
              interpolation_method=interpolation_method, acc=[1.0])

    sim = SimulationParametersInteracting(L, N, dx, x,
                                          q1, m10, x10, v10,
                                          q2, m20, x20, v20, cfl)
    alg = SSPRK54()
    ode = ODEProblem{true}(ODE.interacting_coupled_rhs!, statevector, tspan,
                           params; saveat=t,
                           #callback=callbacks
                           )

    sol = solve(ode, alg; dt=dt)#, callback=cb)
    res = Simulation(sim, sol)
    return res
end

function coupled_system_fractional(L, N, tf, cfl, sf=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    grid = UniformStaggeredGrid1D(Float64[-L, L], N)
    x = GridFunctions.Grids.coords(grid) # cell centers
    tspan = (0.0, tf)
    dx = spacing(grid)

    dt = cfl * dx / c
    @show dx, dt
    t = 0.0:dt:tf

    Π = zeros(Float64, length(grid))
    Ψ = zeros(Float64, length(grid))

    field_statevector = hcat(Π, Ψ)

    ## Particle 1
    q1 = 0.04
    m10 = 1.0
    x10 = -L / 4
    v10 = 0.5
    ## Particle 2
    q2 = 0.04
    m20 = 1.0
    x20 = L / 4
    v20 = -0.5

    println("Particle 1:         Particle 2:")
    println("q = $q1             q = $q2")
    println("m = $m10            m = $m20")
    println("x = $x10            x = $x20")
    println("v = $v10            v = $v20")

    particle_statevector = [m10, x10, v10, m20, x20, v20]
    statevector = ArrayPartition(field_statevector,
                                 particle_statevector)

    numerical_flux = NumericalFluxes.flux_godunov
    slope_limiter = Limiters.upwind
    flux_limiter = Limiters.minmod
    interpolation_method = QuadraticInterpolation # LinearInterpolation # QuadraticInterpolation

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L,
              dt=dt, x1=x10, x2=x20,
              q1=q1, q2=q2, sf=sf,
              interpolation_method=interpolation_method)
    sim = SimulationParameters2(L, N, dx, x,
                                q1, m10, x10, v10,
                                q2, m20, x20, v20, cfl)
    # alg = SSPRK54()
    # alg = ImplicitMidpoint()

    alg = KenCarp4()
    alg = IMEXEuler()
    ode = SplitODEProblem{true}(ODE.particle_rhs!, ODE.field_rhs!, statevector, tspan,
                                params; saveat=t)

    # return statevector, params
    sol = solve(ode, alg; adaptive=false, dt=dt)
    # PotentialEnergy = q1 * f(0.0, x10)
    # KineticEnergyNR = 0.5 * m10 * v10^2
    # KineticEnergyRel = (sqrt(1 / (1 - v10^2)) - 1) * m10 # (γ-1)m
    # Energy = KineticEnergyRel + PotentialEnergy
    # @show Energy, PotentialEnergy, KineticEnergyNR, KineticEnergyRel
    # # @assert Energy < 0
    res = Simulation(sim, sol)
    return res
end

end # end of module
