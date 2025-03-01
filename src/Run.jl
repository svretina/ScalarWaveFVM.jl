module Run

using ..Equations
using ..ODE
using ..InitialData
using ..Reconstructions
using ..NumericalFluxes
using ..Limiters
# using ..FractionalStepMethods

using GridFunctions
using RecursiveArrayTools
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqLowOrderRK
using DataInterpolations
using DiffEqCallbacks

function save()
    # SavingCallback
end

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
    q2 = -1
    direction1 = +1
    direction2 = +1

    # check for signs
    # either qs have same sign and directions not
    # or qs have different sign and directions have same sign
    cond1 = (q1 * q2 > 0 && direction1 * direction2 < 0)
    cond2 = (q1 * q2 < 0 && direction1 * direction2 > 0)

    @assert cond1 || cond2

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L, f₀=2, vmax=0.9, A=1000 / 6,
              dt=dt, x1=500, x2=-500,
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

function particle_given_potential(tf, dt)
    tspan = (0.0, tf)
    t = 0.0:dt:tf

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

    statevector = [m10, x10, v10] #, m20, x20, v20]
    params = (q1=q1, pifield=Π, psifield=Ψ)
    alg = SSPRK54()
    ode = ODEProblem{true}(ODE.particle_rhs_functions!, statevector,
                           tspan, params; saveat=t)
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

function coupled_system(L, N, tf, cfl, sf=true)
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
    # Ψ = zeros(Float64, length(grid))

    n = 5
    A = 1 / (2π)
    λ = 80 # L / n
    @show 2A * π / λ
    @show λ
    Ψ = InitialData.dxSineWave.(0.0, x, A, λ, c)
    # Π = InitialData.dtSineWave.(0.0, x, 5, c, 1.0, 50)

    field_statevector = hcat(Π, Ψ)

    ## Particle 1
    q1 = 0.08
    # @assert abs(q1) > dt
    m10 = 1.0
    x10 = 0.0
    v10 = 0.0

    println("Particle:")
    println("q = $q1")
    println("m = $m10")
    println("x = $x10")
    println("v = $v10")

    particle_statevector = [m10, x10, v10]#, m20, x20, v20]
    statevector = ArrayPartition(field_statevector,
                                 particle_statevector)

    numerical_flux = NumericalFluxes.flux_godunov
    slope_limiter = Limiters.upwind
    flux_limiter = Limiters.minmod
    interpolation_method = QuadraticInterpolation # LinearInterpolation # QuadraticInterpolation

    function condition(u, t, integrator)
        # Check if either particle is outside the domain
        particle1_out = u.x[2][2] >= L || u.x[2][2] < -L
        # particle2_out = u.x[2][5] >= L || u.x[2][5] < -L
        return particle1_out #|| particle2_out
    end

    function affect!(integrator)
        # Check and adjust first particle
        if integrator.u.x[2][2] >= L
            integrator.u.x[2][2] -= 2L
        elseif integrator.u.x[2][2] < -L
            integrator.u.x[2][2] += 2L
        end

        # Check and adjust second particle
        if integrator.u.x[2][5] >= L
            integrator.u.x[2][5] -= 2L
        elseif integrator.u.x[2][5] < -L
            integrator.u.x[2][5] += 2L
        end
    end

    # Define and apply the callback
    cb = DiscreteCallback(condition, affect!)

    params = (equation=equation, numerical_flux=numerical_flux,
              slope_limiter=slope_limiter, flux_limiter=flux_limiter,
              h=dx, N=N, x=x, cfl=cfl, L=L,
              dt=dt, x1=x10,
              q1=q1, sf=sf,
              pifield=Π, psifield=Ψ,
              interpolation_method=interpolation_method)

    alg = SSPRK54()
    ode = ODEProblem{true}(ODE.coupled_rhs!, statevector, tspan,
                           params; saveat=t)

    sol = solve(ode, alg; adaptive=false, dt=dt)#, callback=cb)
    return GridFunctions.coords(grid), sol
end

# function coupled_system_fractional(L, N, tf, cfl)
#     c = 1.0
#     equation = LinearScalarWaveEquation1D(c)

#     grid = UniformStaggeredGrid1D(Float64[-L, L], N)
#     x = coords(grid) # cell centers

#     Π = zeros(Float64, length(grid))
#     Ψ = zeros(Float64, length(grid))

#     field_statevector = hcat(Π, Ψ)
#     q1 = -1.0
#     q2 = -1.0
#     m10 = 1.0
#     m20 = 1.0
#     x10 = -L / 2
#     x20 = L / 2
#     v10 = 0.1
#     v20 = -0.1
#     cond1 = (q1 * q2 > 0 && sign(v10) * sign(v20) < 0)
#     cond2 = (q1 * q2 < 0 && sign(v10) * sign(v20) > 0)
#     @assert cond1 || cond2
#     # check for signs
#     # either qs have same sign and directions not
#     # or qs have different sign and directions have same sign

#     particle_statevector = [m10, x10, v10, m20, x20, v20]
#     statevector = ArrayPartition(field_statevector, particle_statevector)

#     tspan = (0.0, tf)
#     dx = spacing(grid)

#     dt = cfl * dx / c
#     @show dx, dt
#     t = 0.0:dt:tf
#     numerical_flux = NumericalFluxes.flux_godunov
#     slope_limiter = Limiters.minmod
#     flux_limiter = Limiters.minmod

#     params = (equation=equation, numerical_flux=numerical_flux,
#               slope_limiter=slope_limiter, flux_limiter=flux_limiter,
#               h=dx, N=N, x=x, cfl=cfl, L=L,
#               dt=dt, x1=L / 2, x2=-L / 2,
#               direction1=direction1, direction2=direction2,
#               q1=q1, q2=q2, pifield=Π, psifield=Ψ)

#     alg1 = SSPRK54()
#     alg2 = Euler()
#     alg = FractionalStepMethods.MyCustonSolver(alg1, alg2; splitting=Strang)

#     ## TODO
#     ## implement source! function in ODE.jl
#     ## take care of dependencies!

#     ode = SplitODEProblem{true}(ODE.muscl!, ODE.source!, statevector, tspan,
#                                 params; saveat=t)

#     # return statevector, params
#     sol = solve(ode, alg; adaptive=false, dt=dt)
#     # PotentialEnergy = q1 * f(0.0, x10)
#     # KineticEnergyNR = 0.5 * m10 * v10^2
#     # KineticEnergyRel = (sqrt(1 / (1 - v10^2)) - 1) * m10 # (γ-1)m
#     # Energy = KineticEnergyRel + PotentialEnergy
#     # @show Energy, PotentialEnergy, KineticEnergyNR, KineticEnergyRel
#     # # @assert Energy < 0
#     return GridFunctions.coords(grid), sol
# end

end # end of module
