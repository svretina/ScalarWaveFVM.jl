module Run

using RecursiveArrayTools
using ..Equations
using ..ODE
using ..InitialData
using ..Reconstructions
using ..NumericalFluxes
using ..Limiters
using ..ScalarField

using GridFunctions
using RecursiveArrayTools
using OrdinaryDiffEq
using DataInterpolations
# using DiffEqCallbacks

abstract type SimulationParameters{T} end

struct SimulationParameters1D{T} <: SimulationParameters{T}
    L::Int64
    N::Int64
    dx::T
    dt::T
    x::Vector{T}
    q::T
    m::T
    z0::T
    v0::T
    vmax::T
    A::T
    ω::T
    k::T
    direction::Int
    Nosc::T
    sf::Bool
    cfl::T
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

function print_simulation_info(title; params...)
    println("\n" * "="^50)
    println("Running Simulation: $title")
    println("-"^50)
    for (key, value) in params
        println(rpad(string(key), 20) * ": " * string(value))
    end
    println("="^50 * "\n")
end

function print_simulation_info(sim::SimulationParameters1D; title="")
    println("\n" * "="^50)
    if isempty(title)
        if sim.k == 0.0
            println("Running Simulation: Forced Motion")
            SimType = "Forced Motion"
        else
            println("Running Simulation: Particle in Potential")
            SimType = "Harmonic Potential"
        end
    else
        println("Running Simulation: $title")
        SimType = title
    end
    println("-"^50)

    T = 2π / sim.ω
    tf = sim.Nosc * T

    println(rpad("SimType", 20) * ": " * SimType)
    println(rpad("L", 20) * ": " * string(sim.L))
    println(rpad("N", 20) * ": " * string(sim.N))
    println(rpad("dx", 20) * ": " * string(sim.dx))
    println(rpad("dt", 20) * ": " * string(sim.dt))
    println(rpad("cfl", 20) * ": " * string(sim.cfl))
    println(rpad("Nosc", 20) * ": " * string(sim.Nosc))
    println(rpad("tf", 20) * ": " * string(tf))
    println(rpad("T", 20) * ": " * string(T))
    println(rpad("vmax", 20) * ": " * string(sim.vmax))
    println(rpad("q", 20) * ": " * string(sim.q))
    println(rpad("A", 20) * ": " * string(sim.A))
    println(rpad("omega", 20) * ": " * string(sim.ω))
    if sim.k != 0.0
        println(rpad("k", 20) * ": " * string(sim.k))
        println(rpad("m", 20) * ": " * string(sim.m))
    end

    println("="^50 * "\n")
end

function forced_motion(Nosc, dx, cfl, vmax, q, sf=true, muscl=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    A = 0.01
    ω = round(vmax / A; digits=9)
    T = 2π / ω

    max_retries = 5
    for attempt in 1:max_retries
        println("==========================================")
        println("Attempt $attempt / $max_retries (cfl = $cfl)")

        L = ceil(2Nosc * T)
        N = 2Int64(L)
        N = Int64(ceil(N / dx))

        grid = UniformStaggeredGrid1D(Float64[-L, L], N)
        x = coords(grid) # cell centers

        Φ = zeros(Float64, length(grid))
        Π = zeros(Float64, length(grid))
        Ψ = zeros(Float64, length(grid))

        statevector = hcat(Π, Ψ, Φ)
        dx = spacing(grid)

        dt = cfl * dx / c

        numerical_flux = NumericalFluxes.flux_godunov

        slope_limiter = Limiters.MC
        flux_limiter = Limiters.upwind

        z0 = 0.0
        direction = +1

        tf = Nosc * T
        L0 = ceil(2Nosc * T)

        tspan = (0.0, tf)
        t = 0.0:dt:(tf + dt)
        sim = SimulationParameters1D(Int64(L), Int64(N), float(dx), float(dt),
                                     convert(Vector{Float64}, x), float(q), 1.0, float(z0),
                                     float(vmax), float(vmax), float(A), float(ω),
                                     0.0, direction, float(Nosc), sf, float(cfl))

        params = (equation=equation, numerical_flux=numerical_flux,
                  slope_limiter=slope_limiter, flux_limiter=flux_limiter,
                  h=dx, N=N, x=x, cfl=cfl, L=L0, L0=L0, vmax=vmax, A=A,
                  dt=dt, z0=z0, ω=ω, q=q)

        print_simulation_info(sim; title="Forced Motion")
        if muscl
            alg = SSPRK54()
            if sf
                println("Self Force = True")
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
        if sol.retcode == :Success || sol.retcode == ReturnCode.Success
            res = Simulation(sim, sol)
            return res
        else
            println("Warning: Instability detected (retcode: $(sol.retcode)). Retrying")
        end
    end
    error("Simulation failed after $max_retries attempts")
end

function particle_in_potential(Nosc, dx, cfl, vmax, q, sf=true)
    c = 1.0
    equation = LinearScalarWaveEquation1D(c)

    max_retries = 5
    for attempt in 1:max_retries
        println("==========================================")
        println("Attempt $attempt / $max_retries (cfl = $cfl)")

        A = 0.01
        ω = round(vmax / A; digits=9)
        T = 2π / ω

        tf = Nosc * T
        L0 = 2tf
        L = Int(ceil(L0))
        tf = Nosc * T
        L0 = 2tf
        L = Int(ceil(L0))
        @assert L >= L0

        tspan = (0.0, tf)
        dt = cfl * dx / c
        t = 0.0:dt:tf
        tspan = (0.0, tf)
        dt = cfl * dx / c
        t = 0.0:dt:tf

        N = Int64(2L / dx)

        grid = UniformStaggeredGrid1D(Float64[-L, L], N)
        x = coords(grid) # cell centers

        @assert spacing(grid) == dx

        N = length(x)

        Φ = zeros(Float64, length(x))
        Π = zeros(Float64, length(x))
        Ψ = zeros(Float64, length(x))

        field_statevector = hcat(Π, Ψ, Φ)

        ## Particle 1
        m0 = 1.0
        z0 = 0.0
        v0 = vmax
        # z0 = 0.0
        # v0 = 0.6283185307179586

        k = (m0 / q) * ω^2

        V(x) = InitialData.Harmonic(x, k)
        dxV(x) = InitialData.dxHarmonic(x, k)

        particle_statevector = [m0, z0, v0]
        statevector = ArrayPartition(field_statevector, particle_statevector)

        numerical_flux = NumericalFluxes.flux_rusanov
        slope_limiter = Limiters.MC
        flux_limiter = Limiters.MC
        interpolation_method = AkimaInterpolation # QuadraticInterpolation # LinearInterpolation # QuadraticInterpolation

        params = (equation=equation, numerical_flux=numerical_flux,
                  slope_limiter=slope_limiter, flux_limiter=flux_limiter,
                  h=dx, N=N, x=x, dx=dx, cfl=cfl, L=L,
                  dt=dt, z0=z0,
                  q=q, sf=sf,
                  Potential=V,
                  dxPotential=dxV,
                  interpolation_method=interpolation_method,
                  m0=m0,
                  L0=L0)#, acc=Number[0.0])

        sim = SimulationParameters1D(Int64(L), Int64(N), float(dx), float(dt),
                                     convert(Vector{Float64}, x), float(q), float(m0),
                                     float(z0),
                                     float(v0), float(vmax), float(A), float(ω),
                                     float(k), 0, float(Nosc), sf, float(cfl))

        print_simulation_info(sim; title="Particle in Potential")

        alg = SSPRK54()

        ode = ODEProblem{true}(ODE.radiation_reaction_rhs!, statevector, tspan,
                               params; saveat=t)
        integrator = init(ode, alg; adaptive=false, dt=dt, p=params)
        for (uprev, tprev, u, t) in intervals(integrator)
            if iszero(t % 10)
                @show t
            end
        end

        if integrator.sol.retcode == :Success ||
           integrator.sol.retcode == ReturnCode.Success
            # sim defined earlier

            res = Simulation(sim, integrator.sol)
            return res
        else
            println("Warning: Instability detected (retcode: $(integrator.sol.retcode)). Retrying...")
        end
    end
    error("Simulation failed after $max_retries attempts")
end

# function particle_given_potential(tf, dt)
#     tspan = (0.0, tf)
#     t = 0.0:dt:tf

#     n = 50
#     c = 1
#     A = 1 / (2π)
#     λ = 10
#     Lwave = λ * n

#     # Ψ(tt, xx) = InitialData.dxSineWave(tt, xx, n, c, A, Lwave)
#     # Π(tt, xx) = InitialData.dtSineWave(tt, xx, n, c, A, Lwave)
#     a = 0.25
#     Ψ(tt, xx) = InitialData.dxHarmonic(tt, xx, a)
#     Π(tt, xx) = InitialData.dtHarmonic(tt, xx, a)
#     q1 = -1.0
#     m10 = 1.0
#     x10 = 0.0
#     v10 = 0.0

#     statevector = [m10, x10, v10] #, m20, x20, v20]
#     params = (q1=q1, pifield=Π, psifield=Ψ)
#     alg = RK4()
#     ode = ODEProblem{true}(ODE.particle_rhs_functions!, statevector,
#                            tspan, params; saveat=t,
#                            #callback=cb, save_everystep=false
#                            )
#     sol = solve(ode, alg; adaptive=false, dt=dt)
#     return sol
# end

# function particle_in_fixed_potential_interpolation(L, dx, dt, Nosc, k)
#     # n = 50
#     # c = 1
#     # A = 1 / (2π)
#     # λ = 10
#     # Lwave = λ * n

#     #     k = 1.0
#     T = 2π / sqrt(k)
#     tf = Nosc * T
#     tspan = (0.0, tf)
#     t = 0.0:dt:tf

#     x = (-L):dx:L

#     Ψ(tt, xx) = InitialData.dxHarmonic(tt, xx, k)
#     Π(tt, xx) = InitialData.dtHarmonic(tt, xx, k)

#     q1 = 1.0
#     m10 = 1.0
#     x10 = 0.1
#     v10 = 0.0
#     interpolation_method = QuadraticInterpolation
#     statevector = [m10, x10, v10] #, m20, x20, v20]
#     params = (q1=q1, pifield=Π, psifield=Ψ,
#               interpolation_method=interpolation_method, x=x)

#     alg = SSPRK54() #
#     ode = ODEProblem{true}(ODE.particle_rhs_interpolation!, statevector,
#                            tspan, params; saveat=t)
#     sol = solve(ode, alg; adaptive=false, dt=dt)
#     sim = SimulationParameters1D(Int64(L), Int64(length(x)), float(dx), float(dt), 
#                                  convert(Vector{Float64}, collect(x)), float(q1), float(m10), float(x10), 
#                                  float(v10), float(v10), 0.0, 0.0, 
#                                  float(k), 0, float(Nosc), true, 0.0)
#     res = Simulation(sim, sol)
#     return res
# end

# function particle_in_potential_split(L, dx, dt, Nosc, k, sf=true)
#     c = 1.0
#     equation = LinearScalarWaveEquation1D(c)

#     T = 2π / sqrt(abs(k))
#     tf = Nosc * T
#     L0 = 2tf
#     @show tf
#     @assert L > L0
#     @show L
#     tspan = (0.0, tf)
#     t = 0.0:dt:tf
#     N = Int64(2L / dx)
#     grid = UniformStaggeredGrid1D(Float64[-L, L], N)
#     x = coords(grid) # cell centers
#     @assert spacing(grid) == dx
#     N = length(x)
#     cfl = dt / dx
#     @show cfl

#     Φ = zeros(Float64, length(x))
#     Π = zeros(Float64, length(x))
#     Ψ = zeros(Float64, length(x))

#     V(x) = InitialData.Harmonic(x, k)
#     dxV(x) = InitialData.dxHarmonic(x, k)
#     field_statevector = hcat(Π, Ψ, Φ)

#     ## Particle 1
#     q = 0.1
#     m10 = 1.0
#     x10 = 1.0
#     v10 = 0.0

#     println("Particle:")
#     println("q = $q")
#     println("m = $m10")
#     println("x = $x10")
#     println("v = $v10")

#     particle_statevector = [m10, x10, v10] #, m20, x20, v20]
#     statevector = ArrayPartition(field_statevector, particle_statevector)

#     numerical_flux = NumericalFluxes.flux_rusanov
#     slope_limiter = Limiters.minmod
#     flux_limiter = Limiters.minmod
#     interpolation_method = AkimaInterpolation # QuadraticInterpolation # LinearInterpolation # QuadraticInterpolation

#     params = (equation=equation, numerical_flux=numerical_flux,
#               slope_limiter=slope_limiter, flux_limiter=flux_limiter,
#               h=dx, N=N, x=x, dx=dx, cfl=cfl, L=L,
#               dt=dt, x1=x10,
#               q=q, sf=sf,
#               Potential=V,
#               dxPotential=dxV,
#               interpolation_method=interpolation_method, L0=L0, acc=Number[0.0])

#     # alg = ImplicitMidpoint()
#     # alg = ImplicitEuler()
#     # alg = RK4() ## loses energy
#     # alg = Tsit5() ## fails
#     # alg = TsitPap8() ## ok
#     # alg = Euler()
#     # alg = Vern7() ## ok
#     # alg = SSPRK54()
#     # alg = AutoTsit5(Rosenbrock23())
#     # alg = SSPSDIRK2()
#     # alg = Trapezoid()

#     field_alg = Euler()
#     particle_alg = ImplicitMidpoint()
#     # ODE.particle_with_potential_rhs!
#     particle_ode = ODEProblem{true}(ODE.particle_with_potential_rhs!,
#                                     statevector, tspan,
#                                     params; saveat=t)

#     field_ode = ODEProblem{true}(ODE.field_rhs!,
#                                  statevector, tspan,
#                                  params; saveat=t)

#     particle_integrator = init(particle_ode, particle_alg; adaptive=false, dt=dt, p=params)
#     field_integrator = init(field_ode, field_alg; adaptive=false, dt=dt, p=params)
#     for (i, t) in enumerate(t)
#         @show t
#         set_u!(particle_integrator, field_integrator.u)
#         step!(particle_integrator, dt, false)

#         set_u!(field_integrator, particle_integrator.u)
#         step!(field_integrator, dt, false)
#     end
#     # for (uprev, tprev, u, t) in intervals(integrator)
#     #     # if t > 0
#     #     #     dt_new = (atanh(1.0) - abs(uprev.x[2][3])) / abs(integrator.p.acc[1])
#     #     # end
#     #     @show t
#     #     a = integrator.p.acc[1]
#     #     if abs(a) > 0.1
#     #         @show t, a
#     #     end
#     #     # set_proposed_dt!(integrator, dt_new)
#     # end
#     # sol = solve(ode, alg; dt=dt, adaptive=false)
#     sim = SimulationParameters1D(Int64(L), Int64(length(x)), float(dx), float(dt), 
#                                  convert(Vector{Float64}, x), float(q), float(m10), float(x10), 
#                                  float(v10), float(v10), 0.0, 0.0, 
#                                  float(k), 0, float(Nosc), sf, float(cfl))

#     res = Simulation(sim, field_integrator.sol)
#     return res
# end

# ####
# function interacting_particles(L, N, tf, cfl, v0, q0, same_q=true, sf=true)
#     c = 1.0
#     equation = LinearScalarWaveEquation1D(c)

#     grid = UniformStaggeredGrid1D(Float64[-L, L], N)
#     x = GridFunctions.Grids.coords(grid) # cell centers

#     dx = spacing(grid)

#     dt = cfl * dx / c
#     @show dx, dt
#     t = 0.0:dt:(tf + dt)
#     tspan = (0.0, tf)
#     Φ = zeros(Float64, length(grid))
#     Π = zeros(Float64, length(grid))
#     Ψ = zeros(Float64, length(grid))

#     field_statevector = hcat(Π, Ψ, Φ)

#     ## Particle 1
#     q1 = Float64(q0)
#     m10 = 1.0
#     x10 = -200.0
#     v10 = v0
#     ## Particle 2
#     if same_q
#         q2 = q0
#     else
#         q2 = -q0
#     end
#     m20 = 1.0
#     x20 = -x10
#     v20 = -v10

#     # println("Particle 1:         Particle 2:")
#     # println("q = $q1             q = $q2")
#     # println("m = $m10            m = $m20")
#     # println("x = $x10            x = $x20")
#     # println("v = $v10            v = $v20")

#     particle_statevector = [m10, x10, v10, m20, x20, v20]
#     statevector = ArrayPartition(field_statevector,
#                                  particle_statevector)

#     numerical_flux = NumericalFluxes.flux_rusanov
#     slope_limiter = Limiters.minmod
#     flux_limiter = Limiters.minmod
#     interpolation_method = PCHIPInterpolation # QuadraticInterpolation

#     params = (equation=equation, numerical_flux=numerical_flux,
#               slope_limiter=slope_limiter, flux_limiter=flux_limiter,
#               h=dx, N=N, x=x, cfl=cfl, L=L,
#               dt=dt, x1=x10, x2=x20,
#               q1=q1, q2=q2, sf=sf,
#               interpolation_method=interpolation_method)#, acc=[1.0])

#     info_params = Dict(:L => L, :N => N, :dx => dx, :cfl => cfl, :dt => dt, :tf => tf,
#                        :q1 => q1, :m1 => m10, :x1 => x10, :v1 => v10,
#                        :q2 => q2, :m2 => m20, :x2 => x20, :v2 => v20,
#                        :SimType => "Interacting Particles")
#     print_simulation_info("Interacting Particles", info_params)

#     sim = SimulationParametersInteracting(L, N, dx, x,
#                                           q1, m10, x10, v10,
#                                           q2, m20, x20, v20, cfl)
#     alg = SSPRK54()
#     # alg = SSPRK53_2N1()
#     # alg = ImplicitEuler()
#     ode = ODEProblem{true}(ODE.interacting_coupled_rhs!, statevector, tspan,
#                            params; saveat=t,
#                            #callback=callbacks
#                            )

#     sol = solve(ode, alg; dt=dt)#, callback=cb)
#     res = Simulation(sim, sol)
#     return res
# end

# using ..FractionalStepMethods
# function forced_motion_fractional(Nosc, dx, cfl, vmax, q)
#     c = 1.0
#     equation = LinearScalarWaveEquation1D(c)

#     # vmax = 0.999
#     # ω = 2π
#     # T = 1.0
#     # A = vmax / ω

#     A = 1.0
#     ω = round(vmax / A; digits=9)
#     T = (2π / ω)

#     L = 2ceil(4Nosc * T / 3)
#     @show ω
#     @show A
#     @show L
#     @show T
#     N = 2Int64(L)
#     N = Int64(ceil(N / dx))
#     @show N
#     grid = UniformStaggeredGrid1D(Float64[-L, L], N)
#     x = coords(grid) # cell centers
#     @show length(x)
#     Φ = zeros(Float64, length(grid))
#     Π = zeros(Float64, length(grid))
#     Ψ = zeros(Float64, length(grid))

#     statevector = hcat(Π, Ψ, Φ)

#     dx = spacing(grid)

#     dt = cfl * dx / c
#     @show dx, dt

#     numerical_flux = NumericalFluxes.flux_rusanov
#     slope_limiter = Limiters.minmod
#     flux_limiter = Limiters.minmod

#     # q = +1
#     x10 = 0.0
#     direction = +1

#     tf = Nosc * T
#     L0 = 2tf
#     @show tf
#     @assert L > L0

#     tspan = (0.0, tf)
#     t = 0.0:dt:(tf + dt)
#     @show length(t)
#     sim = SimulationParameters1D(Int64(L), Int64(N), float(dx), float(dt), 
#                                  convert(Vector{Float64}, x), float(q), 1.0, float(x10), 
#                                  float(vmax), float(vmax), float(A), float(ω), 
#                                  0.0, direction, float(Nosc), true, float(cfl))

#     params = (equation=equation, numerical_flux=numerical_flux,
#               slope_limiter=slope_limiter, flux_limiter=flux_limiter,
#               h=dx, N=N, x=x, cfl=cfl, L=L, vmax=vmax, A=A,
#               dt=dt, x1=x10, ω=ω, q=q)

#     ode = SplitODEProblem(ODE.muscl!, ODE.source_rhs_forced!, statevector, tspan, params;
#                           saveat=t)

#     # if muscl
#     #     alg = SSPRK54()
#     #     # alg = SSPRK22()
#     #     if sf
#     #         ode = ODEProblem{true}(ODE.field_rhs_forced!,
#     #                                statevector, tspan, params;
#     #                                saveat=t)
#     #     else
#     #         ode = ODEProblem{true}(ODE.muscl!, statevector, tspan, params; saveat=t)
#     #     end
#     # else
#     #     alg = Euler()
#     #     ode = ODEProblem{true}(ODE.LxW!, statevector, tspan, params)
#     # end
#     alg = FractionalStepMethods.MyCustomSolver(Euler(), Euler())
#     # return statevector, params
#     sol = solve(ode, alg; adaptive=false, dt=dt)
#     res = Simulation(sim, sol)
#     return res
# end

end # end of module
