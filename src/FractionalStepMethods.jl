module FractionalStepMethods

using UnPack
using DiffEqBase           # Core interfaces and problem types
using OrdinaryDiffEqSSPRK  # SSPRK solvers for flux (hyperbolic)
using OrdinaryDiffEq       # Implicit solvers for source, includes solve and ODEProblem
import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, alg_cache  # Explicitly import the type

# Define splitting method
@enum SplittingMethod LieTrotter Strang Yoshida

# Custom solver struct, subtyping OrdinaryDiffEqAlgorithm
struct MyCustomSolver{S1,S2} <: OrdinaryDiffEqAlgorithm
    solver_flux::S1    # Solver for flux term (L)
    solver_source::S2  # Solver for source term (S)
    splitting::SplittingMethod
end

# Constructor with default splitting
function MyCustomSolver(solver_flux, solver_source; splitting=Strang)
    MyCustomSolver(solver_flux, solver_source, splitting)
end

# Define a minimal cache for MyCustomSolver
struct MyCustomSolverCache{C1,C2}
    cache_flux::C1    # Cache for flux solver
    cache_source::C2  # Cache for source solver
end

# Implement alg_cache for MyCustomSolver
function OrdinaryDiffEq.alg_cache(alg::MyCustomSolver, u, rate_prototype,
                                  ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
                                  ::Type{tTypeNoUnits}, uprev, f::SciMLBase.SplitFunction,
                                  t, dt, reltol, p, calck,
                                  ::Val{inplace}) where {uEltypeNoUnits,
                                                         uBottomEltypeNoUnits,tTypeNoUnits,
                                                         inplace}
    # Delegate to sub-solvers' alg_cache, using f.f1 and f.f2
    cache_flux = alg_cache(alg.solver_flux, u, rate_prototype, uEltypeNoUnits,
                           uBottomEltypeNoUnits, tTypeNoUnits, uprev, f.f1, t, dt, reltol,
                           p, calck, Val(inplace))
    cache_source = alg_cache(alg.solver_source, u, rate_prototype, uEltypeNoUnits,
                             uBottomEltypeNoUnits, tTypeNoUnits, uprev, f.f2, t, dt, reltol,
                             p, calck, Val(inplace))
    MyCustomSolverCache(cache_flux, cache_source)
end

# Extend SciML solve interface, allocating u and du
function DiffEqBase.__solve(prob::SplitODEProblem, alg::MyCustomSolver; dt=nothing,
                            saveat=nothing, kwargs...)
    @unpack f1, f2, u0, tspan, p = prob  # f1 = flux_rhs, f2 = source_rhs
    t_start, t_end = tspan

    # Determine time points to save
    if saveat === nothing
        if dt === nothing
            error("Provide dt or saveat for consistent stepping")
        end
        N = Int(ceil((t_end - t_start) / dt))
        t = range(t_start; stop=t_end, length=N + 1)
    else
        t = saveat isa Number ? collect(t_start:saveat:t_end) : collect(saveat)
        N = length(t) - 1
        dt = dt === nothing ? (t[2] - t[1]) : dt  # Use dt as hint if provided
    end

    # Allocate u and du once, efficiently
    u = [i == 1 ? copy(u0) : similar(u0) for i in 1:length(t)]  # Solution history
    du = similar(u0)  # Working array for RHS and intermediate states

    # Check if functions are in-place
    inplace = DiffEqBase.isinplace(prob)
    if !inplace
        error("This solver requires in-place RHS functions for non-allocating behavior")
    end

    # Select splitting function
    if alg.splitting == Strang
        step_func = strang_step_inplace
    elseif alg.splitting == Yoshida
        w1_yoshida = 1 / (2 - 2^(1 / 3))  # ≈ 0.6756035959798289
        w2_yoshida = (1 - 2 * w1_yoshida) / 2  # ≈ -0.3512017919596578
        step_func = (u, du, f1, f2, p, t, dt, s_flux, s_source) -> yoshida_step_inplace(u,
                                                                                        du,
                                                                                        f1,
                                                                                        f2,
                                                                                        p,
                                                                                        t,
                                                                                        dt,
                                                                                        s_flux,
                                                                                        s_source,
                                                                                        w1_yoshida,
                                                                                        w2_yoshida)
    else  # LieTrotter
        step_func = lie_trotter_step_inplace
    end

    # Time-stepping loop (no allocations here)
    for i in 1:N
        step_func(u[i], du, f1, f2, p, t[i], t[i + 1] - t[i], alg.solver_flux,
                  alg.solver_source)
        u[i + 1] .= du  # Update next state from du
    end

    # Build and return solution
    return DiffEqBase.build_solution(prob, alg, t, u; retcode=:Success)
end

# Strang splitting (in-place, non-allocating)
function strang_step_inplace(u, du, f1, f2, p, t, dt, solver_flux, solver_source)
    # Step 1: Source for dt/2
    prob_s1 = ODEProblem{true}(f2, u, (t, t + dt / 2), p)
    integrator_s1 = init(prob_s1, solver_source; dt=dt / 2)  # dt provided
    integrator_s1.u .= u
    step!(integrator_s1, dt / 2, true)
    du .= integrator_s1.u

    # Step 2: Flux for dt
    prob_f = ODEProblem{true}(f1, du, (t, t + dt), p)
    integrator_f = init(prob_f, solver_flux; dt=dt)  # dt provided
    integrator_f.u .= du
    step!(integrator_f, dt, true)
    du .= integrator_f.u

    # Step 3: Source for dt/2
    prob_s2 = ODEProblem{true}(f2, du, (t + dt / 2, t + dt), p)
    integrator_s2 = init(prob_s2, solver_source; dt=dt / 2)  # dt provided
    integrator_s2.u .= du
    step!(integrator_s2, dt / 2, true)
    du .= integrator_s2.u  # Final state in du
end

# Lie-Trotter splitting (in-place, non-allocating)
function lie_trotter_step_inplace(u, du, f1, f2, p, t, dt, solver_flux, solver_source)
    # Step 1: Flux for dt
    prob_f = ODEProblem{true}(f1, u, (t, t + dt), p)
    integrator_f = init(prob_f, solver_flux; dt=dt)  # dt provided
    integrator_f.u .= u
    step!(integrator_f, dt, true)
    du .= integrator_f.u

    # Step 2: Source for dt
    prob_s = ODEProblem{true}(f2, du, (t, t + dt), p)
    integrator_s = init(prob_s, solver_source; dt=dt)  # dt provided
    integrator_s.u .= du
    step!(integrator_s, dt, true)
    du .= integrator_s.u  # Final state in du
end

# Yoshida splitting (in-place, non-allocating)
function yoshida_step_inplace(u, du, f1, f2, p, t, dt, solver_flux, solver_source,
                              w1_yoshida, w2_yoshida)
    t0 = t
    # Step 1: S for w1 * dt
    prob_s1 = ODEProblem{true}(f2, u, (t0, t0 + w1_yoshida * dt), p)
    integrator_s1 = init(prob_s1, solver_source; dt=w1_yoshida * dt)
    integrator_s1.u .= u
    step!(integrator_s1, w1_yoshida * dt, true)
    du .= integrator_s1.u
    t0 += w1_yoshida * dt

    # Step 2: L for w1 * dt
    prob_l1 = ODEProblem{true}(f1, du, (t0, t0 + w1_yoshida * dt), p)
    integrator_l1 = init(prob_l1, solver_flux; dt=w1_yoshida * dt)
    integrator_l1.u .= du
    step!(integrator_l1, w1_yoshida * dt, true)
    du .= integrator_l1.u
    t0 += w1_yoshida * dt

    # Step 3: S for w2 * dt
    prob_s2 = ODEProblem{true}(f2, du, (t0, t0 + w2_yoshida * dt), p)
    integrator_s2 = init(prob_s2, solver_source; dt=w2_yoshida * dt)
    integrator_s2.u .= du
    step!(integrator_s2, w2_yoshida * dt, true)
    du .= integrator_s2.u
    t0 += w2_yoshida * dt

    # Step 4: L for w2 * dt
    prob_l2 = ODEProblem{true}(f1, du, (t0, t0 + w2_yoshida * dt), p)
    integrator_l2 = init(prob_l2, solver_flux; dt=w2_yoshida * dt)
    integrator_l2.u .= du
    step!(integrator_l2, w2_yoshida * dt, true)
    du .= integrator_l2.u
    t0 += w2_yoshida * dt

    # Step 5: S for w2 * dt
    prob_s3 = ODEProblem{true}(f2, du, (t0, t0 + w2_yoshida * dt), p)
    integrator_s3 = init(prob_s3, solver_source; dt=w2_yoshida * dt)
    integrator_s3.u .= du
    step!(integrator_s3, w2_yoshida * dt, true)
    du .= integrator_s3.u
    t0 += w2_yoshida * dt

    # Step 6: L for w1 * dt
    prob_l3 = ODEProblem{true}(f1, du, (t0, t0 + w1_yoshida * dt), p)
    integrator_l3 = init(prob_l3, solver_flux; dt=w1_yoshida * dt)
    integrator_l3.u .= du
    step!(integrator_l3, w1_yoshida * dt, true)
    du .= integrator_l3.u
    t0 += w1_yoshida * dt

    # Step 7: S for w1 * dt
    prob_s4 = ODEProblem{true}(f2, du, (t0, t0 + w1_yoshida * dt), p)
    integrator_s4 = init(prob_s4, solver_source; dt=w1_yoshida * dt)
    integrator_s4.u .= du
    step!(integrator_s4, w1_yoshida * dt, true)
    du .= integrator_s4.u  # Final state in du
end

# sol = solve(split_ode, MyCustomSolver(Tsit5(), RK4(), splitting=Strang), dt=dt, saveat=0.1)

end # end of module
