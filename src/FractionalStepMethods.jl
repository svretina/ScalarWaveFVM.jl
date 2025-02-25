module FractionalStepMethods

using DiffEqBase

# Define splitting method
@enum SplittingMethod LieTrotter Strang

# Custom solver struct
struct MyCustomSolver{S1,S2}
    solver_flux::S1    # Solver for flux term (L)
    solver_source::S2  # Solver for source term (S)
    splitting::SplittingMethod
end

# Constructor with default splitting
function MyCustomSolver(solver_flux, solver_source; splitting=Strang)
    MyCustomSolver(solver_flux, solver_source, splitting)
end

# Extend SciML solve interface for SplitODEProblem
function DiffEqBase.__solve(prob::SplitODEProblem, alg::MyCustomSolver; dt, saveat=nothing,
                            kwargs...)
    @unpack f1, f2, u0, tspan, p = prob  # f1 = flux_rhs, f2 = source_rhs
    t_start, t_end = tspan
    N = Int(ceil((t_end - t_start) / dt))  # Number of full steps

    # Determine time points to save (passed to solve, not ODEProblem)
    if saveat === nothing
        t = range(t_start; stop=t_end, length=N + 1)
    else
        t = saveat isa Number ? collect(t_start:saveat:t_end) : collect(saveat)
        N = length(t) - 1  # Adjust N based on saveat
        dt = t[2] - t[1]   # Assume uniform saveat for simplicity; adjust if needed
    end

    # Preallocate solution (u0 can be a matrix)
    u = Vector{typeof(u0)}(undef, length(t))
    u[1] = copy(u0)  # Works for scalars, vectors, or matrices

    # Check if functions are in-place
    inplace = DiffEqBase.isinplace(prob)
    step_func = alg.splitting == Strang ?
                (inplace ? strang_step_inplace : strang_step_outofplace) :
                (inplace ? lie_trotter_step_inplace : lie_trotter_step_outofplace)

    # Time-stepping loop
    for i in 1:N
        u[i + 1] = step_func(u[i], f1, f2, p, t[i], dt, alg.solver_flux, alg.solver_source)
    end

    # Build and return solution
    return DiffEqBase.build_solution(prob, alg, t, u; retcode=:Success)
end

# Strang splitting (in-place)
function strang_step_inplace(u, f1, f2, p, t, dt, solver_flux, solver_source)
    # Step 1: Source for dt/2
    prob_s1 = ODEProblem{true}(f2, u, (t, t + dt / 2), p)
    u_ss = solve(prob_s1, solver_source; dt=dt / 2, saveat=[t + dt / 2]).u[end]

    # Step 2: Flux for dt
    prob_f = ODEProblem{true}(f1, u_ss, (t, t + dt), p)
    u_s = solve(prob_f, solver_flux; dt=dt, saveat=[t + dt]).u[end]

    # Step 3: Source for dt/2
    prob_s2 = ODEProblem{true}(f2, u_s, (t + dt / 2, t + dt), p)
    u_next = solve(prob_s2, solver_source; dt=dt / 2, saveat=[t + dt]).u[end]

    return u_next
end

# Strang splitting (out-of-place)
function strang_step_outofplace(u, f1, f2, p, t, dt, solver_flux, solver_source)
    prob_s1 = ODEProblem{false}(f2, u, (t, t + dt / 2), p)
    u_ss = solve(prob_s1, solver_source; dt=dt / 2, saveat=[t + dt / 2]).u[end]

    prob_f = ODEProblem{false}(f1, u_ss, (t, t + dt), p)
    u_s = solve(prob_f, solver_flux; dt=dt, saveat=[t + dt]).u[end]

    prob_s2 = ODEProblem{false}(f2, u_s, (t + dt / 2, t + dt), p)
    u_next = solve(prob_s2, solver_source; dt=dt / 2, saveat=[t + dt]).u[end]

    return u_next
end

# Lie-Trotter splitting (in-place)
function lie_trotter_step_inplace(u, f1, f2, p, t, dt, solver_flux, solver_source)
    prob_f = ODEProblem{true}(f1, u, (t, t + dt), p)
    u_f = solve(prob_f, solver_flux; dt=dt, saveat=[t + dt]).u[end]

    prob_s = ODEProblem{true}(f2, u_f, (t, t + dt), p)
    u_next = solve(prob_s, solver_source; dt=dt, saveat=[t + dt]).u[end]

    return u_next
end

# Lie-Trotter splitting (out-of-place)
function lie_trotter_step_outofplace(u, f1, f2, p, t, dt, solver_flux, solver_source)
    prob_f = ODEProblem{false}(f1, u, (t, t + dt), p)
    u_f = solve(prob_f, solver_flux; dt=dt, saveat=[t + dt]).u[end]

    prob_s = ODEProblem{false}(f2, u_f, (t, t + dt), p)
    u_next = solve(prob_s, solver_source; dt=dt, saveat=[t + dt]).u[end]

    return u_next
end

# Example usage with a matrix u0
function flux_rhs!(du, u, p, t)  # Example: simple wave flux
    c, dx = p
    n = size(u, 1)
    du[:, 1] .= c^2 .* (circshift(u[:, 2], -1) .- circshift(u[:, 2], 1)) ./ (2 * dx)  # v' = c^2 w_x
    du[:, 2] .= (circshift(u[:, 1], -1) .- circshift(u[:, 1], 1)) ./ (2 * dx)         # w' = v_x
end

function source_rhs!(du, u, p, t)  # Example: discontinuous source
    du[:, 1] .= t > 0.5 ? 1.0 : 0.0  # Step source on v
    du[:, 2] .= 0.0                  # No source on w
end

# sol = solve(split_ode, MyCustomSolver(Tsit5(), RK4(), splitting=Strang), dt=dt, saveat=0.1)

end # end of module
