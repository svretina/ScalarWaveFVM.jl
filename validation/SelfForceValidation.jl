module SelfForceValidation

using ScalarWaveFVM
using OrdinaryDiffEq

include("RetardedFieldReference.jl")
using .RetardedFieldReference

export DAMPED_PARTICLE_CONFIG, run_damped_particle,
       particle_history, trajectory_history, turning_point_amplitudes,
       interpolate_at, reference_particle_trajectory,
       particle_validation_metrics, compare_remainder_field

const DAMPED_PARTICLE_CONFIG = (
    resolutions=(0.0500, 0.0250, 0.0125),
    validation_times=(2.0, 4.0, 6.0, 8.0, 9.0),
    n_oscillations=4.75,
    cfl=0.4,
    maximum_velocity=0.30,
    charge=0.35,
    amplitude=0.10,
    rest_mass=1.0,
)

function run_damped_particle(dx; config=DAMPED_PARTICLE_CONFIG)
    return ScalarWaveFVM.Run.particle_in_potential_relativistic(
        config.n_oscillations, dx, config.cfl,
        config.maximum_velocity, config.charge;
        amplitude=config.amplitude,
        initial_position=config.amplitude,
        initial_velocity=0.0,
        rest_mass=config.rest_mass,
    )
end

function particle_history(sim)
    particle = [state.x[2] for state in sim.sol.u]
    return (; t=collect(sim.sol.t),
            mass=getindex.(particle, 1),
            position=getindex.(particle, 2),
            velocity=getindex.(particle, 3),
            proper_time=getindex.(particle, 4))
end

function trajectory_history(sim)
    history = particle_history(sim)
    return TrajectoryHistory(history.t, history.position,
                             history.velocity, history.proper_time)
end

function turning_point_amplitudes(position)
    return [abs(position[i]) for i in 2:(length(position) - 1)
            if abs(position[i]) >= abs(position[i - 1]) &&
               abs(position[i]) > abs(position[i + 1])]
end

function interpolate_at(target_time, time, values)
    result = similar(target_time)
    for (i, target) in pairs(target_time)
        right = searchsortedfirst(time, target)
        if right <= 1
            result[i] = first(values)
        elseif right > length(time)
            result[i] = last(values)
        else
            left = right - 1
            weight = (target - time[left]) / (time[right] - time[left])
            result[i] = (1 - weight) * values[left] + weight * values[right]
        end
    end
    return result
end

"""Independent pure-retarded-self-field particle reference."""
function reference_particle_trajectory(final_time; config=DAMPED_PARTICLE_CONFIG)
    spring_constant = config.rest_mass *
                      (config.maximum_velocity / config.amplitude)^2

    function rhs!(du, u, p, t)
        mass, position, velocity, proper_time = u
        inverse_gamma = sqrt(1 - velocity^2)
        du[1] = -0.5 * config.charge^2 * inverse_gamma
        du[2] = velocity
        du[3] = -spring_constant * position * inverse_gamma^3 / mass
        du[4] = inverse_gamma
        return nothing
    end

    initial = [config.rest_mass, config.amplitude, 0.0, 0.0]
    problem = ODEProblem(rhs!, initial, (0.0, final_time))
    return solve(problem, Vern9(); abstol=1e-13, reltol=1e-13)
end

function particle_validation_metrics(sim, history, reference_solution;
                                     config=DAMPED_PARTICLE_CONFIG)
    spring_constant = config.rest_mass *
                      (config.maximum_velocity / config.amplitude)^2
    peaks = turning_point_amplitudes(history.position)
    reference = hcat(reference_solution.(history.t)...)'
    mass_identity = @. history.mass +
                       0.5 * config.charge^2 * history.proper_time
    energy_identity = @. history.mass / sqrt(1 - history.velocity^2) +
                         0.5 * spring_constant * history.position^2 +
                         0.5 * config.charge^2 * history.t

    regular_fields = [ScalarWaveFVM.ODE.regular_fields_at_particle(
                          state.x[1], state.x[2][2], sim.sol.prob.p)
                      for state in sim.sol.u]
    regular_phi = getindex.(regular_fields, 1)
    regular_pi = getindex.(regular_fields, 2)
    regular_psi = getindex.(regular_fields, 3)
    force_field = @. regular_psi + history.velocity * regular_pi
    mass_field = @. regular_pi + history.velocity * regular_psi

    return (minimum_mass=minimum(history.mass),
            maximum_speed=maximum(abs, history.velocity),
            first_peak=first(peaks), last_peak=last(peaks),
            fractional_peak_decay=1 - last(peaks) / first(peaks),
            maximum_position_error=maximum(abs,
                history.position .- reference[:, 2]),
            maximum_velocity_error=maximum(abs,
                history.velocity .- reference[:, 3]),
            maximum_mass_error=maximum(abs,
                history.mass .- reference[:, 1]),
            maximum_mass_residual=maximum(abs.(mass_identity .- config.rest_mass)),
            maximum_energy_residual=maximum(abs.(energy_identity .-
                                                  first(energy_identity))),
            maximum_regular_phi=maximum(abs, regular_phi),
            maximum_force_field=maximum(abs, force_field),
            maximum_mass_field=maximum(abs, mass_field))
end

function _state_at_saved_time(sim, time)
    index = findfirst(t -> isapprox(t, time; atol=32eps(time)), sim.sol.t)
    index === nothing && error("time $time is not present in the saved solution")
    return sim.sol.u[index]
end

function compare_remainder_field(sim, history, time;
                                 config=DAMPED_PARTICLE_CONFIG)
    params = sim.sol.prob.p
    h, x, L = params.h, params.x, params.L
    numerical = _state_at_saved_time(sim, time).x[1]
    initial = trajectory_state(history, first(history.t))

    boundary_safe = @. L - abs(x) > time + 4h
    front_safe = @. abs(abs(x - initial.z) - time) > 4h
    indices = findall(boundary_safe .& front_safe)
    reference = [cell_average_exact_remainder_1d(
                     history, time, x[i], h, config.charge)
                 for i in indices]
    reference_pi = getproperty.(reference, :Π)
    reference_psi = getproperty.(reference, :Ψ)
    reference_phi = getproperty.(reference, :Φ)
    error_pi = numerical[indices, 1] .- reference_pi
    error_psi = numerical[indices, 2] .- reference_psi
    error_phi = numerical[indices, 3] .- reference_phi

    metrics = (L1_pi=h * sum(abs, error_pi),
               L1_psi=h * sum(abs, error_psi),
               L1_phi=h * sum(abs, error_phi),
               Linf_pi=maximum(abs, error_pi),
               Linf_psi=maximum(abs, error_psi),
               Linf_phi=maximum(abs, error_phi),
               cells=length(indices))
    profile = (; x=x[indices],
               numerical_pi=numerical[indices, 1], reference_pi, error_pi,
               numerical_psi=numerical[indices, 2], reference_psi, error_psi,
               numerical_phi=numerical[indices, 3], reference_phi, error_phi)
    return metrics, profile
end

end
