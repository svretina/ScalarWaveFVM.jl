using ScalarWaveFVM
using CairoMakie
using DelimitedFiles
using LaTeXStrings
using OrdinaryDiffEq

const RESOLUTIONS = (0.0500, 0.0250, 0.0125)
const N_OSCILLATIONS = 4.75
const CFL = 0.4
const V_MAX = 0.30
const CHARGE = 0.35
const AMPLITUDE = 0.10
const REST_MASS = 1.0
const SPRING_CONSTANT = REST_MASS * (V_MAX / AMPLITUDE)^2

function particle_history(sim)
    particle = [u.x[2] for u in sim.sol.u]
    return (; t=sim.sol.t,
            mass=getindex.(particle, 1),
            position=getindex.(particle, 2),
            velocity=getindex.(particle, 3),
            proper_time=getindex.(particle, 4))
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
function reference_particle_trajectory(final_time)
    function rhs!(du, u, p, t)
        mass, position, velocity, proper_time = u
        inverse_gamma = sqrt(1 - velocity^2)
        du[1] = -0.5 * CHARGE^2 * inverse_gamma
        du[2] = velocity
        du[3] = -SPRING_CONSTANT * position * inverse_gamma^3 / mass
        du[4] = inverse_gamma
        return nothing
    end

    problem = ODEProblem(rhs!, [REST_MASS, AMPLITUDE, 0.0, 0.0],
                         (0.0, final_time))
    return solve(problem, Vern9(); abstol=1e-13, reltol=1e-13)
end

sims = [ScalarWaveFVM.Run.particle_in_potential_relativistic(
            N_OSCILLATIONS, dx, CFL, V_MAX, CHARGE;
            amplitude=AMPLITUDE,
            initial_position=AMPLITUDE,
            initial_velocity=0.0,
            rest_mass=REST_MASS,
        ) for dx in RESOLUTIONS]
histories = particle_history.(sims)
reference_solution = reference_particle_trajectory(last(first(histories).t))
reference_time = range(0.0, last(first(histories).t); length=4001)
reference_states = reference_solution.(reference_time)
reference = (; t=reference_time,
             position=getindex.(reference_states, 2))

ScalarWaveFVM.enable_plotting!()
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
mkpath(output_dir)
stem = joinpath(output_dir, "damped_trajectory_resolution_comparison")
labels = [L"\Delta x=0.05", L"\Delta x=0.025", L"\Delta x=0.0125"]
figure = ScalarWaveFVM.ParticleHistoryPlots.plot_particle_trajectory_resolutions(
    sims, labels; reference, savepath=stem * ".png")
save(stem * ".pdf", figure)

finest = last(histories)
summary = zeros(length(sims), 17)
for (i, (dx, sim, history)) in enumerate(zip(RESOLUTIONS, sims, histories))
    peaks = turning_point_amplitudes(history.position)
    finest_position = interpolate_at(history.t, finest.t, finest.position)
    maximum_difference = maximum(abs, history.position .- finest_position)
    reference_values = hcat(reference_solution.(history.t)...)'
    mass_identity = @. history.mass + 0.5 * CHARGE^2 * history.proper_time
    energy_identity = @. history.mass / sqrt(1 - history.velocity^2) +
                         0.5 * SPRING_CONSTANT * history.position^2 +
                         0.5 * CHARGE^2 * history.t

    regular_fields = [ScalarWaveFVM.ODE.regular_fields_at_particle(
                          state.x[1], state.x[2][2], sim.sol.prob.p)
                      for state in sim.sol.u]
    regular_phi = getindex.(regular_fields, 1)
    regular_pi = getindex.(regular_fields, 2)
    regular_psi = getindex.(regular_fields, 3)
    force_field = @. regular_psi + history.velocity * regular_pi
    mass_field = @. regular_pi + history.velocity * regular_psi

    summary[i, :] .= (dx, sim.sol.prob.p.h, sim.sol.prob.p.dt,
                      minimum(history.mass), maximum(abs, history.velocity),
                      first(peaks), last(peaks),
                      1 - last(peaks) / first(peaks), maximum_difference,
                      maximum(abs, history.position .- reference_values[:, 2]),
                      maximum(abs, history.velocity .- reference_values[:, 3]),
                      maximum(abs, history.mass .- reference_values[:, 1]),
                      maximum(abs.(mass_identity .- REST_MASS)),
                      maximum(abs.(energy_identity .- first(energy_identity))),
                      maximum(abs, regular_phi), maximum(abs, force_field),
                      maximum(abs, mass_field))

    history_data = hcat(history.t, history.mass, history.position,
                        history.velocity, history.proper_time)
    filename = stem * "_dx_$(replace(string(dx), '.' => 'p')).csv"
    open(filename, "w") do io
        println(io, "t,mass,position,velocity,proper_time")
        writedlm(io, history_data, ',')
    end
end

open(stem * "_summary.csv", "w") do io
    println(io, "requested_dx,actual_dx,dt,min_mass,max_abs_velocity," *
                "first_peak,last_peak,fractional_peak_decay," *
                "max_position_difference_vs_finest,max_position_error_vs_reference," *
                "max_velocity_error_vs_reference,max_mass_error_vs_reference," *
                "max_mass_identity_residual,max_energy_identity_residual," *
                "max_abs_regular_phi,max_abs_force_field,max_abs_mass_field")
    writedlm(io, summary, ',')
end

println("SAVED=$(stem).png")
println("SUMMARY_ROWS=")
display(summary)
