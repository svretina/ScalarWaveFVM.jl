using ScalarWaveFVM
using CairoMakie
using DelimitedFiles
using LaTeXStrings

include(joinpath(pkgdir(ScalarWaveFVM), "validation", "SelfForceValidation.jl"))
using .SelfForceValidation

config = DAMPED_PARTICLE_CONFIG
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
mkpath(output_dir)
stem = joinpath(output_dir, "damped_trajectory_resolution_comparison")

angular_frequency = config.maximum_velocity / config.amplitude
final_time = config.n_oscillations * 2π / angular_frequency
reference_solution = reference_particle_trajectory(final_time; config)
histories = NamedTuple[]
validation = NamedTuple[]

for dx in config.resolutions
    sim = run_damped_particle(dx; config)
    history = particle_history(sim)
    metrics = particle_validation_metrics(sim, history, reference_solution; config)
    push!(histories, history)
    push!(validation, metrics)

    history_data = hcat(history.t, history.mass, history.position,
                        history.velocity, history.proper_time)
    filename = stem * "_dx_$(replace(string(dx), '.' => 'p')).csv"
    open(filename, "w") do io
        println(io, "t,mass,position,velocity,proper_time")
        writedlm(io, history_data, ',')
    end

    sim = nothing
    GC.gc()
end

reference_time = range(0.0, minimum(last(history.t) for history in histories);
                       length=4001)
reference_states = reference_solution.(reference_time)
reference = (; t=reference_time, position=getindex.(reference_states, 2))

ScalarWaveFVM.enable_plotting!()
labels = [L"\Delta x=0.05", L"\Delta x=0.025", L"\Delta x=0.0125"]
figure = ScalarWaveFVM.ParticleHistoryPlots.plot_particle_trajectory_resolutions(
    histories, labels; reference, savepath=stem * ".png")
save(stem * ".pdf", figure)

finest = last(histories)
summary = zeros(length(histories), 17)
for (i, (dx, history, metrics)) in
    enumerate(zip(config.resolutions, histories, validation))
    finest_position = interpolate_at(history.t, finest.t, finest.position)
    maximum_difference = maximum(abs, history.position .- finest_position)
    summary[i, :] .= (dx, dx, config.cfl * dx,
                      metrics.minimum_mass, metrics.maximum_speed,
                      metrics.first_peak, metrics.last_peak,
                      metrics.fractional_peak_decay, maximum_difference,
                      metrics.maximum_position_error,
                      metrics.maximum_velocity_error,
                      metrics.maximum_mass_error,
                      metrics.maximum_mass_residual,
                      metrics.maximum_energy_residual,
                      metrics.maximum_regular_phi,
                      metrics.maximum_force_field,
                      metrics.maximum_mass_field)
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
