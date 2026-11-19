using ScalarWaveFVM
using CairoMakie
using DelimitedFiles

output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
particle_path = joinpath(output_dir,
    "damped_trajectory_resolution_comparison_summary.csv")
field_path = joinpath(output_dir,
    "retarded_field_from_coupled_trajectory_summary.csv")
stem = joinpath(output_dir, "end_to_end_self_force_validation")

particle = readdlm(particle_path, ',', Float64; skipstart=1)
field_all = readdlm(field_path, ',', Float64; skipstart=1)
field = field_all[field_all[:, 3] .== maximum(field_all[:, 3]), :]

h = particle[:, 2]
field_h = field[:, 2]
field_errors = (field[:, 5], field[:, 6], field[:, 7])
force_error = particle[:, 16]
particle_errors = (particle[:, 10], particle[:, 11], particle[:, 12])

strictly_decreases(values) = all(diff(values) .< 0)
all(strictly_decreases, field_errors) ||
    error("retarded-field errors do not all decrease under refinement")
strictly_decreases(force_error) ||
    error("worldline spatial-force residual does not decrease under refinement")
all(strictly_decreases, particle_errors) ||
    error("reduced-particle errors do not all decrease under refinement")

orders(values) = log2.(values[1:(end - 1)] ./ values[2:end])
metrics = ("field_Pi_L1", "field_Psi_L1", "field_Phi_L1",
           "worldline_force", "particle_z", "particle_v", "particle_m")
values = (field_errors..., force_error, particle_errors...)
open(stem * "_summary.csv", "w") do io
    println(io, "metric,coarse,medium,fine,order_coarse_medium,order_medium_fine")
    for (metric, error_values) in zip(metrics, values)
        convergence = orders(error_values)
        println(io, join((metric, error_values..., convergence...), ','))
    end
end

ScalarWaveFVM.enable_plotting!()
figure = ScalarWaveFVM.ParticleHistoryPlots.plot_self_force_validation(
    field_h, field_errors, h, force_error, particle_errors;
    savepath=stem * ".png")
save(stem * ".pdf", figure)

println("PASS: retarded fields, spatial self-force, and reduced dynamics converge")
println("SAVED=$(stem).png")
