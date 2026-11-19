using ScalarWaveFVM
using CairoMakie
using DelimitedFiles

include(joinpath(pkgdir(ScalarWaveFVM), "validation", "SelfForceValidation.jl"))
using .SelfForceValidation

config = DAMPED_PARTICLE_CONFIG
plot_only = "--plot-only" in ARGS
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
stem = joinpath(output_dir, "retarded_field_from_coupled_trajectory")
mkpath(output_dir)

if plot_only
    data = readdlm(stem * "_fine_profile.csv", ',', Float64; skipstart=1)
    final_profile = (x=data[:, 1],
                     numerical_pi=data[:, 2], reference_pi=data[:, 3],
                     error_pi=data[:, 4],
                     numerical_psi=data[:, 5], reference_psi=data[:, 6],
                     error_psi=data[:, 7],
                     numerical_phi=data[:, 8], reference_phi=data[:, 9],
                     error_phi=data[:, 10])
    summary = readdlm(stem * "_summary.csv", ',', Float64; skipstart=1)
else
    rows = Vector{Vector{Float64}}()
    final_profile = nothing
    for (resolution_index, dx) in enumerate(config.resolutions)
        sim = run_damped_particle(dx; config)
        history = trajectory_history(sim)
        for time in config.validation_times
            metrics, profile = compare_remainder_field(sim, history, time; config)
            push!(rows, [dx, sim.sol.prob.p.h, time, metrics.cells,
                         metrics.L1_pi, metrics.L1_psi, metrics.L1_phi,
                         metrics.Linf_pi, metrics.Linf_psi, metrics.Linf_phi])
            if resolution_index == length(config.resolutions) &&
               time == last(config.validation_times)
                global final_profile = profile
            end
        end
        sim = nothing
        history = nothing
        GC.gc()
    end
    summary = reduce(vcat, permutedims.(rows))
    open(stem * "_summary.csv", "w") do io
        println(io, "requested_dx,actual_dx,time,compared_cells," *
                    "L1_Pi,L1_Psi,L1_Phi,Linf_Pi,Linf_Psi,Linf_Phi")
        writedlm(io, summary, ',')
    end

    profile_data = hcat(final_profile.x,
                        final_profile.numerical_pi, final_profile.reference_pi,
                        final_profile.error_pi,
                        final_profile.numerical_psi, final_profile.reference_psi,
                        final_profile.error_psi,
                        final_profile.numerical_phi, final_profile.reference_phi,
                        final_profile.error_phi)
    open(stem * "_fine_profile.csv", "w") do io
        println(io, "x,numerical_Pi,reference_Pi,error_Pi," *
                    "numerical_Psi,reference_Psi,error_Psi," *
                    "numerical_Phi,reference_Phi,error_Phi")
        writedlm(io, profile_data, ',')
    end
end

ScalarWaveFVM.enable_plotting!()
figure = ScalarWaveFVM.RetardedPotentialPlots.plot_remainder_reference_comparison(
    final_profile; savepath=stem * ".png")
save(stem * ".pdf", figure)

println("SAVED=$(stem).png")
println("SUMMARY_ROWS=")
display(summary)
