using ScalarWaveFVM
using CairoMakie
using DelimitedFiles
using LaTeXStrings

include(joinpath(pkgdir(ScalarWaveFVM), "validation", "SelfForceValidation.jl"))
using .SelfForceValidation

config = DAMPED_PARTICLE_CONFIG
control_volume = (-1.0, 1.0)
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
stem = joinpath(output_dir, "remainder_energy_balance_resolutions")
mkpath(output_dir)

balances = NamedTuple[]
summary = zeros(length(config.resolutions), 13)
for (i, dx) in enumerate(config.resolutions)
    sim = run_damped_particle(dx; config)
    balance = ScalarWaveFVM.Diagnostics.remainder_energy_balance(
        sim; domain=control_volume)
    push!(balances, balance)

    energy_scale = maximum(abs, balance.energy_change)
    maximum_residual = maximum(abs, balance.residual)
    summary[i, :] .= (dx, sim.sol.prob.p.h, sim.sol.prob.p.dt,
                      maximum_residual, maximum_residual / energy_scale,
                      last(balance.energy_change), last(balance.predicted_change),
                      last(balance.integrated_source),
                      -last(balance.integrated_boundary_flux_left),
                      -last(balance.integrated_boundary_flux_right),
                      -last(balance.integrated_boundary_flux),
                      -last(balance.integrated_dissipation),
                      last(balance.integrated_reconstruction))

    data = hcat(balance.t, balance.energy, balance.energy_change,
                balance.integrated_source, balance.integrated_boundary_flux_left,
                balance.integrated_boundary_flux_right,
                balance.integrated_boundary_flux, balance.integrated_dissipation,
                balance.integrated_reconstruction, balance.predicted_change,
                balance.residual, balance.energy_rate, balance.source_power,
                balance.boundary_flux_left, balance.boundary_flux_right,
                balance.dissipation, balance.reconstruction,
                balance.predicted_rate, balance.power_residual)
    filename = stem * "_dx_$(replace(string(dx), '.' => 'p')).csv"
    open(filename, "w") do io
        println(io, "t,energy,energy_change,integrated_source," *
                    "integrated_boundary_flux_left," *
                    "integrated_boundary_flux_right,integrated_boundary_flux," *
                    "integrated_dissipation,integrated_reconstruction," *
                    "predicted_change,residual,energy_rate,source_power," *
                    "boundary_flux_left,boundary_flux_right,dissipation," *
                    "reconstruction,predicted_rate,power_residual")
        writedlm(io, data, ',')
    end

    sim = nothing
    GC.gc()
end

open(stem * "_summary.csv", "w") do io
    println(io, "requested_dx,actual_dx,dt,max_abs_residual," *
                "relative_max_residual,final_energy_change," *
                "final_predicted_change,final_integrated_source," *
                "final_negative_integrated_boundary_flux_left," *
                "final_negative_integrated_boundary_flux_right," *
                "final_negative_integrated_boundary_flux," *
                "final_negative_integrated_dissipation," *
                "final_integrated_reconstruction")
    writedlm(io, summary, ',')
end

ScalarWaveFVM.enable_plotting!()
labels = [L"\Delta x=0.05", L"\Delta x=0.025", L"\Delta x=0.0125"]
figure = ScalarWaveFVM.ParticleHistoryPlots.plot_remainder_energy_balance_resolutions(
    balances, labels; savepath=stem * ".png")
save(stem * ".pdf", figure)

power_stem = joinpath(output_dir, "remainder_energy_power_resolutions")
power_figure = ScalarWaveFVM.ParticleHistoryPlots.plot_remainder_energy_power_resolutions(
    balances, labels; savepath=power_stem * ".png")
save(power_stem * ".pdf", power_figure)

println("SAVED=$(stem).png")
println("SAVED=$(power_stem).png")
println("SUMMARY_ROWS=")
display(summary)
