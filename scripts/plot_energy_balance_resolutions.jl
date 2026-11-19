using ScalarWaveFVM
using CairoMakie
using DelimitedFiles
using LaTeXStrings

include(joinpath(pkgdir(ScalarWaveFVM), "validation", "SelfForceValidation.jl"))
using .SelfForceValidation

config = DAMPED_PARTICLE_CONFIG
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
stem = joinpath(output_dir, "remainder_energy_balance_resolutions")
mkpath(output_dir)

balances = NamedTuple[]
summary = zeros(length(config.resolutions), 10)
for (i, dx) in enumerate(config.resolutions)
    sim = run_damped_particle(dx; config)
    balance = ScalarWaveFVM.Diagnostics.remainder_energy_balance(sim)
    push!(balances, balance)

    energy_scale = maximum(abs, balance.energy_change)
    maximum_residual = maximum(abs, balance.residual)
    summary[i, :] .= (dx, sim.sol.prob.p.h, sim.sol.prob.p.dt,
                      maximum_residual, maximum_residual / energy_scale,
                      last(balance.energy_change), last(balance.predicted_change),
                      last(balance.integrated_source),
                      -last(balance.integrated_dissipation),
                      last(balance.integrated_reconstruction))

    data = hcat(balance.t, balance.energy, balance.energy_change,
                balance.integrated_source, balance.integrated_dissipation,
                balance.integrated_reconstruction, balance.predicted_change,
                balance.residual)
    filename = stem * "_dx_$(replace(string(dx), '.' => 'p')).csv"
    open(filename, "w") do io
        println(io, "t,energy,energy_change,integrated_source," *
                    "integrated_dissipation,integrated_reconstruction," *
                    "predicted_change,residual")
        writedlm(io, data, ',')
    end

    sim = nothing
    GC.gc()
end

open(stem * "_summary.csv", "w") do io
    println(io, "requested_dx,actual_dx,dt,max_abs_residual," *
                "relative_max_residual,final_energy_change," *
                "final_predicted_change,final_integrated_source," *
                "final_negative_integrated_dissipation," *
                "final_integrated_reconstruction")
    writedlm(io, summary, ',')
end

ScalarWaveFVM.enable_plotting!()
labels = [L"\Delta x=0.05", L"\Delta x=0.025", L"\Delta x=0.0125"]
figure = ScalarWaveFVM.ParticleHistoryPlots.plot_remainder_energy_balance_resolutions(
    balances, labels; savepath=stem * ".png")
save(stem * ".pdf", figure)

println("SAVED=$(stem).png")
println("SUMMARY_ROWS=")
display(summary)
