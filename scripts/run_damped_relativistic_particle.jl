using ScalarWaveFVM
using CairoMakie
using DelimitedFiles

sim = ScalarWaveFVM.Run.particle_in_potential_relativistic(
    6.25, 0.025, 0.4, 0.30, 0.20;
    amplitude=0.10,
    initial_position=0.10,
    initial_velocity=0.0,
    rest_mass=1.0,
)

ScalarWaveFVM.enable_plotting!()
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
mkpath(output_dir)
stem = joinpath(output_dir, "damped_relativistic_particle")

history = ScalarWaveFVM.ParticleHistoryPlots.plot_particle_history(
    sim; savepath=stem * "_history.png")
save(stem * "_history.pdf", history)
balance_figure = ScalarWaveFVM.ParticleHistoryPlots.plot_remainder_energy_balance(
    sim; savepath=stem * "_energy_balance.png")
save(stem * "_energy_balance.pdf", balance_figure)

particle = [u.x[2] for u in sim.sol.u]
mass = getindex.(particle, 1)
position = getindex.(particle, 2)
velocity = getindex.(particle, 3)
data = hcat(sim.sol.t,
            mass,
            position,
            velocity,
            getindex.(particle, 4))
open(stem * "_history.csv", "w") do io
    println(io, "t,mass,position,velocity,proper_time")
    writedlm(io, data, ',')
end

balance = ScalarWaveFVM.Diagnostics.remainder_energy_balance(sim)
balance_data = hcat(balance.t, balance.energy, balance.energy_change,
                    balance.integrated_source, balance.integrated_dissipation,
                    balance.integrated_reconstruction, balance.predicted_change,
                    balance.residual)
open(stem * "_energy_balance.csv", "w") do io
    println(io, "t,energy,energy_change,integrated_source,integrated_dissipation," *
                "integrated_reconstruction,predicted_change,residual")
    writedlm(io, balance_data, ',')
end

println("SAVED=$(stem)_history.png")
println("SAVED=$(stem)_energy_balance.png")
peak_amplitudes = [abs(position[i]) for i in 2:(length(position) - 1)
                   if abs(position[i]) >= abs(position[i - 1]) &&
                      abs(position[i]) > abs(position[i + 1])]
println("MASS_RANGE=$(extrema(mass))")
println("VELOCITY_RANGE=$(extrema(velocity))")
println("POSITION_PEAKS=$(peak_amplitudes)")
println("MAX_BALANCE_RESIDUAL=$(maximum(abs, balance.residual))")
