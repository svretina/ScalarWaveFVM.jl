using ScalarWaveFVM
using CairoMakie
using DelimitedFiles

sim = ScalarWaveFVM.Run.particle_in_potential_relativistic(
    3.25, 0.02, 0.4, 0.3, 0.05;
    amplitude=0.1,
    initial_position=0.1,
    initial_velocity=0.0,
    rest_mass=1.0,
)

ScalarWaveFVM.enable_plotting!()
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
mkpath(output_dir)
stem = joinpath(output_dir, "relativistic_particle_history_3osc")
figure = ScalarWaveFVM.ParticleHistoryPlots.plot_particle_history(
    sim; savepath=stem * ".png")
save(stem * ".pdf", figure)

particle = [u.x[2] for u in sim.sol.u]
data = hcat(sim.sol.t,
            getindex.(particle, 1),
            getindex.(particle, 2),
            getindex.(particle, 3),
            getindex.(particle, 4))
open(stem * ".csv", "w") do io
    println(io, "t,mass,position,velocity,proper_time")
    writedlm(io, data, ',')
end

println("SAVED=$(stem).png")
println("SAVED=$(stem).pdf")
println("SAVED=$(stem).csv")
