module PreparePlots

using ..ForcedPlots
using ..InteractingPlots
using ..ParticlePotentialPlots
using ..Run

function run_potential_sims(L, N, tf, cfl, q, sf)
    sim1 = Run.particle_in_potential(L, N, tf, cfl, q, sf)
    sim2 = Run.particle_in_potential(L, 2N, tf, cfl, q, sf)
    sim4 = Run.particle_in_potential(L, 4N, tf, cfl, q, sf)
    return sim1, sim2, sim4
end

function run_forced_sims(Nosc, dx0, cfl, vmax, q, sf=true)
    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q, sf)
    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q, sf)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q, sf)
    return sim1, sim2, sim4
end

function run_interacting_sims(L, N, tf, cfl, v0, q0; same_q=true, sf=true)
    sim1 = Run.interacting_particles(L, N, tf, cfl, v0, q0, same_q, sf)
    sim2 = Run.interacting_particles(L, 2N, tf, cfl, v0, q0, same_q, sf)
    sim4 = Run.interacting_particles(L, 4N, tf, cfl, v0, q0, same_q, sf)
    return sim1, sim2, sim4
end

function prepare_forced_plots(Nosc, dx0, cfl, vmax, q, sf=true)
    sim1, sim2, sim4 = run_forced_sims(Nosc, dx0, cfl, vmax, q, sf)
    tidx_end = length(sim1.sol.t) - 1
    _ = ForcedPlots.plot_pifield_resolutions(tidx_end, sim1, sim2, sim4)
    _ = ForcedPlots.plot_psifield_resolutions(tidx_end, sim1, sim2, sim4)
    _ = ForcedPlots.plot_pifield_convord(sim1, sim2, sim4)
    _ = ForcedPlots.plot_psifield_convord(sim1, sim2, sim4)
    return nothing
end

function prepare_interacting_plots(L, N, tf, cfl, v0, q0; same_q=true, sf=true)
    sim1, sim2, sim4 = run_interacting_sims(L, N, tf, cfl, v0, q0; same_q=same_q, sf=sf)
    tidx_end = length(sim1.sol.t) - 1
    _ = InteractingPlots.plot_pifield_resolutions(tidx_end, sim1, sim2, sim4)
    _ = InteractingPlots.plot_psifield_resolutions(tidx_end, sim1, sim2, sim4)
    _ = InteractingPlots.plot_pifield_convord(sim1, sim2, sim4)
    _ = InteractingPlots.plot_psifield_convord(sim1, sim2, sim4)
    _ = InteractingPlots.plot_particle_positions_resolutions(sim1, sim2, sim4)
    _ = InteractingPlots.plot_particle_velocities_resolutions(sim1, sim2, sim4)
    _ = InteractingPlots.plot_particle_masses_resolutions(sim1, sim2, sim4)
    return nothing
end

function prepare_potential_plots(L, N, cfl, q)
    tf = 0.75L
    sim1sf, sim2sf, sim4sf = run_potential_sims(L, N, tf, cfl, q, true)
    sim1, sim2, sim4 = run_potential_sims(L, N, tf, cfl, q, false)

    tidx_end = length(sim1.sol.t) - 1
    _ = ParticlePotentialPlots.plot_pifield_resolutions(tidx_end, sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_psifield_resolutions(tidx_end, sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_pifield_convord(sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_psifield_convord(sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_particle_positions_resolutions(sim1, sim2, sim4,
                                                                   sim1sf, sim2sf, sim4sf)
    return nothing
end

end # end of module
