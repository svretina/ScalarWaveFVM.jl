module PreparePlots

using ..ForcedPlots
using ..InteractingPlots
using ..ParticlePotentialPlots
using ..Run

using JLD2

function run_potential_sims(L, N, tf, cfl, q, sf)
    sim1 = Run.particle_in_potential(L, N, tf, cfl, q, sf)
    sim2 = Run.particle_in_potential(L, 2N, tf, cfl, q, sf)
    sim4 = Run.particle_in_potential(L, 4N, tf, cfl, q, sf)
    return sim1, sim2, sim4
end

function run_fixed_potential_sims(L, dx, dt, Nosc, k)
    sim1 = Run.particle_in_fixed_potential_interpolation(L, dx, dt, Nosc, k)
    sim2 = Run.particle_in_fixed_potential_interpolation(L, dx, dt / 2, Nosc, k)
    sim4 = Run.particle_in_fixed_potential_interpolation(L, dx, dt / 4, Nosc, k)
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

function prepare_forced_plots(Nosc, dx0, cfl, vmax, q; paper=true, sf=true)
    sim1, sim2, sim4 = run_forced_sims(Nosc, dx0, cfl, vmax, q, sf)
    tidx_end = length(sim1.sol.t)
    _ = ForcedPlots.plot_pifield_resolutions(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_pifield_resolutions(Int64(floor(tidx_end/2)), sim1, sim2, sim4, 4, paper)
    _ = ForcedPlots.plot_pifield_resolutions(Int64(floor(tidx_end/4)), sim1, sim2, sim4, 8, paper)
    
    _ = ForcedPlots.plot_psifield_resolutions(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_psifield_resolutions(Int64(floor(tidx_end/2)), sim1, sim2, sim4,4, paper)
    _ = ForcedPlots.plot_psifield_resolutions(Int64(floor(tidx_end/4)), sim1, sim2, sim4,8, paper)

    _ = ForcedPlots.plot_pifield_diff(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_pifield_diff(Int64(floor(tidx_end/2)), sim1, sim2, sim4, 4, paper)
    _ = ForcedPlots.plot_pifield_diff(Int64(floor(tidx_end/4)), sim1, sim2, sim4, 8, paper)
    
    _ = ForcedPlots.plot_psifield_diff(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_psifield_diff(Int64(floor(tidx_end/2)), sim1, sim2, sim4, 4, paper)
    _ = ForcedPlots.plot_psifield_diff(Int64(floor(tidx_end/4)), sim1, sim2, sim4, 8, paper)
  
    _ = ForcedPlots.plot_pifield_convord(sim1, sim2, sim4, paper)
    _ = ForcedPlots.plot_psifield_convord(sim1, sim2, sim4, paper)
    return nothing
end

function prepare_forced_plots(sim1, sim2, sim4; paper=true, sf=true)    
    tidx_end = length(sim1.sol.t)
    _ = ForcedPlots.plot_pifield_resolutions(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_pifield_resolutions(Int64(floor(tidx_end/2)), sim1, sim2, sim4, 4, paper)
    _ = ForcedPlots.plot_pifield_resolutions(Int64(floor(tidx_end/4)), sim1, sim2, sim4, 8, paper)
    
    _ = ForcedPlots.plot_psifield_resolutions(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_psifield_resolutions(Int64(floor(tidx_end/2)), sim1, sim2, sim4, 4, paper)
    _ = ForcedPlots.plot_psifield_resolutions(Int64(floor(tidx_end/4)), sim1, sim2, sim4, 8, paper)

    _ = ForcedPlots.plot_pifield_diff(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_pifield_diff(Int64(floor(tidx_end/2)), sim1, sim2, sim4, 4, paper)
    _ = ForcedPlots.plot_pifield_diff(Int64(floor(tidx_end/4)), sim1, sim2, sim4, 8, paper)
    
    _ = ForcedPlots.plot_psifield_diff(tidx_end, sim1, sim2, sim4, 2, paper)
    _ = ForcedPlots.plot_psifield_diff(Int64(floor(tidx_end/2)), sim1, sim2, sim4, 4, paper)
    _ = ForcedPlots.plot_psifield_diff(Int64(floor(tidx_end/4)), sim1, sim2, sim4, 8, paper)
  
    
    _ = ForcedPlots.plot_pifield_convord(sim1, sim2, sim4, paper)
    _ = ForcedPlots.plot_psifield_convord(sim1, sim2, sim4, paper)
    return nothing
end



function prepare_interacting_plots(L, N, tf, cfl, v0, q0; same_q=true, sf=true, paper=true)
    sim1, sim2, sim4 = run_interacting_sims(L, N, tf, cfl, v0, q0; same_q=same_q, sf=sf, paper=paper)
    tidx_end = length(sim1.sol.t)
    _ = InteractingPlots.plot_pifield_resolutions(tidx_end, sim1, sim2, sim4, paper=paper)
    _ = InteractingPlots.plot_psifield_resolutions(tidx_end, sim1, sim2, sim4, paper=paper)
    _ = InteractingPlots.plot_pifield_convord(sim1, sim2, sim4, paper=paper)
    _ = InteractingPlots.plot_psifield_convord(sim1, sim2, sim4, paper=paper)
    _ = InteractingPlots.plot_particle_positions_resolutions(sim1, sim2, sim4, paper=paper)
    _ = InteractingPlots.plot_particle_velocities_resolutions(sim1, sim2, sim4, paper=paper)
    _ = InteractingPlots.plot_particle_masses_resolutions(sim1, sim2, sim4, paper=paper)
    return nothing
end

function prepare_potential_plots(L, N, cfl, q)
    tf = 0.75L
    sim1sf, sim2sf, sim4sf = run_potential_sims(L, N, tf, cfl, q, true)
    sim1, sim2, sim4 = run_potential_sims(L, N, tf, cfl, q, false)

    tidx_end = length(sim1.sol.t)
    _ = ParticlePotentialPlots.plot_pifield_resolutions(tidx_end, sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_psifield_resolutions(tidx_end, sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_pifield_convord(sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_psifield_convord(sim1sf, sim2sf, sim4sf)
    _ = ParticlePotentialPlots.plot_particle_positions_resolutions(sim1, sim2, sim4,
                                                                   sim1sf, sim2sf, sim4sf)
    return nothing
end

function prepare_fixed_potential_plots(L, dx, dt, Nosc, k; paper=true)
    tf = 0.75L
    sim1, sim2, sim4 = run_fixed_potential_sims(L, dx, dt, Nosc, k)

    tidx_end = length(sim1.sol.t)
    _ = ParticlePotentialPlots.plot_particle_position_analytic(sim1; invert=false, paper=paper)
    _ = ParticlePotentialPlots.plot_particle_velocity_analytic(sim1; paper=paper)
    _ = ParticlePotentialPlots.plot_particle_mass_analytic(sim1; paper=paper)
    _ = ParticlePotentialPlots.plot_particle_position_analytic_diff(sim1, sim2, sim4; paper=paper)
    _ = ParticlePotentialPlots.plot_particle_velocity_analytic_diff(sim1, sim2, sim4; paper=paper)
    _ = ParticlePotentialPlots.plot_particle_mass_analytic_diff(sim1, sim2, sim4; paper=paper)

    return nothing
end

function prepare_all_plots(;paper=true)
    ## forced motion case
    println("Running Forced Motion")
    q=1.0
    vmax=0.5
    cfl=0.25
    Nosc = 3
    
    dx0=0.5
    prepare_forced_plots(Nosc, dx0, cfl, vmax, q; paper=paper, sf=true)
    dx0=0.25
    prepare_forced_plots(Nosc, dx0, cfl, vmax, q; paper=paper, sf=true)
    dx0=0.125/4
    prepare_forced_plots(Nosc, dx0, cfl, vmax, q; paper=paper, sf=true)

    ## interacting particles
    L = 2000
    N = 2000
    cfl=0.25
    q=0.1
    
    tf=270
    v0=0.5
    prepare_interacting_plots(L, N, tf, cfl, v0, q0; same_q=true, sf=true, paper=paper)
    tf=100
    v0=0.0
    prepare_interacting_plots(L, N, tf, cfl, v0, q0; same_q=true, sf=true, paper=paper)
    tf=130
    v0=0.0
    prepare_interacting_plots(L, N, tf, cfl, v0, q0; same_q=false, sf=true, paper=paper)

    ## particle in potential
    println("Running Particle in harmonic oscillator potential")
    println("k=1")
    k = 1
    dt=0.1
    q=1
    L=100
    dx=0.5
    Nosc=3
    prepare_fixed_potential_plots(L,dx,dt,Nosc, k; paper=paper)

    println("k=50")
    k = 50
    dt=0.001
    q=1
    L=100
    dx=0.1
    Nosc=3
    prepare_fixed_potential_plots(L,dx,dt,Nosc, k; paper=paper)
    return nothing
end

end # end of module
