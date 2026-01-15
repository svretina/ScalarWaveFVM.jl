module PreparePlots

using ..ParticleMotion
using ..ForcedPlots
using ..InteractingPlots
using ..ParticlePotentialPlots
using ..PlottingUtils
using ..Run

using JLD2

function run_radiation_sims(Nosc, dx0, cfl, vmax, q)
    sim1 = Run.particle_in_potential(Nosc, dx0, cfl, vmax, q)
    sim2 = Run.particle_in_potential(Nosc, dx0 / 2, cfl, vmax, q)
    sim4 = Run.particle_in_potential(Nosc, dx0 / 4, cfl, vmax, q)
    return sim1, sim2, sim4
end

function run_forced_sims(Nosc, dx0, cfl, vmax, q)
    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q)
    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q)
    return sim1, sim2, sim4
end

function prepare_forced_plots(sim1, sim2, sim4; paper=true)
    tf = sim1.sol.t[end]
    target_time = round(tf; digits=1)
    tidx_end = argmin(abs.(sim1.sol.t .- target_time))
    println("Plotting at t = $(sim1.sol.t[tidx_end]) (target: $target_time, index: $tidx_end)")

    tidx_end2 = 2tidx_end - 1
    tidx_end4 = 4tidx_end - 3

    @assert sim1.sol.t[tidx_end] == sim2.sol.t[tidx_end2]
    @assert sim1.sol.t[tidx_end] == sim4.sol.t[tidx_end4]

    Nosc = sim1.params.Nosc
    dx0 = sim1.params.dx
    cfl = sim1.params.cfl
    vmax = sim1.params.vmax
    q = sim1.params.q

    ## 1) Πr + Retarded field
    ForcedPlots.plot_pifield_resolutions_vs_retarded(tidx_end2, sim2, sim4; paper=paper)
    ## 2) Ψr + Retarded field
    ForcedPlots.plot_psifield_resolutions_vs_retarded(tidx_end2, sim2, sim4; paper=paper)
    ## 3) Convergence of Πr vs Retarded field
    ForcedPlots.plot_pifield_error_convergence(tidx_end, Nosc, dx0, cfl, vmax, q;
                                               paper=paper)
    ## 4) Convergence of Ψr vs Retarded field
    ForcedPlots.plot_psifield_error_convergence(tidx_end, Nosc, dx0, cfl, vmax, q;
                                                paper=paper)
    ## 5) Self Convergence Πr
    ForcedPlots.plot_pifield_self_convergence_order(sim1, sim2, sim4; paper=paper)
    ## 6) Self Convergence Ψr
    ForcedPlots.plot_psifield_self_convergence_order(sim1, sim2, sim4; paper=paper)
    return nothing
end

function prepare_potential_plots(sim1, sim2, sim4; paper=true)
    tf = sim1.sol.t[end]
    target_time = round(tf; digits=1)
    tidx_end = argmin(abs.(sim1.sol.t .- target_time))
    println("Plotting at t = $(sim1.sol.t[tidx_end]) (target: $target_time, index: $tidx_end)")

    tidx_end2 = 2tidx_end - 1
    tidx_end4 = 4tidx_end - 3

    @assert sim1.sol.t[tidx_end] == sim2.sol.t[tidx_end2]
    @assert sim1.sol.t[tidx_end] == sim4.sol.t[tidx_end4]

    ## 1) Πr Resolutions
    ParticlePotentialPlots.plot_pifield_resolutions(tidx_end, sim1, sim2, sim4;
                                                    paper=paper)
    ## 2) Ψr Resolutions
    ParticlePotentialPlots.plot_psifield_resolutions(tidx_end, sim1, sim2, sim4;
                                                     paper=paper)
    ## 3) Self Convergence Πr
    ParticlePotentialPlots.plot_pifield_self_convergence_order(sim1, sim2, sim4;
                                                               paper=paper)
    ## 4) Self Convergence Ψr
    ParticlePotentialPlots.plot_psifield_self_convergence_order(sim1, sim2, sim4;
                                                                paper=paper)
    ## 5) Position
    ParticlePotentialPlots.plot_particle_position(sim1; paper=paper)
    ## 6) Velocity
    ParticlePotentialPlots.plot_particle_velocity(sim1; paper=paper)
    ## 7) Position Resolutions
    ParticlePotentialPlots.plot_particle_position_resolutions(sim1, sim2, sim4;
                                                              paper=paper)
    ## 8) Velocity Resolutions
    ParticlePotentialPlots.plot_particle_velocity_resolutions(sim1, sim2, sim4;
                                                              paper=paper)
    ## 9) Position Self Convergence
    ParticlePotentialPlots.plot_particle_position_self_convergence_order(sim1, sim2, sim4;
                                                                         paper=paper)
    ## 10) Velocity Self Convergence
    ParticlePotentialPlots.plot_particle_velocity_self_convergence_order(sim1, sim2, sim4;
                                                                         paper=paper)
    ## 11) Energy
    ParticlePotentialPlots.plot_energy(sim4; paper=paper)
    ## 12) Scaled Differences (Position)
    ParticlePotentialPlots.plot_particle_position_diff_scaled(sim1, sim2, sim4; paper=paper)
    ## 13) Scaled Differences (Velocity)
    ParticlePotentialPlots.plot_particle_velocity_diff_scaled(sim1, sim2, sim4; paper=paper)
    return nothing
end

function prepare_all_plots(; paper=true)
    ## forced motion case
    println("\n" * "="^50)
    println("PREPARING PLOTS: Forced Motion")
    println("="^50)
    q = 1.0
    vmax = 0.1
    cfl = 0.1
    Nosc = 4

    dx = 0.1 / 2
    # Create experiment folder
    PlottingUtils.set_experiment_path!("forced_motion", q, vmax, Nosc, dx, cfl)

    sim1, sim2, sim4 = run_forced_sims(Nosc, dx, cfl, vmax, q)
    prepare_forced_plots(sim1, sim2, sim4; paper=paper)

    ## particle in potential
    println("\n" * "="^50)
    println("PREPARING PLOTS: Particle in harmonic potential")
    println("="^50)
    q = 1
    # Create experiment folder
    PlottingUtils.set_experiment_path!("radiation_reaction", q, vmax, Nosc, dx, cfl)

    sim1, sim2, sim4 = run_radiation_sims(Nosc, dx, cfl, vmax, q)
    prepare_potential_plots(sim1, sim2, sim4; paper=paper)

    # Reset path (optional, but good practice)
    PlottingUtils.reset_fig_path!()
    return nothing
end

end # end of module
