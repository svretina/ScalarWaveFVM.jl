module PreparePlots

using ..ParticleMotion
using ..ForcedPlots
using ..InteractingPlots
using ..ParticlePotentialPlots
using ..PlottingUtils
using ..Run
using LaTeXStrings
using CairoMakie
using GLMakie
using LinearAlgebra
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
    plot_source_kink_harmonic(; paper=paper)
    plot_trajectory_makie(; paper=paper)
    # Reset path (optional, but good practice)
    PlottingUtils.reset_fig_path!()
    return nothing
end

function plot_source_kink_harmonic(; n_crossings=4, paper=false)
    CairoMakie.activate!()

    mytheme = mytheme_aps()
    set_theme!(mytheme)

    # --- Physics Parameters ---
    h = 1.0
    x_cell = 0.0
    q = 1.0

    # Motion Setup
    Ampl = 1.2
    ω = 0.1 / Ampl
    period = 2π / ω
    num_periods = ceil(n_crossings / 4)
    t_end = num_periods * period * 1.05
    dt = 0.005
    t_range = 0:dt:t_end
    nt = length(t_range)

    # --- Pre-calculate Physics ---
    source_val = zeros(nt)
    source_deriv = zeros(nt)
    # We only need z_hist to determine shading regions
    z_hist = zeros(nt)

    # Cell Geometry
    half_h = h / 2
    _h = 1.0 / h
    left_face = x_cell - half_h
    right_face = x_cell + half_h

    for (i, t) in enumerate(t_range)
        z = Ampl * sin(ω * t)
        v = Ampl * ω * cos(ω * t)
        acc = -Ampl * ω^2 * sin(ω * t)

        z_hist[i] = z

        # Logic
        base_val = 0.5 * q * acc
        val_L = -base_val / (1.0 + v)^2
        val_R = base_val / (1.0 - v)^2
        term_to_subtract = 0.0

        if z >= right_face
            term_to_subtract = val_L
        elseif z <= left_face
            term_to_subtract = val_R
        else
            frac_L = (z - left_face) * _h
            frac_R = 1.0 - frac_L
            term_to_subtract = frac_L * val_L + frac_R * val_R
        end

        source_val[i] = term_to_subtract
    end

    source_deriv[1:(end - 1)] = diff(source_val) ./ dt
    source_deriv[end] = source_deriv[end - 1]

    # --- CALCULATE SHADING INTERVALS ---
    # Identify time indices where particle is INSIDE the cell
    inside_mask = (z_hist .>= left_face) .& (z_hist .<= right_face)

    # Pad to detect edges correctly
    padded_mask = [false; inside_mask; false]
    changes = diff(padded_mask)
    starts = findall(x -> x == 1, changes)
    stops = findall(x -> x == -1, changes) .- 1

    # --- Visualization ---
    # Fig 1: Source Term
    fig_source = Figure(; size=(253, 200), figure_padding=(15, 15, 10, 10))
    ax_source = Axis(fig_source[1, 1]; xlabel=L"t", ylabel=L"S_{\Pi}(t)",
                     title=L"\textrm{Source Term}")

    # Fig 2: Time Derivative
    fig_deriv = Figure(; size=(253, 200), figure_padding=(15, 15, 10, 10))
    ax_deriv = Axis(fig_deriv[1, 1]; xlabel=L"t", ylabel=L"\partial_t S_{\Pi}",
                    title=L"\textrm{Time Derivative}")

    # 1. Add Shaded Regions (Background)
    for k in 1:length(starts)
        t_s = t_range[starts[k]]
        t_e = t_range[stops[k]]
        # Use gray shading for "Inside Cell" regions
        vspan!(ax_source, t_s, t_e; color=(:gray, 0.2))
        vspan!(ax_deriv, t_s, t_e; color=(:gray, 0.2))
    end

    # 2. Plot Curves
    lines!(ax_source, t_range, source_val; linewidth=1.7)
    lines!(ax_deriv, t_range, source_deriv; linewidth=1.7)

    # Save
    try
        savename_source = simple_name("source_kink_harmonic")
        savename_deriv = simple_name("source_kink_harmonic_deriv")

        # Save Source Term
        save(joinpath(get_fig_path(), savename_source), fig_source)
        # Save Derivative
        save(joinpath(get_fig_path(), savename_deriv), fig_deriv)

        if paper
            save(joinpath(get_paper_fig_path(), savename_source), fig_source)
            save(joinpath(get_paper_fig_path(), savename_deriv), fig_deriv)
        end
    catch
        @warn "Save paths not found."
    end
    return fig_source, fig_deriv
end

function get_line_points(P1, P2; len=100.0)
    center = P1
    dir = (P2 - P1) / norm(P2 - P1)
    start_pt = center - len * dir
    end_pt = center + len * dir
    return [start_pt[1], end_pt[1]], [start_pt[2], end_pt[2]]
end

function get_parabola_geometry(Focus, DirP1, DirP2)
    d_vec = DirP2 - DirP1
    n = [-d_vec[2], d_vec[1]] / norm(d_vec)
    dist = dot(n, Focus - DirP1)
    V = Focus - (dist / 2) * n
    p = abs(dist) / 2
    axis_dir = n * sign(dist)
    tangent_dir = [-axis_dir[2], axis_dir[1]]
    return V, p, axis_dir, tangent_dir
end

function get_parabola_pt_tangent(t, V, p, axis_dir, tangent_dir)
    pt = V + t .* tangent_dir + (t^2 / (4 * p)) .* axis_dir
    # Tangent vector: T(t) = tangent_dir + (2t / 4p) * axis_dir = tangent_dir + (t / 2p) * axis_dir
    tangent_vec = tangent_dir + (t / (2 * p)) .* axis_dir
    return pt, tangent_vec
end

function get_parabola_points(Focus, DirP1, DirP2; t_range)
    V, p, axis_dir, tangent_dir = get_parabola_geometry(Focus, DirP1, DirP2)
    xs, ys = Float64[], Float64[]
    for t in t_range
        pt, _ = get_parabola_pt_tangent(t, V, p, axis_dir, tangent_dir)
        push!(xs, pt[1])
        push!(ys, pt[2])
    end
    return xs, ys
end

function plot_trajectory_makie(; paper=false)
    CairoMakie.activate!()

    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Coordinates
    A = [0.0, 0.0]
    B = [-5.640519, 33.116689]
    C = [-39.379917, -10.392570]
    F = [-15.580886, -15.031356]
    G = [-16.320616, 0.287527]
    D = [0.0, 15.397437]
    E = [11.714639, 21.040898]

    # 1. Create Figure (Explicitly using CairoMakie.Figure)
    # Size: ~3.4 inches * 72 pt/inch ≈ 245 pt. Using 300 for safe margins.
    fig = Figure(; size=(253, 200), fontsize=12)

    # 2. Create Axis (Explicitly using CairoMakie.Axis)
    # DataAspect() ensures 1 unit x = 1 unit y (crucial for physics plots)
    ax = Axis(fig[1, 1]; aspect=DataAspect())

    # Hide decorations
    hidedecorations!(ax)
    hidespines!(ax)

    # -- Draw Trajectory (Thin) --
    V, p, axis_dir, tangent_dir = get_parabola_geometry(A, B, C)
    px, py = get_parabola_points(A, B, C; t_range=range(-65, 65; length=300))
    lines!(ax, px, py; color=:black, linewidth=1.0)

    # -- Calculate Points and Tangents for Velocity Vectors --
    # We choose t values that roughly correspond to the original F and D locations
    t_F = -25.0
    t_D = 10.0

    pt_F, tan_F = get_parabola_pt_tangent(t_F, V, p, axis_dir, tangent_dir)
    pt_D, tan_D = get_parabola_pt_tangent(t_D, V, p, axis_dir, tangent_dir)

    # Normalize tangents and scale for visualization if needed, 
    # but here we use them for the dashed line direction too.
    unit_tan_F = tan_F / norm(tan_F)
    unit_tan_D = tan_D / norm(tan_D)

    # -- Draw Gray Dashed Lines (Tangents) --
    len_dash = 50.0
    hx, hy = [pt_F[1] - len_dash * unit_tan_F[1], pt_F[1] + len_dash * unit_tan_F[1]],
             [pt_F[2] - len_dash * unit_tan_F[2], pt_F[2] + len_dash * unit_tan_F[2]]
    lines!(ax, hx, hy; color=:gray, linestyle=:dash, linewidth=0.8)

    gx, gy = [pt_D[1] - len_dash * unit_tan_D[1], pt_D[1] + len_dash * unit_tan_D[1]],
             [pt_D[2] - len_dash * unit_tan_D[2], pt_D[2] + len_dash * unit_tan_D[2]]
    lines!(ax, gx, gy; color=:gray, linestyle=:dash, linewidth=0.8)

    # -- Draw Velocity Vectors (Thick) --
    # Scale tangents for the arrow length
    v_scale = 15.0
    arrows!(ax, [pt_F[1]], [pt_F[2]], [v_scale * unit_tan_F[1]], [v_scale * unit_tan_F[2]];
            color=:black,
            linewidth=1.0,
            arrowsize=5)

    arrows!(ax, [pt_D[1]], [pt_D[2]], [v_scale * unit_tan_D[1]], [v_scale * unit_tan_D[2]];
            color=:black,
            linewidth=1.0,
            arrowsize=5)

    # -- Draw Points (Small) --
    scatter!(ax, [pt_F[1], pt_D[1]], [pt_F[2], pt_D[2]];
             color=:black,
             markersize=5)

    # -- Annotations --
    text!(ax, pt_F[1] - 7.5, pt_F[2] + 8; text=L"v_p^i(t_i\,)",
          align=(:center, :center), color=:black)

    text!(ax, pt_D[1] + 0, pt_D[2] + 8; text=L"v_p^{i}(t_{i+1})",
          align=(:center, :center), color=:black)

    # Set Limits
    xlims!(ax, -50, 20)
    ylims!(ax, -45, 45)

    # Save
    try
        savename = simple_name("trajectory_makie")
        save(joinpath(get_fig_path(), savename), fig)
        if paper
            save(joinpath(get_paper_fig_path(), savename), fig)
        end
    catch
        @warn "Save paths not found."
    end
    return fig
end

end # end of module
