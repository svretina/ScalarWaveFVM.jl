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
    plot_source_kink_harmonic()
    # Reset path (optional, but good practice)
    PlottingUtils.reset_fig_path!()
    return nothing
end

"""
    plot_source_kink_demonstration(; paper=false)

Plots the time-evolution of the 'Exact Integration' source term for a single cell 
as a particle traverses it. Visualizes the C0 'kink' and the C-1 derivative discontinuity
that degrades RK45 convergence.
"""
function plot_source_kink_const_vel(; paper=false)
    CairoMakie.activate!()
    try
        set_theme!(mytheme_aps())
    catch
        @warn "mytheme_aps() not found, using default theme."
        set_theme!(theme_latexfonts()) # Fallback close to APS style
    end

    # --- Simulation Parameters ---
    h = 1.0           # Grid spacing
    x_cell = 0.0      # Cell centered at 0 (Faces at -0.5, 0.5)
    v = 0.5           # Particle velocity
    q = 1.0           # Charge
    a = 1.0           # Amplitude factor

    # Time setup: Particle moves from z = -1.0 to z = 1.0
    # Crosses Left Face (-0.5) at t=1.0, Right Face (0.5) at t=3.0
    z_start = -1.0
    dt = 0.001
    t_end = 4.0
    t_range = 0:dt:t_end

    # Pre-allocate arrays
    nt = length(t_range)
    source_val = zeros(nt)
    source_deriv = zeros(nt)

    # --- Physics Logic (Replicating your _calculate_singular_sources) ---
    # We track the source term *for the cell at x_cell* as z changes
    for (i, t) in enumerate(t_range)
        z = z_start + v * t

        # Exact copy of your logic
        left_face = x_cell - h / 2
        right_face = x_cell + h / 2

        Δr = right_face - z
        Δl = left_face - z

        # The KINK generators: abs() and sign()
        absΔr = abs(Δr)
        absΔl = abs(Δl)
        σr = sign(Δr)
        σl = sign(Δl)

        den_r = (1.0 - v * σr)^2
        den_l = (1.0 - v * σl)^2

        frac_r = absΔr / den_r
        frac_l = absΔl / den_l

        # We plot term_Π (Momentum Source) as the example
        term_Π = -0.5 * q * (frac_r - frac_l)

        source_val[i] = a * term_Π / h
    end

    # Calculate Finite Difference Derivative (what the RK solver "feels")
    source_deriv[1:(end - 1)] = diff(source_val) ./ dt
    source_deriv[end] = source_deriv[end - 1] # Pad last

    # --- Plotting ---
    fig = Figure(; size=(2 * 253, 300), figure_padding=(15, 15, 10, 10)) # APS column width approx

    # Top Axis: The Source Term
    ax1 = Axis(fig[1, 1];
               xlabel=L"t", # Shared x-axis, label hidden
               ylabel=L"S_{\Pi}(t)",
               title=L"\textrm{Exact Source Term (Piecewise } {C}^0)")

    # Bottom Axis: The Time Derivative
    ax2 = Axis(fig[1, 2];
               xlabel=L"t",
               ylabel=L"\partial_t S_{\Pi}",
               title=L"\textrm{Time Derivative (Discontinuous)}",
               limits=(nothing, (-1.5, 1.5)))

    # Link x-axes
    linkxaxes!(ax1, ax2)

    # Plot Data
    lines!(ax1, t_range, source_val; color=:blue, linewidth=1.5, label=L"Source")
    lines!(ax2, t_range, source_deriv; color=:red, linewidth=1.5,
           label=L"\partial_t Source")

    # Mark Grid Crossings (The "Kinks")
    # Face 1 crossing: -0.5 - (-1.0) / 0.5 = 1.0s
    # Face 2 crossing:  0.5 - (-1.0) / 0.5 = 3.0s
    # vlines!(ax1, [1.0, 3.0]; color=:black, linestyle=:dash, linewidth=0.8, alpha=0.5)
    # vlines!(ax2, [1.0, 3.0]; color=:black, linestyle=:dash, linewidth=0.8, alpha=0.5)

    colgap!(fig.layout, 15)
    # Layout adjustments
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    # colsize!(fig.layout, 1, Relative(0.9))

    # Save logic matching your function
    try
        savename = simple_name("source_kink_demonstration")
        save(joinpath(get_fig_path(), savename), fig)
        if paper
            save(joinpath(get_paper_fig_path(), savename), fig)
        end
    catch
        @warn "Save paths/functions not defined. Displaying figure only."
    end

    return fig
end

function plot_source_kink_const_accel(; paper=false)
    CairoMakie.activate!()
    try
        set_theme!(mytheme_aps())
    catch
        @warn "Theme not found, using fallback."
        set_theme!(theme_latexfonts())
    end

    # --- Simulation Parameters ---
    h = 1.0
    x_cell = 0.0
    q = 1.0
    a = 1.0

    # Kinematics: Safe Acceleration (Subsonic)
    # We want to cross the cell (-0.5 to 0.5) but keep v < 1.0
    z_start = -1.0
    v0 = 0.4        # Start with some velocity
    acc = 0.1       # Gentle acceleration
    # At t=4, v = 0.4 + 0.1*4 = 0.8 (Safe from singularity)

    # Time setup
    dt = 0.001
    t_end = 4.0
    t_range = 0:dt:t_end

    nt = length(t_range)
    source_val = zeros(nt)
    source_deriv = zeros(nt)

    # --- Physics Logic ---
    for (i, t) in enumerate(t_range)
        # Kinematics: z = z0 + v0*t + 0.5*a*t^2
        z = z_start + v0 * t + 0.5 * acc * t^2
        v = v0 + acc * t

        left_face = x_cell - h / 2
        right_face = x_cell + h / 2

        Δr = right_face - z
        Δl = left_face - z

        absΔr = abs(Δr)
        absΔl = abs(Δl)
        σr = sign(Δr)
        σl = sign(Δl)

        # Denominator (Doppler factor)
        den_r = (1.0 - v * σr)^2
        den_l = (1.0 - v * σl)^2

        frac_r = absΔr / den_r
        frac_l = absΔl / den_l

        term_Π = -0.5 * q * (frac_r - frac_l)
        source_val[i] = a * term_Π / h
    end

    # Finite Difference
    source_deriv[1:(end - 1)] = diff(source_val) ./ dt
    source_deriv[end] = source_deriv[end - 1]

    # --- Plotting ---
    fig = Figure(; size=(506, 250), figure_padding=(15, 15, 10, 10))

    # Left Axis: Source Term
    ax1 = Axis(fig[1, 1];
               xlabel=L"t",
               ylabel=L"S_{\Pi}(t)",
               title=L"\textrm{Accelerating Source } (C^0)",
               limits=((-0.2, 4.2), nothing)) # Let Y auto-scale

    # Right Axis: Derivative
    # We auto-scale Y here too because the jump height varies with velocity
    ax2 = Axis(fig[1, 2];
               xlabel=L"t",
               ylabel=L"\partial_t S_{\Pi}",
               title=L"\textrm{Time Derivative (Discontinuous)}",
               limits=((-0.2, 4.2), nothing))

    linkxaxes!(ax1, ax2)

    lines!(ax1, t_range, source_val; color=:blue, linewidth=1.5)
    lines!(ax2, t_range, source_deriv; color=:red, linewidth=1.5)

    # Calculate exact crossing times for vertical lines
    # Quadratic equation: 0.5*acc*t^2 + v0*t + (z_start - boundary) = 0
    function solve_crossing(target_z)
        c = z_start - target_z
        discriminant = v0^2 - 4 * (0.5 * acc) * c
        if discriminant >= 0
            return (-v0 + sqrt(discriminant)) / acc
        end
        return NaN
    end

    t1 = solve_crossing(-0.5) # Enter cell
    t2 = solve_crossing(0.5)  # Exit cell

    if !isnan(t1)
        vlines!(ax1, [t1]; color=:darkgray, linestyle=:dash, linewidth=0.8)
        vlines!(ax2, [t1]; color=:darkgray, linestyle=:dash, linewidth=0.8)
    end
    if !isnan(t2)
        vlines!(ax1, [t2]; color=:darkgray, linestyle=:dash, linewidth=0.8)
        vlines!(ax2, [t2]; color=:darkgray, linestyle=:dash, linewidth=0.8)
    end

    # --- Layout ---
    colsize!(fig.layout, 1, Relative(0.5))
    colsize!(fig.layout, 2, Relative(0.5))
    colgap!(fig.layout, 25)

    try
        savename = simple_name("source_kink_accel")
        save(joinpath(get_fig_path(), savename), fig)
        if paper
            save(joinpath(get_paper_fig_path(), savename), fig)
        end
    catch
        @warn "Save paths not found."
    end

    return fig
end

function interactive_source_field_kink_accel(; n_crossings=4)
    GLMakie.activate!()

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
    z_hist = zeros(nt)
    v_hist = zeros(nt)

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
        v_hist[i] = v

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
    # Find all time intervals where particle is INSIDE the cell
    inside_mask = (z_hist .>= left_face) .& (z_hist .<= right_face)

    # We detect changes in the mask to find start/stop indices
    # Pad with false to handle edge cases (starting or ending inside)
    padded_mask = [false; inside_mask; false]
    changes = diff(padded_mask)
    starts = findall(x -> x == 1, changes) # Index in padded, which maps to t_range index
    stops = findall(x -> x == -1, changes) .- 1 # Adjust for padding shift

    # --- Visualization ---
    fig = Figure(; size=(900, 900), fontsize=18, figure_padding=30)

    # Slider
    sg = SliderGrid(fig[3, 1],
                    (label="Time", range=t_range, format="{:.3f} s", startvalue=0.0))
    time_slider = sg.sliders[1]
    idx = lift(t -> searchsortedfirst(t_range, t), time_slider.value)

    # Observables
    z_cur = lift(i -> z_hist[i], idx)
    v_cur = lift(i -> v_hist[i], idx)

    # Top Plot: Field Context
    x_space = range(-1.8, 1.8; length=800)
    field_profile = lift((z, v) -> begin
                             map(x -> abs(x - z) / (1.0 - v * sign(x - z))^2, x_space)
                         end, z_cur, v_cur)

    pt_source = lift(i -> Point2f(t_range[i], source_val[i]), idx)
    pt_deriv = lift(i -> Point2f(t_range[i], source_deriv[i]), idx)

    # Plot 1: Spatial Domain
    ax_field = Axis(fig[1, 1];
                    title="Physical Context", xlabel="Position (x)",
                    ylabel="Field Strength",
                    limits=(-1.8, 1.8, -0.5, 15.0))

    ax_motion = Axis(fig[1, 1];
                     yaxisposition=:right, ylabel="Particle Track",
                     yticklabelcolor=:gray, ylabelcolor=:gray,
                     limits=(-1.8, 1.8, -1, 1), yticks=[0])
    hidespines!(ax_motion)
    hidexdecorations!(ax_motion)

    # Shade the cell in the spatial plot
    vspan!(ax_field, -0.5, 0.5; color=(:gray, 0.2))
    vlines!(ax_field, [-0.5, 0.5]; color=:black, linestyle=:dash, linewidth=1)

    lines!(ax_field, x_space, field_profile; color=:blue, linewidth=3)
    scatter!(ax_motion, lift(z -> Point2f(z, 0), z_cur); color=:red, markersize=18)

    # --- PLOTS 2 & 3: Source & Derivative ---
    g_bottom = fig[2, 1] = GridLayout()
    colgap!(g_bottom, 30)

    ax_source = Axis(g_bottom[1, 1]; xlabel="t", ylabel="Source Term Value",
                     title="Source Term")
    ax_deriv = Axis(g_bottom[1, 2]; xlabel="t", ylabel="Time Derivative",
                    title="Derivative")

    # 1. ADD SHADED REGIONS (Background)
    # Loop through calculated intervals and add gray vspans
    for k in 1:length(starts)
        t_s = t_range[starts[k]]
        t_e = t_range[stops[k]]
        vspan!(ax_source, t_s, t_e; color=(:gray, 0.2)) # Light gray area
        vspan!(ax_deriv, t_s, t_e; color=(:gray, 0.2))
    end

    # 2. Plot Standard Solid Curves
    lines!(ax_source, t_range, source_val; color=:black, linewidth=2)
    lines!(ax_deriv, t_range, source_deriv; color=:red, linewidth=2)

    # Indicators
    vlines!(ax_source, lift(i -> t_range[i], idx); color=:gray)
    vlines!(ax_deriv, lift(i -> t_range[i], idx); color=:gray)

    scatter!(ax_source, pt_source; color=:gold, markersize=15, strokewidth=1,
             strokecolor=:black)
    scatter!(ax_deriv, pt_deriv; color=:red, markersize=15, strokewidth=1,
             strokecolor=:black)

    # Layout Priority
    rowsize!(fig.layout, 1, Relative(0.6))

    display(fig)
    return fig
end

# function save_source_kink_gif(; n_crossings=4, filename="source_kink.gif")
#     GLMakie.activate!()

#     # --- Physics Parameters ---
#     h = 1.0
#     x_cell = 0.0
#     q = 1.0

#     # Motion Setup: z(t) = Ampl * sin(ωt)
#     Ampl = 1.2
#     ω = 0.1 / Ampl # Keep max velocity < 1.0

#     period = 2π / ω
#     num_periods = ceil(n_crossings / 4)
#     t_end = num_periods * period * 1.05

#     dt = 0.005
#     t_range = 0:dt:t_end
#     nt = length(t_range)

#     # --- Pre-calculate Physics ---
#     source_val = zeros(nt)
#     source_deriv = zeros(nt)
#     z_hist = zeros(nt)
#     v_hist = zeros(nt)
#     a_hist = zeros(nt)

#     # Cell Geometry
#     half_h = h / 2
#     _h = 1.0 / h
#     left_face = x_cell - half_h
#     right_face = x_cell + half_h

#     for (i, t) in enumerate(t_range)
#         # Kinematics
#         z = Ampl * sin(ω * t)
#         v = Ampl * ω * cos(ω * t)
#         acc = -Ampl * ω^2 * sin(ω * t)

#         z_hist[i] = z
#         v_hist[i] = v
#         a_hist[i] = acc

#         # Source Logic (Your function)
#         base_val = 0.5 * q * acc
#         val_L = -base_val / (1.0 + v)^2
#         val_R = base_val / (1.0 - v)^2

#         term_to_subtract = 0.0
#         if z >= right_face
#             term_to_subtract = val_L
#         elseif z <= left_face
#             term_to_subtract = val_R
#         else
#             frac_L = (z - left_face) * _h
#             frac_R = 1.0 - frac_L
#             term_to_subtract = frac_L * val_L + frac_R * val_R
#         end
#         source_val[i] = term_to_subtract
#     end

#     # Derivative
#     source_deriv[1:(end - 1)] = diff(source_val) ./ dt
#     source_deriv[end] = source_deriv[end - 1]

#     # --- Visualization Setup ---
#     fig = Figure(; size=(900, 900), fontsize=18, figure_padding=30)

#     # Observable Index for Animation
#     idx = Observable(1)

#     # Derived Observables
#     t_cur = lift(i -> t_range[i], idx)
#     z_cur = lift(i -> z_hist[i], idx)
#     v_cur = lift(i -> v_hist[i], idx)

#     # Field Profile
#     x_space = range(-1.8, 1.8; length=800)
#     field_profile = lift((z, v) -> begin
#                              map(x -> abs(x - z) / (1.0 - v * sign(x - z))^2, x_space)
#                          end, z_cur, v_cur)

#     # PLOT 1: Spatial Domain
#     ax_field = Axis(fig[1, 1]; title="Physical Context: Analytic Field & Cell",
#                     xlabel="Position (x)", ylabel="Field Strength Ψ(x)",
#                     limits=(-1.8, 1.8, -0.5, 15.0))

#     vspan!(ax_field, -0.5, 0.5; color=(:gray, 0.1))
#     vlines!(ax_field, [-0.5, 0.5]; color=:black, linestyle=:dash)
#     lines!(ax_field, x_space, field_profile; color=:blue, linewidth=3)

#     # Secondary Axis for Particle
#     ax_motion = Axis(fig[1, 1]; yaxisposition=:right, limits=(-1.8, 1.8, -1, 1), yticks=[0])
#     hidespines!(ax_motion)
#     hidexdecorations!(ax_motion)
#     scatter!(ax_motion, lift(z -> Point2f(z, 0), z_cur); color=:red, markersize=18)

#     # PLOT 2 & 3: Time Series
#     g_bottom = fig[2, 1] = GridLayout()
#     colgap!(g_bottom, 30)
#     ax_source = Axis(g_bottom[1, 1]; xlabel="t", ylabel="Source Term Value")
#     ax_deriv = Axis(g_bottom[1, 2]; xlabel="t", ylabel="Time Derivative")

#     lines!(ax_source, t_range, source_val; color=:black)
#     lines!(ax_deriv, t_range, source_deriv; color=:red)

#     vlines!(ax_source, t_cur; color=:gray)
#     vlines!(ax_deriv, t_cur; color=:gray)

#     pt_source = lift(i -> Point2f(t_range[i], source_val[i]), idx)
#     pt_deriv = lift(i -> Point2f(t_range[i], source_deriv[i]), idx)
#     scatter!(ax_source, pt_source; color=:gold, markersize=15, strokewidth=1)
#     scatter!(ax_deriv, pt_deriv; color=:red, markersize=15, strokewidth=1)

#     rowsize!(fig.layout, 1, Relative(0.6))

#     # --- Recording Loop ---
#     # Record every 4th frame to keep file size reasonable while smooth
#     step_size = 4
#     frames = 1:step_size:nt

#     record(fig, filename, frames; framerate=10) do i
#         idx[] = i
#     end
#     println("Saved animation to $filename")
# end

function plot_source_kink_harmonic(; n_crossings=4, paper=false)
    CairoMakie.activate!()

    set_theme!(mytheme_aps())

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
    # fig = Figure(; size=(900, 450), fontsize=18)
    fig = Figure(; size=(600, 250), figure_padding=(15, 15, 10, 10))
    # Grid layout for the two panels
    ax_source = Axis(fig[1, 1]; xlabel="t", title=L"\textrm{Source Term}")
    ax_deriv = Axis(fig[1, 2]; xlabel="t", title=L"\textrm{Time Derivative}")

    # 1. Add Shaded Regions (Background)
    for k in 1:length(starts)
        t_s = t_range[starts[k]]
        t_e = t_range[stops[k]]
        # Use gray shading for "Inside Cell" regions
        vspan!(ax_source, t_s, t_e; color=(:gray, 0.2))
        vspan!(ax_deriv, t_s, t_e; color=(:gray, 0.2))
    end

    # 2. Plot Curves
    lines!(ax_source, t_range, source_val; color=:black, linewidth=1.5)
    lines!(ax_deriv, t_range, source_deriv; color=:red, linewidth=1.5)
    # Save
    try
        savename = simple_name("source_kink_harmonic")
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
