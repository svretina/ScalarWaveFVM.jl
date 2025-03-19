module ParticlePotentialPlots

using ..PlottingUtils
using CairoMakie
# using GLMakie
using LinearAlgebra

@inline name(pre, sim, post=nothing) = name_potential(pre, sim, post)
@inline name_analytic(pre, sim, post=nothing) = name_analytic_potential(pre, sim, post)

function plot_pifield_resolutions(i, sim1, sim2, sim4)
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 1])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp = sim1.sol.u[i].x[2][2]

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("particle_potential/pifield", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_resolutions(i, sim1, sim2, sim4)
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @show t1, t2, t4
    @assert t1 == t2
    @assert t2 == t4
    t = t1

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 2]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 2])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 2]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp = sim1.sol.u[i].x[2][2]

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("particle_potential/psifield", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield_convord(sim1, sim2, sim4)
    set_theme!(mytheme_aps())
    nt = length(sim1.sol.t)
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_raw = sim1.sol.u[i].x[1][:, 1]  # Matches x1_filtered length
        y1_filtered = y1_raw[mask1]

        y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 1])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    fig = Figure(; ize=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Order of \Pi_\mathrm{r} field}",
              limits=(0, sim1.sol.t[end], 0, 4))
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("particle_potential/pifield_convord", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_convord(sim1, sim2, sim4)
    set_theme!(mytheme_aps())
    nt = length(sim1.sol.t)
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_raw = sim1.sol.u[i].x[1][:, 2]  # Matches x1_filtered length
        y1_filtered = y1_raw[mask1]

        y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 2])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 2]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    fig = Figure(; ize=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Order of \Psi_\mathrm{r} field}",
              limits=(0, sim1.sol.t[end], 0, 4))
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("particle_potential/psifield_convord", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_position(sim1; invert=true)
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][2]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{particle 1}")
        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle 1}")
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name("particle_potential/particle_position", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_positions_resolutions(sim1, sim2, sim4, sim1sf, sim2sf, sim4sf;
                                             invert=true)
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1l = zeros(nt)
    xs1m = zeros(nt)
    xs1h = zeros(nt)

    xs2l = zeros(nt)
    xs2m = zeros(nt)
    xs2h = zeros(nt)

    for i in 1:nt
        xs1l[i] = sim1.sol.u[i].x[2][2]
        xs2l[i] = sim1sf.sol.u[i].x[2][2]

        xs1m[i] = sim2.sol.u[2i - 1].x[2][2]
        xs2m[i] = sim2sf.sol.u[2i - 1].x[2][2]

        xs1h[i] = sim4.sol.u[4i - 3].x[2][2]
        xs2h[i] = sim4sf.sol.u[4i - 3].x[2][2]
    end

    fig = Figure(; size=(253, 200))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()
    leg2 = fig[3, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    if invert
        lines!(ax, xs1l, sim1.sol.t; color=color1, linestyle=res_styles[1])
        lines!(ax, xs2l, sim1.sol.t; color=color2, linestyle=res_styles[1])

        lines!(ax, xs1m, sim1.sol.t; color=color1, linestyle=res_styles[2])
        lines!(ax, xs2m, sim1.sol.t; color=color2, linestyle=res_styles[2])

        lines!(ax, xs1h, sim1.sol.t; color=color1, linestyle=res_styles[3])
        lines!(ax, xs2h, sim1.sol.t; color=color2, linestyle=res_styles[3])

        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1)
        lines!(ax, sim1.sol.t, xs2)
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    color_legend = Legend(leg1[1, 1],
                          [LineElement(; color=color1, linestyle=res_styles[1]),
                           LineElement(; color=color2, linestyle=res_styles[1])],
                          [L"\textrm{Particle no SF}", L"\textrm{Particle SF}"];
                          halign=:center,
                          valign=:top,
                          tellwidth=false,
                          tellheight=true,
                          orientation=:horizontal)

    # # 2. Resolution legend (Low, Med, High Res)
    res_legend = Legend(leg2[1, 1],
                        [LineElement(; color=:black, linestyle=res_styles[1]),
                         LineElement(; color=:black, linestyle=res_styles[2]),
                         LineElement(; color=:black, linestyle=res_styles[3])],
                        [L"\textrm{low res}", L"\textrm{med res}", L"\textrm{high res}"];
                        halign=:center,
                        valign=:top,
                        tellwidth=false,
                        tellheight=true, orientation=:horizontal)

    # Place legends below the axis
    rowgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 2, -10)

    savename = name("particle_potential/particle_positions_resolutions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_position_resolutions(sim1, sim2, sim4; invert=true)
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1l = zeros(nt)
    xs1m = zeros(nt)
    xs1h = zeros(nt)

    for i in 1:nt
        xs1l[i] = sim1.sol.u[i].x[2][2]
        xs1m[i] = sim2.sol.u[2i - 1].x[2][2]
        xs1h[i] = sim4.sol.u[4i - 3].x[2][2]
    end

    fig = Figure(; size=(253, 200))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    if invert
        lines!(ax, xs1l, sim1.sol.t; color=color1, linestyle=res_styles[1])
        lines!(ax, xs1m, sim1.sol.t; color=color1, linestyle=res_styles[2])
        lines!(ax, xs1h, sim1.sol.t; color=color1, linestyle=res_styles[3])

        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1)
        lines!(ax, sim1.sol.t, xs2)
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    # # 2. Resolution legend (Low, Med, High Res)
    res_legend = Legend(leg1[1, 1],
                        [LineElement(; color=:black, linestyle=res_styles[1]),
                         LineElement(; color=:black, linestyle=res_styles[2]),
                         LineElement(; color=:black, linestyle=res_styles[3])],
                        [L"\textrm{low res}", L"\textrm{med res}", L"\textrm{high res}"];
                        halign=:center,
                        valign=:top,
                        tellwidth=false,
                        tellheight=true, orientation=:horizontal)

    # Place legends below the axis
    rowgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 2, -10)

    savename = name("particle_potential/particle_position_resolutions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

### Particle rhs only
###
###

function plot_particle_position_analytic(sim1; invert=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    # Create figure and axis
    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    # vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i][2]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{particle 1}")
        lines!(ax, analytic.(sim1.sol.t), sim1.sol.t; label=L"\textrm{analytic}")
        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle 1}")
        lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
               label=L"\textrm{analytic}")
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_position_analytic_res(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    # Create figure and axis
    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xs2 = zeros(nt)
    xs4 = zeros(nt)

    # vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i][2]
        xs2[i] = sim2.sol.u[2i - 1][2]
        xs4[i] = sim4.sol.u[4i - 3][2]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{particle 1}",
               color=color1, linestyle=res_styles[1])
        lines!(ax, xs2, sim1.sol.t; color=color2, linestyle=res_styles[2])
        lines!(ax, xs4, sim1.sol.t; linestyle=res_styles[3])
        # lines!(ax, analytic.(sim1.sol.t), sim1.sol.t; label=L"\textrm{analytic}")
        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle 1}",
               color=color1, linestyle=res_styles[1])
        lines!(ax, sim1.sol.t, xs2; color=color1, linestyle=res_styles[2])
        lines!(ax, sim1.sol.t, xs4; color=color1, linestyle=res_styles[3])
        lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
               label=L"\textrm{analytic}")
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position_res", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_position_analytic_diff(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)

    err1 = zeros(nt)
    err2 = zeros(nt)
    for i in 1:nt
        err1[i] = sim1.sol.u[i][2] - sim2.sol.u[2i - 1][2]
        err2[i] = sim2.sol.u[2i - 1][2] - sim4.sol.u[4i - 3][2]
    end
    # convfac = @. (xs1 - xs2) / (xs2 - xs4)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    lines!(ax, sim1.sol.t, err1; label=L"\textrm{low-mid}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, 2^4 * err2; label=L"\textrm{2^4\cdot(mid-high)}",
           color=color2, linestyle=res_styles[2])

    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position_diff", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_position_analytic_convfac(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xs2 = zeros(nt)
    xs4 = zeros(nt)
    convfac = zeros(nt)
    # vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i][2]
        xs2[i] = sim2.sol.u[2i - 1][2]
        xs4[i] = sim4.sol.u[4i - 3][2]
        err1 = xs1[i] - xs2[i]
        err2 = xs2[i] - xs4[i]
        if err1 == err2
            convfac[i] = 1.0
        else
            convfac[i] = err1 / err2
        end
    end
    # convfac = @. (xs1 - xs2) / (xs2 - xs4)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position Convergence Factor}")
    # palette = mytheme.palette
    # color1 = palette.color[][1]  # Blue-ish for Particle 1
    # color2 = palette.color[][2]  # Orange-ish for Particle 2
    # res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    lines!(ax, sim1.sol.t, convfac)
    ylims!(ax, 10, 40)
    ax.xlabel = L"t"
    # leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    #rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    #rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position_convfac", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_velocity_analytic(sim1; invert=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        vs1[i] = sim1.sol.u[i][3]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity}")
    # Plot the data
    lines!(ax, sim1.sol.t, vs1; label=L"\textrm{particle 1}")
    lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
           label=L"\textrm{analytic}")
    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_velocity_analytic_res(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    vs1 = zeros(nt)
    vs2 = zeros(nt)
    vs4 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        vs1[i] = sim1.sol.u[i][3]
        vs2[i] = sim2.sol.u[2i - 1][3]
        vs4[i] = sim4.sol.u[4i - 3][3]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, vs1; label=L"\textrm{particle 1}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, vs2; color=color1, linestyle=res_styles[2])
    lines!(ax, sim1.sol.t, vs4; color=color1, linestyle=res_styles[3])
    lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
           label=L"\textrm{analytic}")
    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity_res", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_velocity_analytic_diff(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    err1 = zeros(nt)
    err2 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        err1[i] = sim1.sol.u[i][3] - sim2.sol.u[2i - 1][3]
        err2[i] = sim2.sol.u[2i - 1][3] - sim4.sol.u[4i - 3][3]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, err1; label=L"\textrm{low-mid}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, 2^4 * err2; label=L"\textrm{2^4 \cdot (mid-high)}",
           color=color2, linestyle=res_styles[2])
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity_diff", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_velocity_analytic_convfac(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    convfac = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        err1 = sim1.sol.u[i][3] - sim2.sol.u[2i - 1][3]
        err2 = sim2.sol.u[2i - 1][3] - sim4.sol.u[4i - 3][3]
        if err1 == err2
            convfac[i] = 1.0
        else
            convfac[i] = err1 / err2
        end
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity Convergence Factor}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot
    ylims!(ax, 0, 40)
    # Plot the data
    lines!(ax, sim1.sol.t, convfac; color=color1, linestyle=res_styles[1])
    ax.xlabel = L"t"
    #leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    # rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity_convfac", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_mass_analytic(sim1; invert=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # Create figure and axis
    nt = length(sim1.sol.t)
    ms1 = zeros(nt)
    for i in 1:nt
        ms1[i] = sim1.sol.u[i][1]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    lines!(ax, sim1.sol.t, ms1; label=L"\textrm{particle 1}")
    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_mass", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_mass_analytic_res(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    # Create figure and axis
    nt = length(sim1.sol.t)
    ms1 = zeros(nt)
    ms2 = zeros(nt)
    ms4 = zeros(nt)
    for i in 1:nt
        ms1[i] = sim1.sol.u[i][1]
        ms2[i] = sim2.sol.u[2i - 1][1]
        ms4[i] = sim4.sol.u[4i - 3][1]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, ms1; label=L"\textrm{particle 1}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, ms2; color=color1, linestyle=res_styles[2])
    lines!(ax, sim1.sol.t, ms4; color=color1, linestyle=res_styles[3])

    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_mass_res", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_mass_analytic_diff(sim1, sim2, sim4; invert=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    err1 = zeros(nt)
    err2 = zeros(nt)
    for i in 1:nt
        err1[i] = sim1.sol.u[i][1] - sim2.sol.u[2i - 1][1]
        err2[i] = sim2.sol.u[2i - 1][1] - sim4.sol.u[4i - 3][1]
    end

    fig = Figure()
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, err1; label=L"\textrm{low-mid}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, 2^4 * err2; label=L"\textrm{2^4 \cdot (mid-high)}",
           color=color2, linestyle=res_styles[2])

    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_mass_diff", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

end
