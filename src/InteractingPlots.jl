module InteractingPlots

using ..PlottingUtils
using CairoMakie
using LinearAlgebra
# using GLMakie
using DataInterpolations
using ..ScalarField
using ..Energies

@inline name(pre, sim, post=nothing) = name_interacting(pre, sim, post)

@inline function KineticEnergy(m, v)
    γ = one(v) / sqrt(one(v) - v * v)
    return m * (γ - one(γ))
end

@inline function PotentialEnergy(xp, xs, Φ)
    interpolator_Φ = QuadraticInterpolation(Φ, xs)
    return interpolator_Φ(xp)
end

@inline function TotalEnergy(m, v, xp, xs, Φ)
    return m + KineticEnergy(m, v) + PotentialEnergy(xp, xs, Φ)
end

## More complicated!
## Integrate Hamiltonian of
# @inline function FieldEnergy(Π, Ψ)
#     return 0.5integrate(Π .^ 2 + Ψ .^ 2)
# end

function plot_pifield(i, sim1)
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 1]  # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp1 = sim1.sol.u[i].x[2][2]
    xp2 = sim1.sol.u[i].x[2][5]

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered)
    vlines!(ax, [xp1]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    vlines!(ax, [xp2]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("interacting_particles/pifield", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield(i, sim1)
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 2]  # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp1 = sim1.sol.u[i].x[2][2]
    xp2 = sim1.sol.u[i].x[2][5]

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered)
    vlines!(ax, [xp1]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    vlines!(ax, [xp2]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("interacting_particles/psifield", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield_resolutions(i, sim1, sim2, sim4)
    CairoMakie.activate!()
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

    xp1 = sim1.sol.u[i].x[2][2]
    xp2 = sim1.sol.u[i].x[2][5]

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
    vlines!(ax, [xp1]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    vlines!(ax, [xp2]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("interacting_particles/pifield_res_interacting", sim1, "i=$(i)")
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

    y1_raw = sim1.sol.u[i].x[1][:, 2]  # Matches x1_filtered length
    y1_filtered = y1_raw[mask1]

    y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 2])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 2]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp1 = sim1.sol.u[i].x[2][2]
    xp2 = sim1.sol.u[i].x[2][5]

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
    vlines!(ax, [xp1]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    vlines!(ax, [xp2]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("interacting_particles/psifield_res_interacting", sim1, "i=$(i)")
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

    savename = name("interacting_particles/pifield_convord_interacting", sim1)
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
    for i in 1:nt
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

    savename = name("interacting_particles/psifield_convord_interacting", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_positions(sim1; invert=true)
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xs2 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][2]
        xs2[i] = sim1.sol.u[i].x[2][5]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{particle 1}")
        lines!(ax, xs2, sim1.sol.t; label=L"\textrm{particle 2} ")
        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle 1}")
        lines!(ax, sim1.sol.t, xs2; label=L"\textrm{particle 2}")
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
    savename = name("interacting_particles/particle_positions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_positions_resolutions(sim1, sim2, sim4; invert=true)
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
        xs2l[i] = sim1.sol.u[i].x[2][5]
        xs1m[i] = sim2.sol.u[2i - 1].x[2][2]
        xs2m[i] = sim2.sol.u[2i - 1].x[2][5]

        xs1h[i] = sim4.sol.u[4i - 3].x[2][2]
        xs2h[i] = sim4.sol.u[4i - 3].x[2][5]
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
                          [L"\textrm{Particle 1}", L"\textrm{Particle 2}"];
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

    savename = name("interacting_particles/particle_positions_resolutions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_position_convord_interacting(sim1, sim2, sim4)
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
              title=L"\textrm{Convergence Order of Positions}",
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

    savename = name("interacting_particles/position_convord_interacting", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_velocities(sim1)
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    vs1 = zeros(nt)
    vs2 = zeros(nt)

    for i in 1:nt
        vs1[i] = sim1.sol.u[i].x[2][3]
        vs2[i] = sim1.sol.u[i].x[2][6]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle velocity}",
              limits=(0, sim1.sol.t[end], -1, 1))
    # Plot the data
    lines!(ax, sim1.sol.t, vs1; label=L"\mathrm{particle} 1")
    lines!(ax, sim1.sol.t, vs2; label=L"\mathrm{particle} 2")
    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    # Save the figure
    savename = name("interacting_particles/particle_velocities", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_velocities_resolutions(sim1, sim2, sim4)
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
        xs1l[i] = sim1.sol.u[i].x[2][3]
        xs2l[i] = sim1.sol.u[i].x[2][6]
        xs1m[i] = sim2.sol.u[2i - 1].x[2][3]
        xs2m[i] = sim2.sol.u[2i - 1].x[2][6]

        xs1h[i] = sim4.sol.u[4i - 3].x[2][3]
        xs2h[i] = sim4.sol.u[4i - 3].x[2][6]
    end

    fig = Figure(; size=(253, 200))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()
    leg2 = fig[3, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Velocities}", xlabelpadding=0)
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    lines!(ax, sim1.sol.t, xs1l; color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, xs2l; color=color2, linestyle=res_styles[1])

    lines!(ax, sim1.sol.t, xs1m; color=color1, linestyle=res_styles[2])
    lines!(ax, sim1.sol.t, xs2m; color=color2, linestyle=res_styles[2])

    lines!(ax, sim1.sol.t, xs1h; color=color1, linestyle=res_styles[3])
    lines!(ax, sim1.sol.t, xs2h; color=color2, linestyle=res_styles[3])

    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"

    color_legend = Legend(leg1[1, 1],
                          [LineElement(; color=color1, linestyle=res_styles[1]),
                           LineElement(; color=color2, linestyle=res_styles[1])],
                          [L"\textrm{Particle 1}", L"\textrm{Particle 2}"];
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
    rowgap!(fig.layout, 1, -5)
    rowgap!(fig.layout, 2, -10)

    savename = name("interacting_particles/particle_velocities_resolutions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_masses(sim1)
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    vs1 = zeros(nt)
    vs2 = zeros(nt)

    for i in 1:nt
        vs1[i] = sim1.sol.u[i].x[2][1]
        vs2[i] = sim1.sol.u[i].x[2][4]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}",
              limits=(0, sim1.sol.t[end], -1, 1))
    # Plot the data
    lines!(ax, sim1.sol.t, vs1; label=L"\mathrm{particle} 1")
    lines!(ax, sim1.sol.t, vs2; label=L"\mathrm{particle} 2")
    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    # Save the figure
    savename = name("interacting_particles/particle_masses", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_masses_resolutions(sim1, sim2, sim4)
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
        xs1l[i] = sim1.sol.u[i].x[2][1]
        xs2l[i] = sim1.sol.u[i].x[2][4]
        xs1m[i] = sim2.sol.u[2i - 1].x[2][1]
        xs2m[i] = sim2.sol.u[2i - 1].x[2][4]

        xs1h[i] = sim4.sol.u[4i - 3].x[2][1]
        xs2h[i] = sim4.sol.u[4i - 3].x[2][4]
    end

    fig = Figure(; size=(253, 200))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()
    leg2 = fig[3, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Masses}", xlabelpadding=0)
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    lines!(ax, sim1.sol.t, xs1l; color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, xs2l; color=color2, linestyle=res_styles[1])

    lines!(ax, sim1.sol.t, xs1m; color=color1, linestyle=res_styles[2])
    lines!(ax, sim1.sol.t, xs2m; color=color2, linestyle=res_styles[2])

    lines!(ax, sim1.sol.t, xs1h; color=color1, linestyle=res_styles[3])
    lines!(ax, sim1.sol.t, xs2h; color=color2, linestyle=res_styles[3])

    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"
    ylims!(ax, 0, nothing)
    color_legend = Legend(leg1[1, 1],
                          [LineElement(; color=color1, linestyle=res_styles[1]),
                           LineElement(; color=color2, linestyle=res_styles[1])],
                          [L"\textrm{Particle 1}", L"\textrm{Particle 2}"];
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
    rowgap!(fig.layout, 1, -5)
    rowgap!(fig.layout, 2, -10)

    savename = name("interacting_particles/particle_masses_resolutions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_panel(sim1, ticklabelsize=20)
    GLMakie.activate!()
    t = sim1.sol.t
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    xs1l = zeros(nt)
    vs1l = zeros(nt)
    ms1l = zeros(nt)

    xs2l = zeros(nt)
    vs2l = zeros(nt)
    ms2l = zeros(nt)

    for i in 1:nt
        ms1l[i] = sim1.sol.u[i].x[2][1]
        xs1l[i] = sim1.sol.u[i].x[2][2]
        vs1l[i] = sim1.sol.u[i].x[2][3]

        ms2l[i] = sim1.sol.u[i].x[2][4]
        xs2l[i] = sim1.sol.u[i].x[2][5]
        vs2l[i] = sim1.sol.u[i].x[2][6]
    end
    with_theme(mytheme) do
        fig = Figure(; size=(6 * 253, 4 * 200))
        supertitle_layout = fig[1, 1:2] = GridLayout()

        ax11 = Axis(fig[2, 1]; title=L"\textrm{Field}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)
        ax12 = Axis(fig[2, 2]; title=L"\textrm{Positions}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)
        ax21 = Axis(fig[3, 1]; title=L"\textrm{Masses}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)
        ax22 = Axis(fig[3, 2]; title=L"\textrm{Velocities}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)

        axes = [ax11, ax12, ax21, ax22]

        # Add slider below the plot grid
        slider_layout = fig[4, 1:2] = GridLayout()
        time_slider = Slider(slider_layout[1, 1];
                             range=1:nt,
                             startvalue=1,
                             horizontal=true)
        # Observable for current time
        current_index = time_slider.value

        palette = mytheme.palette
        color1 = palette.color[][1]  # Blue-ish for Particle 1
        color2 = palette.color[][2]  # Orange-ish for Particle 2
        res_styles = [palette.linestyle[][1], palette.linestyle[][2],
                      palette.linestyle[][3]]  # Solid, Dash, Dot

        lines!(ax12, xs1l, t; color=color1, linestyle=res_styles[1])
        lines!(ax12, xs2l, t; color=color2, linestyle=res_styles[1])
        ax12.xlabel = L"x(t)"
        ax12.ylabel = L"t"

        lines!(ax21, t, ms1l; color=color1, linestyle=res_styles[1])
        lines!(ax21, t, ms2l; color=color2, linestyle=res_styles[1])
        ax21.xlabel = L"t"
        ax21.ylabel = L"m(t)"

        lines!(ax22, t, vs1l; color=color1, linestyle=res_styles[1])
        lines!(ax22, t, vs2l; color=color2, linestyle=res_styles[1])
        ax22.xlabel = L"t"
        ax22.ylabel = L"v(t)"

        point_t = lift(current_index) do i
            return [t[i], t[i]]
        end
        point_x = lift(current_index) do i
            return [xs1l[i], xs2l[i]]
        end

        point_m = lift(current_index) do i
            return [ms1l[i], ms2l[i]]
        end
        point_v = lift(current_index) do i
            return [vs1l[i], vs2l[i]]
        end

        scatter!(ax12, point_x, point_t; color=:black, markersize=15)
        scatter!(ax21, point_t, point_m; color=:black, markersize=15)
        scatter!(ax22, point_t, point_v; color=:black, markersize=15)

        q1 = sim1.params.q1
        q2 = sim1.params.q2

        L = sim1.params.L
        x_range = (-L / 8, L / 8)
        mask = (sim1.params.x .>= -L / 8) .& (sim1.params.x .<= L / 8)
        x_filtered = sim1.params.x[mask]

        dynamic_plot = lift(current_index) do i
            empty!(ax11)  # Clear previous plot
            ti = t[i]
            x1 = sim1.sol.u[i].x[2][2]
            v1 = sim1.sol.u[i].x[2][3]
            x2 = sim1.sol.u[i].x[2][5]
            v2 = sim1.sol.u[i].x[2][6]

            Φs1 = ScalarField.Φs.(x_filtered, q1, x1, v1)
            Φs2 = ScalarField.Φs.(x_filtered, q2, x2, v2)
            y_raw = sim1.sol.u[i].x[1][:, 3]
            Φ_tot = y_raw[mask]
            y1 = Φ_tot .+ Φs2
            y2 = Φ_tot .+ Φs1
            interpolator_y1 = QuadraticInterpolation(y1, x_filtered;
                                                     extrapolation=ExtrapolationType.Extension)
            interpolator_y2 = QuadraticInterpolation(y2, x_filtered;
                                                     extrapolation=ExtrapolationType.Extension)

            lines!(ax11, x_filtered, y1; color=color2)
            lines!(ax11, x_filtered, y2; color=color1)
            scatter!(ax11, [x1], [interpolator_y1(x1)]; color=color1, marker=:circle,
                     markersize=15)
            scatter!(ax11, [x2], [interpolator_y2(x2)]; color=color2, marker=:circle,
                     markersize=15)
            autolimits!(ax11)
        end
        supertitle_text = lift(current_index) do i
            L"t=%$(t[i])"
        end

        index_label = lift(current_index) do i
            "Index: $i"
        end
        Label(supertitle_layout[1, 1], supertitle_text;
              fontsize=30,
              halign=:center,
              tellwidth=false,
              tellheight=false)

        Label(slider_layout[1, 2], index_label; tellwidth=false)

        rowsize!(fig.layout, 1, Auto(0.1))  # Small row for supertitle
        rowsize!(fig.layout, 2, Auto(1))    # Equal height for plot rows
        rowsize!(fig.layout, 3, Auto(1))
        rowsize!(fig.layout, 4, Auto(0.2))  # Small row for slider
        colgap!(slider_layout, 10)
        display(fig)
    end
    return nothing
end

function plot_particle_energies(sim1, sim2, sim4, ticklabelsize=20)
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    ke11 = zeros(nt)
    ke12 = zeros(nt)
    ke14 = zeros(nt)

    ke21 = zeros(nt)
    ke22 = zeros(nt)
    ke24 = zeros(nt)

    pe11 = zeros(nt)
    pe12 = zeros(nt)
    pe14 = zeros(nt)

    pe21 = zeros(nt)
    pe22 = zeros(nt)
    pe24 = zeros(nt)

    te11 = zeros(nt)
    te12 = zeros(nt)
    te14 = zeros(nt)
    te21 = zeros(nt)
    te22 = zeros(nt)
    te24 = zeros(nt)
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1],
                  palette.linestyle[][2],
                  palette.linestyle[][3]]  # Solid, Dash, Dot
    x1 = sim1.params.x
    x2 = sim2.params.x
    x4 = sim4.params.x
    q1 = sim1.params.q1
    q2 = sim1.params.q2

    for i in 1:nt
        xp11 = sim1.sol.u[i].x[2][2]
        xp12 = sim2.sol.u[2i - 1].x[2][2]
        xp14 = sim4.sol.u[4i - 3].x[2][2]

        xp21 = sim1.sol.u[i].x[2][5]
        xp22 = sim2.sol.u[2i - 1].x[2][5]
        xp24 = sim4.sol.u[4i - 3].x[2][5]

        mp11 = sim1.sol.u[i].x[2][1]
        mp12 = sim2.sol.u[2i - 1].x[2][1]
        mp14 = sim4.sol.u[4i - 3].x[2][1]

        mp21 = sim1.sol.u[i].x[2][4]
        mp22 = sim2.sol.u[2i - 1].x[2][4]
        mp24 = sim4.sol.u[4i - 3].x[2][4]

        vp11 = sim1.sol.u[i].x[2][3]
        vp12 = sim2.sol.u[2i - 1].x[2][3]
        vp14 = sim4.sol.u[4i - 3].x[2][3]

        vp21 = sim1.sol.u[i].x[2][6]
        vp22 = sim2.sol.u[2i - 1].x[2][6]
        vp24 = sim4.sol.u[4i - 3].x[2][6]

        Φ1 = sim1.sol.u[i].x[1][:, 3]
        Φ2 = sim2.sol.u[2i - 1].x[1][:, 3]
        Φ4 = sim4.sol.u[4i - 3].x[1][:, 3]

        ke11[i] = KineticEnergy(mp11, vp11)
        ke12[i] = KineticEnergy(mp12, vp12)
        ke14[i] = KineticEnergy(mp14, vp14)

        ke21[i] = KineticEnergy(mp21, vp21)
        ke22[i] = KineticEnergy(mp22, vp22)
        ke24[i] = KineticEnergy(mp24, vp24)

        pe11[i] = PotentialEnergy(xp11, x1, Φ1) + ScalarField.Φs(xp11, q2, xp21, vp21)
        pe12[i] = PotentialEnergy(xp12, x2, Φ2) + ScalarField.Φs(xp12, q2, xp22, vp22)
        pe14[i] = PotentialEnergy(xp14, x4, Φ4) + ScalarField.Φs(xp14, q2, xp24, vp24)

        pe21[i] = PotentialEnergy(xp21, x1, Φ1) + ScalarField.Φs(xp21, q1, xp12, vp11)
        pe22[i] = PotentialEnergy(xp22, x2, Φ2) + ScalarField.Φs(xp22, q1, xp12, vp12)
        pe24[i] = PotentialEnergy(xp24, x4, Φ4) + ScalarField.Φs(xp24, q1, xp14, vp14)

        te11[i] = mp11 + ke11[i] + pe11[i]
        te12[i] = mp12 + ke12[i] + pe12[i]
        te14[i] = mp14 + ke14[i] + pe14[i]

        te21[i] = mp21 + ke21[i] + pe21[i]
        te22[i] = mp22 + ke22[i] + pe22[i]
        te24[i] = mp24 + ke24[i] + pe24[i]
    end

    energy1 = Energies.FieldEnergy(sim1)
    energy2 = Energies.FieldEnergy(sim2)[1:2:end]
    energy4 = Energies.FieldEnergy(sim4)[1:4:end]

    fig = Figure(; size=(6 * 253, 4 * 200))
    ax11 = Axis(fig[1, 1]; title=L"\textrm{Kinetic Energy}",
                xticklabelsize=ticklabelsize,
                yticklabelsize=ticklabelsize)
    ax12 = Axis(fig[1, 2]; title=L"\textrm{Potential Energy}",
                xticklabelsize=ticklabelsize,
                yticklabelsize=ticklabelsize)
    ax21 = Axis(fig[2, 1]; title=L"\textrm{Total Energy}",
                xticklabelsize=ticklabelsize,
                yticklabelsize=ticklabelsize)
    ax22 = Axis(fig[2, 2]; title=L"\textrm{Field Energy}",
                xticklabelsize=ticklabelsize,
                yticklabelsize=ticklabelsize)

    ax11.xlabel = L"t"
    ax12.xlabel = L"t"
    ax21.xlabel = L"t"
    ax22.xlabel = L"t"

    lines!(ax11, sim1.sol.t, ke11; color=color1, linestyle=res_styles[1])
    lines!(ax11, sim1.sol.t, ke12; color=color1, linestyle=res_styles[2])
    lines!(ax11, sim1.sol.t, ke14; color=color1, linestyle=res_styles[3])

    lines!(ax11, sim1.sol.t, ke21; color=color2, linestyle=res_styles[1])
    lines!(ax11, sim1.sol.t, ke22; color=color2, linestyle=res_styles[2])
    lines!(ax11, sim1.sol.t, ke24; color=color2, linestyle=res_styles[3])

    lines!(ax12, sim1.sol.t, pe11; color=color1, linestyle=res_styles[1])
    lines!(ax12, sim1.sol.t, pe12; color=color1, linestyle=res_styles[2])
    lines!(ax12, sim1.sol.t, pe14; color=color1, linestyle=res_styles[3])

    lines!(ax12, sim1.sol.t, pe21; color=color2, linestyle=res_styles[1])
    lines!(ax12, sim1.sol.t, pe22; color=color2, linestyle=res_styles[2])
    lines!(ax12, sim1.sol.t, pe24; color=color2, linestyle=res_styles[3])

    lines!(ax21, sim1.sol.t, te11; color=color1, linestyle=res_styles[1])
    lines!(ax21, sim1.sol.t, te12; color=color1, linestyle=res_styles[2])
    lines!(ax21, sim1.sol.t, te14; color=color1, linestyle=res_styles[3])

    lines!(ax21, sim1.sol.t, te21; color=color2, linestyle=res_styles[1])
    lines!(ax21, sim1.sol.t, te22; color=color2, linestyle=res_styles[2])
    lines!(ax21, sim1.sol.t, te24; color=color2, linestyle=res_styles[3])

    lines!(ax22, sim1.sol.t, energy1; color=color1, linestyle=res_styles[1])
    lines!(ax22, sim1.sol.t, energy2; color=color1, linestyle=res_styles[2])
    lines!(ax22, sim1.sol.t, energy4; color=color1, linestyle=res_styles[3])

    display(fig)
    # savename = name("interacting_particles/particle_energies_resolutions", sim1)
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_panel(sim1, sim2, sim4, ticklabelsize=20)
    GLMakie.activate!()
    t = sim1.sol.t
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    xs1l = zeros(nt)
    vs1l = zeros(nt)
    ms1l = zeros(nt)

    xs1m = zeros(nt)
    vs1m = zeros(nt)
    ms1m = zeros(nt)

    xs1h = zeros(nt)
    vs1h = zeros(nt)
    ms1h = zeros(nt)

    xs2l = zeros(nt)
    vs2l = zeros(nt)
    ms2l = zeros(nt)

    xs2m = zeros(nt)
    vs2m = zeros(nt)
    ms2m = zeros(nt)

    xs2h = zeros(nt)
    vs2h = zeros(nt)
    ms2h = zeros(nt)

    for i in 1:nt
        ms1l[i] = sim1.sol.u[i].x[2][1]
        xs1l[i] = sim1.sol.u[i].x[2][2]
        vs1l[i] = sim1.sol.u[i].x[2][3]

        ms2l[i] = sim1.sol.u[i].x[2][4]
        xs2l[i] = sim1.sol.u[i].x[2][5]
        vs2l[i] = sim1.sol.u[i].x[2][6]

        ms1m[i] = sim2.sol.u[2i - 1].x[2][1]
        xs1m[i] = sim2.sol.u[2i - 1].x[2][2]
        vs1m[i] = sim2.sol.u[2i - 1].x[2][3]

        ms2m[i] = sim2.sol.u[2i - 1].x[2][4]
        xs2m[i] = sim2.sol.u[2i - 1].x[2][5]
        vs2m[i] = sim2.sol.u[2i - 1].x[2][6]

        ms1h[i] = sim4.sol.u[4i - 3].x[2][1]
        xs1h[i] = sim4.sol.u[4i - 3].x[2][2]
        vs1h[i] = sim4.sol.u[4i - 3].x[2][3]

        ms2h[i] = sim4.sol.u[4i - 3].x[2][4]
        xs2h[i] = sim4.sol.u[4i - 3].x[2][5]
        vs2h[i] = sim4.sol.u[4i - 3].x[2][6]
    end
    with_theme(mytheme) do
        fig = Figure(; size=(6 * 253, 4 * 200))
        supertitle_layout = fig[1, 1:2] = GridLayout()

        ax11 = Axis(fig[2, 1]; title=L"\textrm{\Phi Field}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)
        ax12 = Axis(fig[2, 2]; title=L"\textrm{Positions}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)
        ax21 = Axis(fig[3, 1]; title=L"\textrm{Masses}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)
        ax22 = Axis(fig[3, 2]; title=L"\textrm{Velocities}",
                    xticklabelsize=ticklabelsize,
                    yticklabelsize=ticklabelsize)

        axes = [ax11, ax12, ax21, ax22]

        # Add slider below the plot grid
        slider_layout = fig[4, 1:2] = GridLayout()
        time_slider = Slider(slider_layout[1, 1];
                             range=1:nt,
                             startvalue=1,
                             horizontal=true)
        # Observable for current time
        current_index = time_slider.value

        palette = mytheme.palette
        color1 = palette.color[][1]  # Blue-ish for Particle 1
        color2 = palette.color[][2]  # Orange-ish for Particle 2
        res_styles = [palette.linestyle[][1], palette.linestyle[][2],
                      palette.linestyle[][3]]  # Solid, Dash, Dot

        lines!(ax12, xs1l, t; color=color1, linestyle=res_styles[1])
        lines!(ax12, xs1m, t; color=color1, linestyle=res_styles[2])
        lines!(ax12, xs1h, t; color=color1, linestyle=res_styles[3])

        lines!(ax12, xs2l, t; color=color2, linestyle=res_styles[1])
        lines!(ax12, xs2m, t; color=color2, linestyle=res_styles[2])
        lines!(ax12, xs2h, t; color=color2, linestyle=res_styles[3])
        ax12.xlabel = L"x(t)"
        ax12.ylabel = L"t"

        lines!(ax21, t, ms1l; color=color1, linestyle=res_styles[1])
        lines!(ax21, t, ms1m; color=color1, linestyle=res_styles[2])
        lines!(ax21, t, ms1h; color=color1, linestyle=res_styles[3])

        lines!(ax21, t, ms2l; color=color2, linestyle=res_styles[1])
        lines!(ax21, t, ms2m; color=color2, linestyle=res_styles[2])
        lines!(ax21, t, ms2h; color=color2, linestyle=res_styles[3])

        ax21.xlabel = L"t"
        ax21.ylabel = L"m(t)"

        lines!(ax22, t, vs1l; color=color1, linestyle=res_styles[1])
        lines!(ax22, t, vs1m; color=color1, linestyle=res_styles[2])
        lines!(ax22, t, vs1h; color=color1, linestyle=res_styles[3])

        lines!(ax22, t, vs2l; color=color2, linestyle=res_styles[1])
        lines!(ax22, t, vs2m; color=color2, linestyle=res_styles[2])
        lines!(ax22, t, vs2h; color=color2, linestyle=res_styles[3])
        ax22.xlabel = L"t"
        ax22.ylabel = L"v(t)"

        point_t = lift(current_index) do i
            return [t[i], t[i]]
        end
        point_x = lift(current_index) do i
            return [xs1l[i], xs2l[i]]
        end

        point_m = lift(current_index) do i
            return [ms1l[i], ms2l[i]]
        end
        point_v = lift(current_index) do i
            return [vs1l[i], vs2l[i]]
        end

        scatter!(ax12, point_x, point_t; color=:black, markersize=15)
        scatter!(ax21, point_t, point_m; color=:black, markersize=15)
        scatter!(ax22, point_t, point_v; color=:black, markersize=15)

        q1 = sim1.params.q1
        q2 = sim1.params.q2

        L = sim1.params.L
        x_range = (-L / 8, L / 8)
        mask = (sim1.params.x .>= -L / 8) .& (sim1.params.x .<= L / 8)
        x_filtered = sim1.params.x[mask]

        dynamic_plot = lift(current_index) do i
            empty!(ax11)  # Clear previous plot
            ti = t[i]
            x1 = sim1.sol.u[i].x[2][2]
            v1 = sim1.sol.u[i].x[2][3]
            x2 = sim1.sol.u[i].x[2][5]
            v2 = sim1.sol.u[i].x[2][6]

            Φs1 = ScalarField.Φs.(x_filtered, q1, x1, v1)
            Φs2 = ScalarField.Φs.(x_filtered, q2, x2, v2)
            y_raw = sim1.sol.u[i].x[1][:, 3]
            Φ_tot = y_raw[mask]
            y1 = Φ_tot .+ Φs2
            y2 = Φ_tot .+ Φs1
            interpolator_y1 = QuadraticInterpolation(y1, x_filtered;
                                                     extrapolation=ExtrapolationType.Extension)
            interpolator_y2 = QuadraticInterpolation(y2, x_filtered;
                                                     extrapolation=ExtrapolationType.Extension)

            lines!(ax11, x_filtered, y1; color=color2)
            lines!(ax11, x_filtered, y2; color=color1)
            scatter!(ax11, [x1], [interpolator_y1(x1)]; color=color1, marker=:circle,
                     markersize=15)
            scatter!(ax11, [x2], [interpolator_y2(x2)]; color=color2, marker=:circle,
                     markersize=15)
            autolimits!(ax11)
        end
        supertitle_text = lift(current_index) do i
            L"t=%$(t[i])"
        end

        index_label = lift(current_index) do i
            "Index: $i"
        end
        Label(supertitle_layout[1, 1], supertitle_text;
              fontsize=30,
              halign=:center,
              tellwidth=false,
              tellheight=false)

        Label(slider_layout[1, 2], index_label; tellwidth=false)

        rowsize!(fig.layout, 1, Auto(0.1))  # Small row for supertitle
        rowsize!(fig.layout, 2, Auto(1))    # Equal height for plot rows
        rowsize!(fig.layout, 3, Auto(1))
        rowsize!(fig.layout, 4, Auto(0.2))  # Small row for slider
        colgap!(slider_layout, 10)
        display(fig)
    end
    return nothing
end

end
