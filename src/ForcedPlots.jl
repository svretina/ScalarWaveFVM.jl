module ForcedPlots

using ..PlottingUtils
using ..ParticleMotion
using ..ScalarField
using CairoMakie
using LinearAlgebra
using DataInterpolations

@inline name(pre, sim, post=nothing) = name_forced(pre, sim, post)

function trapezoidal(y, h)
    n = length(y) - 1
    integral = zero(eltype(y))
    @inbounds for i in 1:n
        integral += y[i]
    end
    @inbounds integral += 0.5 * (y[1] + y[n + 1])  # Adjust endpoints
    integral *= h
    return integral
end

function deriv(y, h)
    n = length(y)
    dy = zeros(eltype(y), n)
    for i in 2:(n - 1)
        dy[i] = (y[i + 1] - y[i - 1]) / (2h)
    end

    dy[1] = (y[1] - y[2]) / h
    dy[n] = (y[n - 1] - y[n]) / h
    return dy
end

function plot_particle_position(sim, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim.sol.t)
    # Create figure and axis

    L = sim.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim.params.x .>= -L / 2) .& (sim.params.x .<= L / 2)
    x_filtered = sim.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    dt = sim.params.dx * sim.params.cfl
    pos = zeros(nt)
    for i in 1:nt
        pos[i], _, _ = ParticleMotion.oscillator(sim.sol.t[i], sim.params.z0, sim.params.A,
                                                 sim.params.ω)
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"t",
              title="Particle Position",)
    # Plot the data

    lines!(ax, pos, sim.sol.t)

    # vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/particle_position", sim)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_pifield(i, sim1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t = t1

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/pifield_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_psifield(i, sim1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t = t1

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/psifield_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_phifield(i, sim1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t = t1

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 3]  # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/phifield_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_pifield_resolutions(i, sim1, sim2, sim4, F=8, paper=false)
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
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/pifield_res_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_total_pifield(i, sim1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis

    L = sim1.params.L
    F = 2
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Pis = ScalarField.Πs.(x_filtered, xp, vp, sim1.params.q)
    xp0, vp0, _ = ParticleMotion.oscillator(0.0, sim1.params.z0, sim1.params.A,
                                            sim1.params.ω)
    Pis0 = ScalarField.Πs.(x_filtered, xp0, vp0, sim1.params.q)
    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = Pis .+ sim1.sol.u[i][mask1, 1] .- Pis0 # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = round(0.1 * (y_max - y_min); digits=4)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=3)

    @show padding
    @show y_limits
    @show y_max
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered)
    # lines!(ax, x_filtered, sim1.sol.u[i][mask1, 1])
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/pifield_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_total_psifield(i, sim1, F=2, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Psis = ScalarField.Ψs.(x_filtered, xp, vp, sim1.params.q)
    xp0, vp0, _ = ParticleMotion.oscillator(0.0, sim1.params.z0, sim1.params.A,
                                            sim1.params.ω)
    Psis0 = ScalarField.Ψs.(x_filtered, xp0, vp0, sim1.params.q)

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = Psis .+ sim1.sol.u[i][mask1, 2] .- Psis0 # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = round(0.1 * (y_max - y_min); digits=4)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=4)

    @show padding
    @show y_limits
    @show y_min, y_max
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    # lines!(ax, x_filtered, y1_filtered)

    lines!(ax, x_filtered, y1_filtered)

    # lines!(ax, x_filtered, Psis0)
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/psifield_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_total_phifield(i, sim1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis
    F = 2
    L = sim1.params.L
    @show L
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Phis = ScalarField.Φs.(x_filtered, Ref(xp), Ref(vp), Ref(sim1.params.q))

    xp0, vp0, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Phis0 = ScalarField.Φs.(x_filtered, Ref(xp0), Ref(vp0), Ref(sim1.params.q))

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]
    y1_filtered = Phis .+ sim1.sol.u[i][mask1, 3] .- Phis0

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = round(0.1 * (y_max - y_min); digits=4)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=4)

    @show padding
    @show y_limits
    @show y_max

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    # Plot the data

    lines!(ax, x_filtered, y1_filtered)

    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/phifield_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_total_pifield_resolutions(i, sim1, sim2, sim4, F=8, paper=false)
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
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Pis = ScalarField.Πs.(x_filtered, xp, vp, sim1.params.q)
    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = Pis .+ sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
    y2_filtered = Pis .+ y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
    y4_filtered = Pis .+ y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = round(0.1 * (y_max - y_min); digits=3)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=3)

    @show padding
    @show y_limits
    @show y_max
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/pifield_res_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_psifield_resolutions(i, sim1, sim2, sim4, F=8, paper=false)
    CairoMakie.activate!()
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
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/psifield_res_forced", sim1, "i=$(i)")
    # save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_pifield_diff(i, sim1, sim2, sim4, paper=false)
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
    F = 2
    L = sim1.params.L
    x_range = (-L / F, L / F)
    @show x_range

    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    diff1 = y1_filtered .- y2_filtered
    diff2 = y2_filtered .- y4_filtered
    # Compute y-limits with padding
    y_min = minimum([minimum(diff1), minimum(diff2)])
    y_max = maximum([maximum(diff1), maximum(diff2)])
    padding = 0.2 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    # Plot the data
    lines!(ax, x_filtered, diff1; label=L"\textrm{low-mid}")
    lines!(ax, x_filtered, 2^2 * diff2; label=L"\textrm{2^2\cdot(mid-high)}",
           linestyle=:dash)
    # lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 10% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/pifield_diff_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_psifield_diff(i, sim1, sim2, sim4, paper=false)
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
    F = 2
    L = sim1.params.L
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    diff1 = y1_filtered .- y2_filtered
    diff2 = y2_filtered .- y4_filtered
    # Compute y-limits with padding
    y_min = minimum([minimum(diff1), minimum(diff2)])
    y_max = maximum([maximum(diff1), maximum(diff2)])
    padding = 0.3 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, diff1; label=L"\textrm{low-mid}")
    lines!(ax, x_filtered, 2^2 * diff2; label=L"\textrm{2^2\cdot(mid-high)}",
           linestyle=:dot)
    # lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/psifield_diff_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_pifield_convfac(sim1, sim2, sim4, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    # Create figure and axis
    F = 2
    L = sim1.params.L

    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # Filter sim2 data (1280 points)
    #mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    # mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 1:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements
        #
        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1) / norm(diff2)
        @show convfac[i]
    end
    # Compute y-limits with padding
    y_min = minimum(convfac)
    y_max = maximum(convfac)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Factor of \Pi_\mathrm{r} field}",
              limits=(sim1.sol.t[2], sim1.sol.t[end], 0, 4),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], convfac[2:end]; linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("forced_motion/pifield_convfac_forced", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_psifield_convfac(sim1, sim2, sim4, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # Filter sim2 data (1280 points)
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 1:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]

        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements
        #
        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1) / norm(diff2)
    end
    # Compute y-limits with padding
    y_min = minimum(convfac)
    y_max = maximum(convfac)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; ize=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Factor of \Psi_\mathrm{r} field}",
              limits=(sim1.sol.t[2], sim1.sol.t[end], 1, 3),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], convfac[2:end]; linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("forced_motion/psifield_convfac_forced", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_pifield_convord(sim1, sim2, sim4, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # Filter sim2 data (1280 points)
    # mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    # mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]

        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    # Compute y-limits with padding
    # y_min = minimum(convfac)
    # y_max = maximum(convfac)
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; ize=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Self Convergence Order of \Pi_\mathrm{r} field}",
              limits=(0, sim1.sol.t[end], 0, 4),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("forced_motion/pifield_convord_forced", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_psifield_convord(sim1, sim2, sim4, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    vmax = sim1.params.vmax
    A = sim1.params.A

    nt = length(sim1.sol.t)
    # Create figure and axis
    A = 1.0
    ω = round(vmax / A; digits=9)
    T = (2π / ω)
    Nosc = sim1.params.Nosc
    tf = Nosc * T

    L = sim1.params.L
    # x_range = (-L / 2, L / 2)
    x_range = (-tf, tf)

    # mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    mask = (sim1.params.x .>= -tf) .& (sim1.params.x .<= tf)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # Filter sim2 data (1280 points)
    # mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    # mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    # Compute y-limits with padding
    # y_min = minimum(convfac)
    # y_max = maximum(convfac)
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; ize=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Order of \Psi_\mathrm{r} field}",
              limits=(0, sim1.sol.t[end], 0, 4),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("forced_motion/psifield_convord_forced", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_convergence_panel(sim1, sim2, sim4, paper=false)
    nt = length(sim1.sol.t)
    fig1 = plot_pifield_resolutions(nt, sim1, sim2, sim4, paper)  # Your existing function
    fig2 = plot_psifield_resolutions(nt, sim1, sim2, sim4, paper)  # Replace with your second function
    fig3 = plot_pifield_convord(sim1, sim2, sim4, paper)  # Replace with your third function
    fig4 = plot_psifield_convord(sim1, sim2, sim4, paper)  # Replace with your fourth function

    fig = Figure(; size=(600, 500))

    # Add supertitle
    supertitle = Label(fig[1, 1:2],
                       L"v_{max}=%$(sim1.params.vmax)";
                       fontsize=20,
                       font=:bold,
                       halign=:center,)

    # Function to recreate Axis content
    function recreate_axis!(fig, src_fig, row, col)
        src_ax = contents(src_fig[1, 1])[1]  # Original Axis

        # Create new Axis with same properties
        new_ax = Axis(fig[row, col];
                      xlabel=src_ax.xlabel[],
                      ylabel=src_ax.ylabel[],
                      title=src_ax.title[],
                      limits=src_ax.limits[],)

        # Copy plot objects
        for plot in src_ax.scene.plots
            if plot isa MakieCore.Scatter
                label = haskey(plot.attributes, :label) ? plot.label[] : nothing
                scatter!(new_ax,
                         plot[1][],
                         plot[2][];
                         label=label,
                         markersize=plot.markersize[],)
            elseif plot isa MakieCore.Lines
                data = plot[1][]  # Single Observable
                label = haskey(plot.attributes, :label) ? plot.label[] : nothing
                if isa(data, Vector{Point{2,Float64}})  # Regular Lines
                    x = [p[1] for p in data]
                    y = [p[2] for p in data]
                    lines!(new_ax, x, y; label=label, linewidth=plot.linewidth[])
                elseif isa(data, Vector{Float64})  # VLines misclassified as Lines
                    @show "ASDFA"
                    vlines!(new_ax,
                            data;
                            color=plot.color[],
                            linestyle=plot.linestyle[],
                            linewidth=plot.linewidth[],)
                end
            end
        end

        # Add Legend if it exists, place in the next row under the column
        leg_content = contents(src_fig[2, 1])
        if !isempty(leg_content)
            Legend(fig[row + 1, col],
                   new_ax;
                   orientation=:horizontal,
                   tellwidth=false,
                   tellheight=true,)
        end
    end

    # Recreate each plot in the grid
    recreate_axis!(fig, fig1, 2, 1)  # Top-left: Π field resolutions (Axis in row 2, Legend in row 3)
    recreate_axis!(fig, fig2, 2, 2)  # Top-right: Ψ field resolutions (Axis in row 2, Legend in row 3)
    recreate_axis!(fig, fig3, 4, 1)  # Bottom-left: Π field convergence (Axis in row 4)
    recreate_axis!(fig, fig4, 4, 2)  # Bottom-right: Ψ field convergence (Axis in row 4)

    # Adjust layout with padding on the right
    rowsize!(fig.layout, 1, Relative(0.1)) # Supertitle row
    rowsize!(fig.layout, 2, Auto())       # Top row of plots (fig1, fig2)
    rowsize!(fig.layout, 3, Relative(0.05)) # Legend row for fig1 and fig2
    rowsize!(fig.layout, 4, Auto())       # Bottom row of plots (fig3, fig4)
    # rowsize!(fig.layout, 5, Relative(0.0)) # Empty row
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.5)) # Left column (50% of width)
    colsize!(fig.layout, 2, Relative(0.45)) # Right column (45% of width, leaves 5% padding)
    # Save the combined figure

    savename = name("forced_motion/convergence_pane_forced", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_balance_law(sim1)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    # F = 6
    L = sim1.params.L
    Nosc = sim1.params.Nosc
    ω = sim1.params.ω
    A = sim1.params.A
    x0 = sim1.params.x0
    T = 2π / ω
    tf = Nosc * T
    L0 = 2tf
    x_range = (-L0 / 4, L0 / 4)
    @show x_range
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    particle_energy = zeros(nt)
    potential = zeros(nt)
    kinetic = zeros(nt)
    self_potential = zeros(nt)
    regular_field_energy = zeros(nt)
    # total_energyr = zeros(nt)
    # total_energys = zeros(nt)
    flux = zeros(nt)
    flux_analytic = zeros(nt)
    bl = zeros(nt)
    dxTi0 = zeros(nt, length(x_filtered))
    Hfr = zeros(nt, length(x_filtered))
    Hfr2 = zeros(nt, length(x_filtered))
    dHpdt_analytic = zeros(nt)
    dHfrdt_num = zeros(nt, length(x_filtered))
    dHfrdt_analytic22 = zeros(nt, length(x_filtered))
    dHfrdt_analytic12 = zeros(nt, length(x_filtered))
    tmp = zeros(nt, length(x_filtered))
    # fluxes2 = zeros(nt)
    # fluxes3 = zeros(nt)
    # fluxes4 = zeros(nt)

    dHfrdt_analytic = zeros(nt)

    q = sim1.params.q
    dx = sim1.params.dx
    xs = sim1.params.x
    cfl = sim1.params.cfl
    dt = cfl * dx

    function deriv(y, h)
        n = length(y)
        dy = zeros(eltype(y), n)
        for i in 2:(n - 1)
            dy[i] = (y[i + 1] - y[i - 1]) / (2h)
        end

        dy[1] = (y[1] - y[2]) / h
        dy[n] = (y[n - 1] - y[n]) / h
        return dy
    end
    function delta(x, z, h)
        x_left = x - h / 2
        x_right = x + h / 2
        return (z >= x_left && z < x_right) ? 1.0 : 0.0
    end

    for i in 1:nt
        ### calculate Acceleration
        ti = sim1.sol.t[i]
        z, v, a = ParticleMotion.oscillator(ti, x0, A, ω)
        ####
        δ = delta.(x_filtered, z, dx)
        ## This is with Φr and z as a generalized coordinate:
        Φr = @view sim1.sol.u[i][mask, 3]
        Πr = @view sim1.sol.u[i][mask, 1]
        Ψr = @view sim1.sol.u[i][mask, 2]

        interpolator_Φ = AkimaInterpolation(Φr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Π = AkimaInterpolation(Πr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Ψ = AkimaInterpolation(Ψr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)

        Φrz = interpolator_Φ(z)
        Πrz = interpolator_Π(z)
        Ψrz = interpolator_Ψ(z)
        # Φrz = Φr .* δ
        # Πrz = Πr .* δ
        # Ψrz = Ψr .* δ

        Φs = ScalarField.Φs.(x_filtered, z, v, q)
        Πs = ScalarField.Πs.(x_filtered, z, v, q)
        Ψs = ScalarField.Ψs.(x_filtered, z, v, q)

        ∂vΦs = ScalarField.Φs001.(x_filtered, z, v, q)
        ∂vΠs = ScalarField.Πs001.(x_filtered, z, v, q)
        ∂vΨs = ScalarField.Ψs001.(x_filtered, z, v, q)
        ∂zΠs = ScalarField.Πs010.(x_filtered, z, v, q)
        ∂xΨs = ScalarField.Ψs100.(x_filtered, z, v, q)

        dxΠr = deriv(Πr, dx)
        dxΨr = deriv(Ψr, dx)

        a2 = 1.0 - v * v
        γ = 1.0 / sqrt(a2)
        signxz = sign.(x_filtered .- z)

        kinetic[i] = γ
        self_potential[i] = γ * q * Φrz  #trapezoidal(sqrt(a2) * q * Φr .* δ, dx)
        dHpdt_analytic[i] = q * Ψrz * v * γ + q * γ^3 * Φrz * v * a +
                            q * γ * Πrz
        Π = Πr .+ Πs
        Ψ = Ψr .+ Ψs
        Φ = Φr .+ Φs

        Hfr[i, :] .= @. 0.5(Πr^2 + Ψr^2) + 0.5q * γ^3 * signxz * Ψr -
                        0.5q * signxz * Πr * v * γ^3 - q * signxz * Ψr * v^2 * γ^3 -
                        0.25 * q^2 * signxz^2 * v * a / (-1 + v^2)^3
        # Hfr[i, :] .-= Hfr[i, 1]
        regular_field_energy[i] = trapezoidal(Hfr[i, :], dx) +
                                  2x_range[2] * 0.25 * q^2 * v * a / (-1 + v^2)^3

        if i == 1
            dHfrdt_num[i, :] .= @. (Hfr[1, :] - Hfr[2, :]) / (dt)
        elseif i == nt
            dHfrdt_num[i, :] .= @. (Hfr[nt - 1, :] - Hfr[nt, :]) / (dt)
        else
            dHfrdt_num[i, :] .= @. (Hfr[i + 1, :] - Hfr[i - 1, :]) / (2dt)
        end
        dHfrdt_analytic12[i, :] .= @. 0.25 * (2q * Πr * δ * v^2 * γ^3 +
                                              2q * Ψr * δ * v * (-1 + 2v^2) * γ^3 -
                                              2q * signxz * Ψr * v * a * γ^5 -
                                              6q * signxz * Πr * v^2 * a * γ^5 -
                                              4q * signxz * Ψr * v^3 * a * γ^5 +
                                              4Ψr * dxΠr + 4Πr * dxΨr +
                                              2q * signxz * dxΠr * γ^5 -
                                              6q * signxz * v^2 * dxΠr * γ^5 +
                                              4q * signxz * v^4 * dxΠr * γ^5 -
                                              2q * signxz * v * dxΨr * γ^5 +
                                              2q * signxz * v^3 * dxΨr * γ^5 +
                                              1 * q^2 * signxz^2 * v * a / (-1 + v^2)^3)

        # dHfrdt_analytic22[i, :] .= @. (2q * Πr * δ * (-1 + v^2)^2 +
        #                                (signxz * (Πr + Ψr * v) - 2Φr * δ * v * a2) * a * q +
        #                                2Ψr * sqrt(a2) * dxΠr - 2Ψr * v^2 * sqrt(a2) * dxΠr -
        #                                q * a2 *
        #                                (v * (Ψr * δ - 2Ψr * a2) - signxz * dxΠr) +
        #                                2Πr * sqrt(a2) * dxΨr - 2Πr * v^2 * sqrt(a2) * dxΨr) /
        #                               denom

        # dHfrdt_analytic22[i, :] .= @. -0.5q * Ψr * δ * v * γ + Ψr * dxΠr + Πr * dxΨr +
        #                               (q * signxz * ((Πr + Ψr * v) * a + a2 * dxΠr)) / denom
        dHfrdt_analytic[i] = trapezoidal(dHfrdt_analytic12[i, :], dx)

        ## dHfdt + dTi0dx
        denom = 1 / (4 * a2^3.5)
        temp = @. -q * (q * signxz^2 * v * sqrt(a2) * a +
                        2signxz * v * a2 * (3Πr * v + Ψr * (1 + 2 * v^2)) *
                        a -
                        2 * (-1 + v^2)^2 *
                        (Πr * δ * (-2 + 3v^2) +
                         v * (Ψr * δ * (-1 + 2 * v^2) -
                              signxz * (v * dxΠr + dxΨr)))) / denom
        bl[i] = trapezoidal(temp, dx)
        Ti0 = -Πr .* Ψ
        dxTi0[i, :] .= deriv(Ti0, dx)

        tmp2 = @. -Ψ * dxΠr - Πr * (dxΨr + ∂xΨs)
        flux_analytic[i] = Ti0[end] - Ti0[begin] # trapezoidal(tmp2, dx)
        flux[i] = trapezoidal(tmp2, dx) # trapezoidal(dxTi0, dx)
    end

    Er = AkimaInterpolation(regular_field_energy, sim1.sol.t)
    dEdtr(t) = DataInterpolations.derivative(Er, t, 1)

    Ep = AkimaInterpolation(self_potential, sim1.sol.t)
    dEdtp(t) = DataInterpolations.derivative(Ep, t, 1)

    function time_derivative(y, h, ti)
        nt, nx = size(y)
        dy = zeros(eltype(y), nx)
        if ti == -1
            ti = nt
        end

        if ti == 1
            dy[:] .= @. (y[1, :] - y[2, :]) / h
        elseif ti == nt
            dy[:] .= @. (y[nt - 1, :] - y[nt, :]) / h
        else
            dy[:] .= @. (y[ti + 1, :] - y[ti - 1, :]) / (2h)
        end
        return dy
    end
    fig = Figure(; size=(2 * 253, 2 * 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Energy}")
    # Plot the data
    ti = 180

    # lines!(ax, x_filtered, Hfr[ti, :];
    #        label=L"\textrm{total_energy}")
    # # lines!(ax, sim1.sol.t, flux_analytic; label=L"\textrm{fluxes}")
    # ax.xlabel = L"x"
    # N = 1300
    # lines!(ax, sim1.sol.t[1:N], flux[1:N];
    #        label=L"\textrm{dxTi0}")
    # lines!(ax, sim1.sol.t[1:N], flux_analytic[1:N];
    #        label=L"\textrm{dxTi0 analytic}")

    # lines!(ax, sim1.sol.t, flux_analytic;
    #        label=L"\textrm{dxTi0 analytic}")

    # lines!(ax, sim1.sol.t, flux;
    #        label=L"\textrm{dxTi0}")
    # lines!(ax, x_filtered, time_derivative(Hfr, dx, ti);
    #        label=L"\textrm{dHdt}")
    # lines!(ax, x_filtered, dHfrdt_analytic12[ti, :];
    #        label=L"\textrm{dHdt an1}")
    # lines!(ax, x_filtered, dHfrdt_analytic22[ti, :];
    #        label=L"\textrm{dHdt an2}")

    # lines!(ax, sim1.sol.t, kinetic;
    #        label=L"\textrm{kinetic}")
    # lines!(ax, sim1.sol.t, self_potential;
    #        label=L"\textrm{self}")

    # lines!(ax, sim1.sol.t, regular_field_energy;
    #        label=L"\textrm{dHdt}")
    # lines!(ax, sim1.sol.t, 0.25dHfrdt_analytic;
    #        label=L"\textrm{dHdt an}")
    # lines!(ax, sim1.sol.t, dEdtr.(sim1.sol.t);
    #        label=L"\textrm{dHdt}")

    # lines!(ax, sim1.sol.t, bl;
    #        label=L"\textrm{bl}")

    # lines!(ax, sim1.sol.t, dHfrdt_analytic;
    #        label=L"\textrm{dHfdt an}")
    # lines!(ax, sim1.sol.t, dEdtr.(sim1.sol.t);
    #        label=L"\textrm{dHdt}")

    # lines!(ax, sim1.sol.t, dEdtp.(sim1.sol.t);
    #        label=L"\textrm{self-energy dt-interp}")
    # lines!(ax, sim1.sol.t, deriv(self_potential, dt);
    #        label=L"\textrm{self-energy dt}")
    # lines!(ax, sim1.sol.t, dHpdt_analytic;
    #        label=L"\textrm{analytic self-energy dt}")

    dtEtot_analytic = dHfrdt_analytic .+ dHpdt_analytic
    dtEtot_num = dEdtr.(sim1.sol.t) .+ dEdtp.(sim1.sol.t)
    # lines!(ax, sim1.sol.t, dtEtot_analytic;
    #        label=L"\textrm{analytic total energy dt}")
    # lines!(ax, sim1.sol.t, dtEtot_num;
    #        label=L"\textrm{numerical total energy}")
    lines!(ax, sim1.sol.t, dtEtot_analytic .- dtEtot_num;
           label=L"\textrm{analytic-num total energy dt}")

    ax.xlabel = L"x"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    return fig
end

end
