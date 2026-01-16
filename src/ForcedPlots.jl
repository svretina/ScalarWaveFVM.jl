module ForcedPlots

using ..PlottingUtils
using ..ParticleMotion
using ..ScalarField
using CairoMakie
using LinearAlgebra
using DataInterpolations

using ..RetardedPotentials
using ..Run

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

## Particle Motion
function get_velocity(sim)
    if hasproperty(sim.sol.prob.p, :equation)
        return sim.sol.prob.p.equation.velocity
    end
    return 1.0
end

function plot_particle_position(sim; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim.sol.t)
    # Create figure and axis

    L = sim.params.L
    c = get_velocity(sim)
    if iszero(sim.sol.t[end])
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * sim.sol.t[end], c * sim.sol.t[end])
    end
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
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
    savename = simple_name("particle_position")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Fields
function plot_pifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    t = sim1.sol.t[i]
    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    y1_filtered = sim1.sol.u[i][mask, 1]

    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)
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
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("pifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t = t1

    # Create figure and axis

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("psifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_phifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t = t1
    q = sim1.params.q
    # Create figure and axis
    xp, vp, ap = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A,
                                           sim1.params.ω)

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]
    y1_filtered = sim1.sol.u[i][mask1, 3] .- 0ScalarField.Πs.(x_filtered, xp, vp, q)

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered;
           label=L"\textrm{low}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("phifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Total field
function plot_total_pifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis

    L = sim1.params.L
    c = get_velocity(sim1)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("pifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_psifield(i, sim1; F=2, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("psifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_phifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis
    c = get_velocity(sim1)
    L = sim1.params.L
    @show L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("phifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Resolutions  
function plot_total_pifield_resolutions(i, sim1, sim2, sim4; F=8, paper=false)
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

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / F, L / F)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("pifield_res_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_resolutions(i, sim1, sim2, sim4; F=8, paper=false)
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

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / F, L / F)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("psifield_res_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_resolutions(i, sim1, sim2, sim4; F=8, paper=false)
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

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("pifield_res_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Differences
function plot_pifield_diff(i, sim1, sim2, sim4; paper=false)
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
    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    @show x_range

    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("pifield_diff_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_diff(i, sim1, sim2, sim4; paper=false)
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
    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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
    savename = simple_name("psifield_diff_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Convergence Factor
function plot_pifield_convfac(sim1, sim2, sim4; paper=false)
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
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
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

    savename = simple_name("pifield_convfac_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_convfac(sim1, sim2, sim4; paper=false)
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
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
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

    savename = simple_name("psifield_convfac_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Convergence Order
function plot_pifield_self_convergence_order(sim1, sim2, sim4; paper=false)
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
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
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

    savename = simple_name("pifield_convord_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_self_convergence_order(sim1, sim2, sim4; paper=false)
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
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
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

    savename = simple_name("psifield_convord_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_convergence_panel(sim1, sim2, sim4; paper=false)
    nt = length(sim1.sol.t)
    fig1 = plot_pifield_resolutions(nt, sim1, sim2, sim4; paper=paper)  # Your existing function
    fig2 = plot_psifield_resolutions(nt, sim1, sim2, sim4; paper=paper)  # Replace with your second function
    fig3 = plot_pifield_self_convergence_order(sim1, sim2, sim4; paper=paper)  # Replace with your third function
    fig4 = plot_psifield_self_convergence_order(sim1, sim2, sim4; paper=paper)  # Replace with your fourth function

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

    savename = simple_name("convergence_pane_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

### Comparison with Retarded
function plot_pifield_vs_retarded(i, sim; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim.sol.t[i]

    dx = sim.params.dx
    L = sim.params.L
    q = sim.params.q
    A = sim.params.A
    ω = sim.params.ω
    z0 = sim.params.z0
    z0 = sim.params.z0

    c = get_velocity(sim)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    x = sim.params.x[mask]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Pis .+ sim.sol.u[i][mask, 1]
    else
        y_sim = sim.sol.u[i][mask, 1] #.- 0Pis0
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret.(x, Ref(sim.params), Ref(t))
    else
        y_lw = Π_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Πs.(x, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x, y_sim; label="Simulation")
    lines!(ax, x, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_pifield_resolutions_vs_retarded(i, sim1, sim2; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim1.sol.t[i]

    dx1 = sim1.params.dx
    dx2 = sim2.params.dx
    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    c = get_velocity(sim1)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])

    x1 = sim1.params.x[mask1]
    x2 = sim2.params.x[mask2]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x1, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Pis .+ sim1.sol.u[i][mask1, 1]
    else
        # xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        # Pis0 = ScalarField.Πs.(x, Ref(xp0), Ref(vp0), Ref(q))
        y_sim1 = sim1.sol.u[i][mask1, 1] #.- 0Pis0
        y_sim2 = sim2.sol.u[2i - 1][mask2, 1]
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret.(x1, Ref(sim1.params), Ref(t))
    else
        y_lw = Π_ret.(x1, Ref(sim1.params), Ref(t)) .- ScalarField.Πs.(x1, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim1)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim1)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_r",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x1, y_sim1; label="dx=$dx1")
    lines!(ax, x2, y_sim2; label="dx=$dx2")
    lines!(ax, x1, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("pifield_retarded_resolutions", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

function plot_psifield_resolutions_vs_retarded(i, sim1, sim2; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim1.sol.t[i]

    dx1 = sim1.params.dx
    dx2 = sim2.params.dx
    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    c = get_velocity(sim1)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])

    x1 = sim1.params.x[mask1]
    x2 = sim2.params.x[mask2]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Psis1 = ScalarField.Ψs.(x1, Ref(z), Ref(v), Ref(q))
    Psis2 = ScalarField.Ψs.(x2, Ref(z), Ref(v), Ref(q))

    if total
        y_sim1 = Psis1 .+ sim1.sol.u[i][mask1, 2]
        y_sim2 = Psis2 .+ sim2.sol.u[2i - 1][mask2, 2]
    else
        xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        Psis0_1 = ScalarField.Ψs.(x1, Ref(xp0), Ref(vp0), Ref(q))
        Psis0_2 = ScalarField.Ψs.(x2, Ref(xp0), Ref(vp0), Ref(q))
        y_sim1 = Psis1 .+ sim1.sol.u[i][mask1, 2] .- Psis0_1
        y_sim2 = Psis2 .+ sim2.sol.u[2i - 1][mask2, 2] .- Psis0_2
    end

    ## Lienard-Wiechert

    if total
        y_lw = Ψ_ret.(x1, Ref(sim1.params), Ref(t))
    else
        y_lw = Ψ_ret.(x1, Ref(sim1.params), Ref(t)) .- ScalarField.Ψs.(x1, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim1)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim1)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_r",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x1, y_sim1; label="dx=$dx1")
    lines!(ax, x2, y_sim2; label="dx=$dx2")
    lines!(ax, x1, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("psifield_retarded_resolutions", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

function plot_psifield_vs_retarded(i, sim; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim.sol.t[i]

    dx = sim.params.dx
    L = sim.params.L
    q = sim.params.q
    A = sim.params.A
    ω = sim.params.ω
    z0 = sim.params.z0

    x_range = (-L / 2, L / 2)
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    x = sim.params.x[mask]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Psis = ScalarField.Ψs.(x, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Psis .+ sim.sol.u[i][mask, 2]
    else
        xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        Psis0 = ScalarField.Ψs.(x, Ref(xp0), Ref(vp0), Ref(q))
        y_sim = Psis .+ sim.sol.u[i][mask, 2] .- Psis0
    end

    ## Lienard-Wiechert

    if total
        y_lw = Ψ_ret.(x, Ref(sim.params), Ref(t))
    else
        y_lw = Ψ_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Ψs.(x, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x, y_sim; label="Simulation")
    lines!(ax, x, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_pifield_diff_vs_retarded(i, sim1, sim2; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    dx1 = sim1.params.dx
    dx2 = sim2.params.dx

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    @assert t1 == t2
    t = t1

    c = get_velocity(sim1)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x = sim1.params.x[mask1]

    c = get_velocity(sim1)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))
    if total
        y_sim1 = Pis .+ sim1.sol.u[i][mask1, 1]
        y_sim2 = Pis .+ mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    else
        y_sim1 = sim1.sol.u[i][mask1, 1]
        y_sim2 = mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret.(x, Ref(sim1.params), Ref(t))
    else
        y_lw = Π_ret.(x, Ref(sim1.params), Ref(t)) .- ScalarField.Πs.(x, z, v, q)
    end

    diff1 = y_sim1 .- y_lw
    diff2 = y_sim2 .- y_lw

    convfac = diff1 ./ diff2
    @show norm(convfac)
    @show norm(diff1) / norm(diff2)
    @show log2(norm(diff1, 1) / norm(diff2, 1))
    # ymin = minimum(convfac)
    # ymax = maximum(convfac)

    # padding = round(0.1(ymax - ymin); digits=3)

    # y_limits = round.((ymin - padding, ymax + padding), digits=3)
    #y_limits = (0.0, 5.0)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              # yscale=log10,
              # ylabel=L"",
              title=L"\textrm{Differences \Pi for t=%$(round(t, digits=3)) }")
    # ,
    # limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    # lines!(ax, x, convfac)
    lines!(ax, x, diff1)
    lines!(ax, x, diff2)

    lines!(ax, x, y_lw; linestyle=:dash)
    lines!(ax, x, y_sim1)
    lines!(ax, x, y_sim2)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_pifield_error_convergence(i, Nosc, dx0, cfl, vmax, q; total=false,
                                        paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q)

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    nt = length(sim1.sol.t)
    convfac = zeros(nt - 1)

    hs = [dx0, dx0 / 2, dx0 / 4, dx0 / 8, dx0 / 16]
    errors = zeros(length(hs))

    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q)
    sim8 = Run.forced_motion(Nosc, dx0 / 8, cfl, vmax, q)
    sim16 = Run.forced_motion(Nosc, dx0 / 16, cfl, vmax, q)

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    t8 = sim8.sol.t[8i - 7]
    t16 = sim16.sol.t[16i - 15]

    @assert t1 == t2
    @assert t2 == t4
    @assert t4 == t8
    @assert t4 == t8
    @assert t8 == t16

    c = get_velocity(sim1)
    if iszero(t1)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t1, c * t1)
    end

    # ### Method 1 (averaging)
    # mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # x1 = sim1.params.x[mask1]

    # z, v, _ = ParticleMotion.oscillator(t1, z0, A, ω)
    # y_lw = Π_ret.(x1, Ref(sim1.params), Ref(t1)) .- ScalarField.Πs.(x1, z, v, q)

    # y_sim1 = sim1.sol.u[i][mask1, 1]
    # y_sim2 = mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    # y_sim4 = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))[mask1]
    # y_sim8 = mymean(mymean(mymean(sim8.sol.u[8i - 7][:, 1])))[mask1]
    # y_sim16 = mymean(mymean(mymean(mymean(sim16.sol.u[16i - 15][:, 1]))))[mask1]

    # mode = 1
    # errors[1] = norm(abs.(y_lw .- y_sim1), mode)
    # errors[2] = norm(abs.(y_lw .- y_sim2), mode)
    # errors[3] = norm(abs.(y_lw .- y_sim4), mode)
    # errors[4] = norm(abs.(y_lw .- y_sim8), mode)
    # errors[5] = norm(abs.(y_lw .- y_sim16), mode)

    ### Method 2
    z, v, _ = ParticleMotion.oscillator(t1, z0, A, ω)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1 = sim1.params.x[mask1]

    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    x2 = sim2.params.x[mask2]

    mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    x4 = sim4.params.x[mask4]

    mask8 = (sim8.params.x .>= x_range[1]) .& (sim8.params.x .<= x_range[2])
    x8 = sim8.params.x[mask8]

    mask16 = (sim16.params.x .>= x_range[1]) .& (sim16.params.x .<= x_range[2])
    x16 = sim16.params.x[mask16]

    y_lw1 = Π_ret.(x1, Ref(sim1.params), Ref(t1)) .- ScalarField.Πs.(x1, z, v, q)
    y_lw2 = Π_ret.(x2, Ref(sim2.params), Ref(t1)) .- ScalarField.Πs.(x2, z, v, q)
    y_lw4 = Π_ret.(x4, Ref(sim4.params), Ref(t1)) .- ScalarField.Πs.(x4, z, v, q)
    y_lw8 = Π_ret.(x8, Ref(sim8.params), Ref(t1)) .- ScalarField.Πs.(x8, z, v, q)
    y_lw16 = Π_ret.(x16, Ref(sim16.params), Ref(t1)) .- ScalarField.Πs.(x16, z, v, q)

    y_sim1 = sim1.sol.u[i][mask1, 1]
    y_sim2 = sim2.sol.u[2i - 1][mask2, 1]
    y_sim4 = sim4.sol.u[4i - 3][mask4, 1]
    y_sim8 = sim8.sol.u[8i - 7][mask8, 1]
    y_sim16 = sim16.sol.u[16i - 15][mask16, 1]

    mode = Inf
    errors[1] = norm(abs.(y_lw1 .- y_sim1), mode)
    errors[2] = norm(abs.(y_lw2 .- y_sim2)[begin:2:end], mode)
    errors[3] = norm(abs.(y_lw4 .- y_sim4)[begin:4:end], mode)
    errors[4] = norm(abs.(y_lw8 .- y_sim8)[begin:8:end], mode)
    errors[5] = norm(abs.(y_lw16 .- y_sim16)[begin:16:end], mode)

    @show hs
    @show errors
    @show t1

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"h",
              # ylabel=L"\Pi",
              yscale=log10,
              xscale=log10,
              xreversed=true,
              title=L"\Vert |\Pi_{\textrm{ret}}-\Pi_\textrm{h}| \Vert_{\infty}")

    lines!(ax, hs, errors)
    scatter!(ax, hs, errors)

    xref = hs[1:3]
    yref1 = hs[1:3] .^ 1
    shift1 = (errors[1] / yref1[1]) * 1.3
    yref1 .*= shift1

    lines!(ax, xref, yref1; linestyle=:dash, label=L"h^1")
    # axislegend(ax; position=:rt)

    # leg = Legend(fig[1, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)
    leg = Legend(fig[1, 1], ax; # Note: Same grid position as the Axis!
                 halign=:right,
                 valign=:top,
                 tellwidth=false,      # Important: Don't shrink the axis to fit the legend
                 tellheight=false,     # Important: Don't shrink the axis to fit the legend
                 margin=(10, 10, 10, 10), # Optional: Add padding from the axis lines
                 framevisible=false,
                 labelsize=12)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("pifield_retarded_convergence", t1)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

function plot_psifield_error_convergence(i, Nosc, dx0, cfl, vmax, q; total=false,
                                         paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q)

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    nt = length(sim1.sol.t)
    convfac = zeros(nt - 1)

    hs = [dx0, dx0 / 2, dx0 / 4, dx0 / 8, dx0 / 16]
    errors = zeros(length(hs))

    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q)
    sim8 = Run.forced_motion(Nosc, dx0 / 8, cfl, vmax, q)
    sim16 = Run.forced_motion(Nosc, dx0 / 16, cfl, vmax, q)

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    t8 = sim8.sol.t[8i - 7]
    t16 = sim16.sol.t[16i - 15]

    @assert t1 == t2
    @assert t2 == t4
    @assert t4 == t8
    @assert t4 == t8
    @assert t8 == t16

    c = get_velocity(sim1)
    if iszero(t1)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t1, c * t1)
    end

    ### Method 2
    z, v, _ = ParticleMotion.oscillator(t1, z0, A, ω)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1 = sim1.params.x[mask1]

    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    x2 = sim2.params.x[mask2]

    mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    x4 = sim4.params.x[mask4]

    mask8 = (sim8.params.x .>= x_range[1]) .& (sim8.params.x .<= x_range[2])
    x8 = sim8.params.x[mask8]

    mask16 = (sim16.params.x .>= x_range[1]) .& (sim16.params.x .<= x_range[2])
    x16 = sim16.params.x[mask16]

    y_lw1 = Ψ_ret.(x1, Ref(sim1.params), Ref(t1)) .- ScalarField.Ψs.(x1, z, v, q)
    y_lw2 = Ψ_ret.(x2, Ref(sim2.params), Ref(t1)) .- ScalarField.Ψs.(x2, z, v, q)
    y_lw4 = Ψ_ret.(x4, Ref(sim4.params), Ref(t1)) .- ScalarField.Ψs.(x4, z, v, q)
    y_lw8 = Ψ_ret.(x8, Ref(sim8.params), Ref(t1)) .- ScalarField.Ψs.(x8, z, v, q)
    y_lw16 = Ψ_ret.(x16, Ref(sim16.params), Ref(t1)) .- ScalarField.Ψs.(x16, z, v, q)

    y_sim1 = sim1.sol.u[i][mask1, 2]
    y_sim2 = sim2.sol.u[2i - 1][mask2, 2]
    y_sim4 = sim4.sol.u[4i - 3][mask4, 2]
    y_sim8 = sim8.sol.u[8i - 7][mask8, 2]
    y_sim16 = sim16.sol.u[16i - 15][mask16, 2]

    mode = Inf
    errors[1] = norm(abs.(y_lw1 .- y_sim1), mode)
    errors[2] = norm(abs.(y_lw2 .- y_sim2)[begin:2:end], mode)
    errors[3] = norm(abs.(y_lw4 .- y_sim4)[begin:4:end], mode)
    errors[4] = norm(abs.(y_lw8 .- y_sim8)[begin:8:end], mode)
    errors[5] = norm(abs.(y_lw16 .- y_sim16)[begin:16:end], mode)

    @show hs
    @show errors
    @show t1

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"h",
              # ylabel=L"\Psi",
              yscale=log10,
              xscale=log10,
              xreversed=true,
              title=L"\Vert |\Psi_{\textrm{ret}}-\Psi_\textrm{h}| \Vert_{\infty}")

    lines!(ax, hs, errors)
    scatter!(ax, hs, errors)

    xref = hs[1:3]
    yref1 = hs[1:3] .^ 1
    shift1 = (errors[1] / yref1[1]) * 1.3
    yref1 .*= shift1

    lines!(ax, xref, yref1; linestyle=:dash, label=L"h^1")
    # axislegend(ax; position=:rt)

    leg = Legend(fig[1, 1], ax; # Note: Same grid position as the Axis!
                 halign=:right,
                 valign=:top,
                 tellwidth=false,      # Important: Don't shrink the axis to fit the legend
                 tellheight=false,     # Important: Don't shrink the axis to fit the legend
                 margin=(10, 10, 10, 10), # Optional: Add padding from the axis lines
                 framevisible=false,
                 labelsize=12)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("psifield_retarded_convergence", t1)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end
end
