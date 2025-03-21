module ForcedPlots

using ..PlottingUtils
using ..ParticleMotion
using CairoMakie
using LinearAlgebra

@inline name(pre, sim, post=nothing) = name_forced(pre, sim, post)

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

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

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
    savename = name("forced_motion/pifield_res_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_resolutions(i, sim1, sim2, sim4, F=8)
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
    savename = name("forced_motion/psifield_res_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield_diff(i, sim1, sim2, sim4, F=8)
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

    diff1 = y1_filtered .- y2_filtered
    diff2 = y2_filtered .- y4_filtered
    # Compute y-limits with padding
    y_min = minimum([minimum(diff1), minimum(diff2)])
    y_max = maximum([maximum(diff1), maximum(diff2)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, diff1; label=L"\textrm{low-mid}")
    lines!(ax, x_filtered, diff2; label=L"\textrm{mid-high}")
    # lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
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
    savename = name("forced_motion/pifield_diff_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_diff(i, sim1, sim2, sim4, F=8)
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
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, diff1; label=L"\textrm{low-mid}")
    lines!(ax, x_filtered, 2^2 * diff2; label=L"\textrm{mid-high}", linestyle=:dot)
    # lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = name("forced_motion/psifield_diff_forced", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield_convfac(sim1, sim2, sim4)
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
              limits=(-1, sim1.sol.t[end], 1, 2))
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
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_convfac(sim1, sim2, sim4)
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
              limits=(-1, sim1.sol.t[end], 1, 3))
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
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield_convord(sim1, sim2, sim4)
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

    savename = name("forced_motion/pifield_convord_forced", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_convord(sim1, sim2, sim4)
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

    savename = name("forced_motion/psifield_convord_forced", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_convergence_panel(sim1, sim2, sim4)
    nt = length(sim1.sol.t)
    fig1 = plot_pifield_resolutions_forced(nt, sim1, sim2, sim4)  # Your existing function
    fig2 = plot_psifield_resolutions_forced(nt, sim1, sim2, sim4)  # Replace with your second function
    fig3 = plot_pifield_convord_forced(sim1, sim2, sim4)  # Replace with your third function
    fig4 = plot_psifield_convord_forced(sim1, sim2, sim4)  # Replace with your fourth function

    fig = Figure(; size=(600, 500))

    # Add supertitle
    supertitle = Label(fig[1, 1:2], L"v_{max}=%$(sim1.params.vmax)";
                       fontsize=20, font=:bold, halign=:center)

    # Function to recreate Axis content
    function recreate_axis!(fig, src_fig, row, col)
        src_ax = contents(src_fig[1, 1])[1]  # Original Axis

        # Create new Axis with same properties
        new_ax = Axis(fig[row, col];
                      xlabel=src_ax.xlabel[],
                      ylabel=src_ax.ylabel[],
                      title=src_ax.title[],
                      limits=src_ax.limits[])

        # Copy plot objects
        for plot in src_ax.scene.plots
            if plot isa MakieCore.Scatter
                label = haskey(plot.attributes, :label) ? plot.label[] : nothing
                scatter!(new_ax, plot[1][], plot[2][]; label=label,
                         markersize=plot.markersize[])
            elseif plot isa MakieCore.Lines
                data = plot[1][]  # Single Observable
                label = haskey(plot.attributes, :label) ? plot.label[] : nothing
                if isa(data, Vector{Point{2,Float64}})  # Regular Lines
                    x = [p[1] for p in data]
                    y = [p[2] for p in data]
                    lines!(new_ax, x, y; label=label, linewidth=plot.linewidth[])
                elseif isa(data, Vector{Float64})  # VLines misclassified as Lines
                    @show "ASDFA"
                    vlines!(new_ax, data; color=plot.color[], linestyle=plot.linestyle[],
                            linewidth=plot.linewidth[])
                end
            end
        end

        # Add Legend if it exists, place in the next row under the column
        leg_content = contents(src_fig[2, 1])
        if !isempty(leg_content)
            Legend(fig[row + 1, col], new_ax; orientation=:horizontal, tellwidth=false,
                   tellheight=true)
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
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

end
