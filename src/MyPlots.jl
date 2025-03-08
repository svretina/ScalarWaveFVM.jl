module MyPlots

using LaTeXStrings
using ..ParticleMotion
using ..Run
using LinearAlgebra
using MakieCore
using Makie
using CairoMakie

const PROJ_PATH = pkgdir(@__MODULE__)
const FIG_PATH = joinpath(PROJ_PATH, "figures")

function mytheme_aps()
    return Theme(
                 # Axis attributes
                 ;
                 Axis=Attributes(; spinewidth=1.1,
                                 xgridvisible=false,
                                 xlabelpadding=-2,
                                 xlabelsize=10,
                                 xminortickalign=1,
                                 xminorticks=IntervalsBetween(5, true),
                                 xminorticksize=3,
                                 xminorticksvisible=true,
                                 xminortickwidth=0.75,
                                 xtickalign=1,
                                 xticklabelsize=8,
                                 xticksize=5,
                                 xticksmirrored=true,
                                 xtickwidth=0.8,
                                 ygridvisible=false,
                                 ylabelpadding=2,
                                 ylabelsize=10,
                                 yminortickalign=1,
                                 yminorticks=IntervalsBetween(5, true),
                                 yminorticksize=3,
                                 yminorticksvisible=true,
                                 yminortickwidth=0.75,
                                 ytickalign=1,
                                 yticklabelsize=8,
                                 yticksize=5,
                                 yticksmirrored=true,
                                 ytickwidth=0.8,
                                 xticklabelfont="cmr10",  # Upright Computer Modern
                                 yticklabelfont="cmr10",  # Upright Computer Modern
                                 xticklabelstyle=Attributes(; italic=false),
                                 yticklabelstyle=Attributes(; italic=false)),
                 # General figure settings
                 colgap=8,
                 figure_padding=0,
                 rowgap=8,
                 size=(243, 165),
                 # Colorbar attributes
                 Colorbar=Attributes(; labelpadding=2,
                                     labelsize=10,
                                     minortickalign=1,
                                     minorticksize=3,
                                     minorticksvisible=true,
                                     minortickwidth=0.75,
                                     size=8,
                                     spinewidth=1.1,
                                     tickalign=1,
                                     ticklabelpad=2,
                                     ticklabelsize=8,
                                     ticksize=5,
                                     tickwidth=0.8),
                 fonts=Attributes(; bold="NewComputerModern10 Bold",
                                  bold_italic="NewComputerModern10 Bold Italic",
                                  italic="NewComputerModern10 Italic",
                                  regular="NewComputerModern Math Regular"),
                 # fonts=Attributes(; bold="ComputerModern Bold",
                 #                  bold_italic="ComputerModern Bold Italic",
                 #                  italic="ComputerModern Italic",
                 #                  regular="ComputerModern Math Regular"),
                 # Legend attributes
                 Legend=Attributes(; colgap=4,
                                   framecolor=(:grey, 0.5),
                                   framevisible=false,
                                   labelsize=7.5,
                                   margin=(0, 0, 0, 0),
                                   nbanks=1,
                                   padding=(2, 2, 2, 2),
                                   rowgap=-10,
                                   #labelfont="cmr10"
                                   ),
                 # Lines attributes
                 Lines=Attributes(;
                                  cycle=Cycle([[:color] => :color, [:marker] => :marker],
                                              true)),
                 # Scatter attributes
                 Scatter=Attributes(;
                                    cycle=Cycle([[:color] => :color, [:marker] => :marker],
                                                true),
                                    markersize=7,
                                    strokewidth=0),
                 markersize=7,
                 # Palette attributes
                 palette=Attributes(;
                                    color=[RGBAf(0.298039, 0.447059, 0.690196, 1.0),
                                           RGBAf(0.866667, 0.517647, 0.321569, 1.0),
                                           RGBAf(0.333333, 0.658824, 0.407843, 1.0),
                                           RGBAf(0.768627, 0.305882, 0.321569, 1.0),
                                           RGBAf(0.505882, 0.447059, 0.701961, 1.0),
                                           RGBAf(0.576471, 0.470588, 0.376471, 1.0),
                                           RGBAf(0.854902, 0.545098, 0.764706, 1.0),
                                           RGBAf(0.54902, 0.54902, 0.54902, 1.0),
                                           RGBAf(0.8, 0.72549, 0.454902, 1.0),
                                           RGBAf(0.392157, 0.709804, 0.803922, 1.0)],
                                    linestyle=[nothing, :dash, :dot, :dashdot, :dashdotdot],
                                    marker=[:circle, :rect, :dtriangle, :utriangle, :cross,
                                            :diamond, :ltriangle, :rtriangle, :pentagon,
                                            :xcross, :hexagon],
                                    markercolor=[RGBAf(0.298039, 0.447059, 0.690196, 1.0),
                                                 RGBAf(0.866667, 0.517647, 0.321569, 1.0),
                                                 RGBAf(0.333333, 0.658824, 0.407843, 1.0),
                                                 RGBAf(0.768627, 0.305882, 0.321569, 1.0),
                                                 RGBAf(0.505882, 0.447059, 0.701961, 1.0),
                                                 RGBAf(0.576471, 0.470588, 0.376471, 1.0),
                                                 RGBAf(0.854902, 0.545098, 0.764706, 1.0),
                                                 RGBAf(0.54902, 0.54902, 0.54902, 1.0),
                                                 RGBAf(0.8, 0.72549, 0.454902, 1.0),
                                                 RGBAf(0.392157, 0.709804, 0.803922, 1.0)],
                                    patchcolor=[RGBAf(0.298039, 0.447059, 0.690196, 1.0),
                                                RGBAf(0.866667, 0.517647, 0.321569, 1.0),
                                                RGBAf(0.333333, 0.658824, 0.407843, 1.0),
                                                RGBAf(0.768627, 0.305882, 0.321569, 1.0),
                                                RGBAf(0.505882, 0.447059, 0.701961, 1.0),
                                                RGBAf(0.576471, 0.470588, 0.376471, 1.0),
                                                RGBAf(0.854902, 0.545098, 0.764706, 1.0),
                                                RGBAf(0.54902, 0.54902, 0.54902, 1.0),
                                                RGBAf(0.8, 0.72549, 0.454902, 1.0),
                                                RGBAf(0.392157, 0.709804, 0.803922, 1.0)]),
                 # PolarAxis attributes
                 PolarAxis=Attributes(; spinewidth=1.1))
end

function mymean(arr)
    arr2 = zeros(eltype(arr), div(length(arr), 2))
    for i in eachindex(arr2)
        j = 2i - 1
        arr2[i] = 0.5(arr[j] + arr[j + 1])
    end
    return arr2
end

function plot_pifield_resolutions_forced(i, sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_res_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_res_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_psifield_resolutions_forced(i, sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_res_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_res_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_pifield_diff_forced(i, sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_diff_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_diff_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_psifield_diff_forced(i, sim1, sim2, sim4)
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
    lines!(ax, x_filtered, diff1; label=L"\textrm{low-mid}")
    lines!(ax, x_filtered, diff2; label=L"\textrm{mid-high}")
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
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_diff_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_diff_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_pifield_convfac_forced(sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_convfac_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_convfac_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).svg"),
         fig)
    return fig
end

function plot_psifield_convfac_forced(sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_convfac_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_convfac_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).svg"),
         fig)
    return fig
end

function plot_pifield_convord_forced(sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_convord_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_convord_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).svg"),
         fig)
    return fig
end

function plot_psifield_convord_forced(sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_convord_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/psifield_convord_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).svg"),
         fig)
    return fig
end

function plot_convergence_panel_forced(sim1, sim2, sim4)
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
    save(joinpath(FIG_PATH,
                  "forced_motion/convergence_panel_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/convergence_panel_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).svg"),
         fig)
    return fig
end

#### Plots for interacting particles

function plot_pifield_interacting(i, sim1)
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
    save(joinpath(FIG_PATH,
                  "interacting_particles/pifield_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "interacting_particles/pifield_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_psifield_interacting(i, sim1)
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
    save(joinpath(FIG_PATH,
                  "interacting_particles/psifield_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "interacting_particles/psifield_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_particles_position(sim1; invert=true)
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
    save(joinpath(FIG_PATH,
                  "interacting_particles/particles_position_potential_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "interacting_particles/particles_position_potential_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).svg"),
         fig)
    return fig
end

function plot_particles_velocity(sim1)
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
    save(joinpath(FIG_PATH,
                  "interacting_particles/particles_velocity_potential_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "interacting_particles/particles_velocity_potential_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl).svg"),
         fig)
    return fig
end

###########
#################3
function plot_pifield_resolutions_potential(i, sim1, sim2, sim4)
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

    # xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

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
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_res_potential_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "forced_motion/pifield_res_potential_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_psifield_resolutions_potential(i, sim1, sim2, sim4)
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

    # xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

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
    save(joinpath(FIG_PATH,
                  "potential/psifield_res_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).pdf"),
         fig)
    save(joinpath(FIG_PATH,
                  "potential/psifield_res_forced_dx=$(sim1.params.dx)_cfl=$(sim1.params.cfl)_i=$(i).svg"),
         fig)
    return fig
end

function plot_field_resolutions(i, sim1, sim2, sim4)
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1
    p = scatter(sim1.params.x, sim1.sol.u[i].x[1][:, 1]; label="low res")
    scatter!(p, sim1.params.x, mymean(sim2.sol.u[2i - 1].x[1][:, 1]); label="mid res")
    scatter!(p, sim1.params.x, mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]));
             label="high res")
    xaxis!(p, "x")
    yaxis!(p, "Π")
    title!(p, "t[$(i)]=$(round(t, digits=3))")
    savefig(p, "pifield_res_i=$(i).pdf")
    return p
end

function plot_pifield(i, sim)
    t = sim.sol.t[i]
    p = scatter(sim.params.x, sim.sol.u[i].x[1][:, 1]; label=false)
    xaxis!(p, "x")
    yaxis!(p, "Π")
    title!(p, "t[$(i)]=$(round(t, digits=3))")
    savefig(p, "pifield_i=$(i).pdf")
    return p
end

function plot_psifield(i, sim)
    t = sim.sol.t[i]
    p = scatter(sim.params.x, sim.sol.u[i].x[1][:, 2]; label=false)
    xaxis!(p, "x")
    yaxis!(p, "Ψ")
    title!(p, "t[$(i)]=$(round(t, digits=3))")
    savefig(p, "psifield_i=$(i).pdf")
    return p
end

function plot_field_errors(i, sim1, sim2, sim4, scale=1.0)
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1
    p = plot()
    scatter!(p, sim1.params.x,
             sim1.sol.u[i].x[1][:, 1] .- mymean(sim2.sol.u[2i - 1].x[1][:, 1]);
             label="low-mid")
    scatter!(p, sim1.params.x,
             scale .* (mymean(sim2.sol.u[2i - 1].x[1][:, 1]) .-
                       mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1])));
             label="mid-high")
    return p
end

# function make_gif_errors(scale=1.0)
#     nt = length(sol1.t)
#     p = plot()
#     ylims!(p, -0.1, 1.1)
#     anim = @animate for i in 1:10:nt
#         plot_errors(i, scale)
#     end
#     gif(anim, "scaled_errors_fps15.gif"; fps=15)
# end

# function make_gif(x, sol)
#     nt = length(sol.t)
#     p = plot()
#     ylims!(p, -0.1, 1.1)
#     anim = @animate for i in 1:10:nt
#         scatter(p, x, sol[:, 1, i])
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

# function make_gif_with_forced_particle(x, sol, params)
#     nt = length(sol.t)

#     x10 = params.x1
#     x20 = params.x2
#     L = params.L
#     direction1 = params.direction1
#     direction2 = params.direction2
#     vmax = params.vmax
#     f0 = params.f0
#     A = params.A
#     anim = @animate for i in 1:10:nt
#         p = plot()
#         ylims!(p, -3.1, 5.5)
#         t = sol.t[i]
#         x1, v1, a1 = ParticleMotion.oscillator2(t, x10, vmax, A, direction1)
#         x2, v2, a2 = ParticleMotion.oscillator2(t, x20, vmax, A, direction2)
#         plot!(p, x, sol[:, 1, i]; label=false)
#         vline!(p, [x1]; label=false)
#         vline!(p, [x2]; label=false)
#         title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

# function make_gif_with_particle(xs, sol)
#     nt = length(sol.t)
#     ymins = zeros(nt)
#     ymaxs = zeros(nt)
#     for i in 1:nt
#         ymins[i] = min(sol.u[i].x[1][:, 1]...)
#         ymaxs[i] = max(sol.u[i].x[1][:, 1]...)
#     end
#     ymin = min(ymins...)
#     ymax = max(ymaxs...)
#     p = plot()
#     anim = @animate for i in 1:10:nt
#         p = plot()
#         # ylims!(p, 1.1ymin, 1.1ymax)
#         t = sol.t[i]
#         x1 = sol.u[i].x[2][2]

#         plot!(p, xs, sol.u[i].x[1][:, 1]; label=false)
#         vline!(p, [x1]; label=false)
#         title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

# function make_gif_with_particles(sim)
#     nt = length(sim.sol.t)
#     # ymins = zeros(nt)
#     # ymaxs = zeros(nt)
#     # for i in 1:nt
#     #     ymins[i] = min(sim.sol.u[i].x[1][:, 1]...)
#     #     ymaxs[i] = max(sim.sol.u[i].x[1][:, 1]...)
#     # end
#     # ymin = min(ymins...)
#     # ymax = max(ymaxs...)
#     anim = @animate for i in 1:10:nt
#         p = plot()
#         # ylims!(p, 1.1ymin, 1.1ymax)
#         t = sim.sol.t[i]
#         x1 = sim.sol.u[i].x[2][2]
#         x2 = sim.sol.u[i].x[2][5]

#         plot!(p, sim.params.x, sim.sol.u[i].x[1][:, 1]; label=false)
#         vline!(p, [x1]; label=false)
#         vline!(p, [x2]; label=false)
#         title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

function plot_position(sim, existing_fig=nothing; invert=false)
    if existing_fig === nothing
        # Create a new figure with 1 row and 2 columns if no figure is provided
        p = plot(; layout=(1, 2))
    else
        # Use the existing figure
        p = existing_fig
        # Get the layout dimensions
        rows, cols = size(p.layout)
        # Check if there's a second row available
        if rows < 2
            error("The existing figure must have at least 2 rows to plot in the second row")
        end
    end
    if existing_fig === nothing
        idx1, idx2 = 1, 2
    else
        idx1 = cols + 1      # Start of second row
        idx2 = cols + 2      # Second position in second row
        # Check if these indices exist in the layout
        if idx2 > rows * cols
            error("Not enough subplots in the existing figure's layout for second row plotting")
        end
    end
    nt = length(sim.sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    for i in 1:nt
        x1[i] = sim.sol.u[i].x[2][2]
        x2[i] = sim.sol.u[i].x[2][5]
    end
    if invert
        plot!(p[idx1], x1, sim.sol.t; title="particle 1", xlabel="x", ylabel="t",
              label=false)
        plot!(p[idx2], x2, sim.sol.t; title="particle 2", xlabel="x", ylabel="",
              label=false)
    else
        plot!(p[idx1], sim.sol.t, x1; title="particle 1", xlabel="t", ylabel="x",
              label=false)
        plot!(p[idx2], sim.sol.t, x2; title="particle 2", xlabel="t", ylabel="",
              label=false)
    end
    return p
end

function plot_position1(sim; invert=false)
    pgfplotsx()
    p = plot()

    nt = length(sim.sol.t)
    x1 = zeros(nt)
    for i in 1:nt
        x1[i] = sim.sol.u[i].x[2][2]
    end
    if invert
        plot!(p, x1, sim.sol.t; title="particle 1", xlabel="x", ylabel="t", label=false)
    else
        plot!(p, sim.sol.t, x1; title="particle 1", xlabel="t", ylabel="x", label=false)
    end
    savefig(p, "particle_location.pdf")
    return p
end

function plot_position_compare_sf(sim1, sim2, invert=false)
    pgfplotsx()

    @assert all(sim1.sol.t .== sim2.sol.t)
    p = plot(; title=L"m=%$(sim1.params.m),
q=%$(sim1.params.q), v_{0}=%$(sim1.params.v0), x_{0}=%$(sim1.params.x0)")

    nt1 = length(sim1.sol.t)
    nt2 = length(sim2.sol.t)
    x1 = zeros(nt1)
    x2 = zeros(nt2)

    for i in 1:nt1
        x1[i] = sim1.sol.u[i].x[2][2]
    end
    for i in 1:nt2
        x2[i] = sim2.sol.u[i].x[2][2]
    end

    if sim2.params.sf
        xtmp = copy(x1)
        x1 .= x2
        x2 .= xtmp
    end
    if invert
        plot!(p, x1, sim1.sol.t; xlabel="x", ylabel="t", label="self force")
        plot!(p, x2, sim2.sol.t; label="no self force")
    else
        plot!(p, sim1.sol.t, x1; xlabel="t", ylabel="x", label="self force")
        plot!(p, sim2.sol.t, x2; label="no self force")
    end
    savefig(p, "comparison_position.pdf")
    return p
end

function plot_position_resolutions(sim1, sim2, sim4; invert=false)
    pgfplotsx()
    p = plot()

    nt = length(sim1.sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    x4 = zeros(nt)
    for i in 1:nt
        x1[i] = sim1.sol.u[i].x[2][2]
        x2[i] = sim2.sol.u[2i - 1].x[2][2]
        x4[i] = sim4.sol.u[4i - 3].x[2][2]
    end
    if sim1.params.sf
        label_low = "low + sf"
        label_mid = "mid + sf"
        label_high = "high + sf"
    else
        label_low = "low"
        label_mid = "mid"
        label_high = "high"
    end
    if invert
        plot!(p, x1, sim1.sol.t; label=label_low)
        plot!(p, x2, sim1.sol.t; label=label_mid)
        plot!(p, x4, sim1.sol.t; label=label_high)
    else
        plot!(p, sim1.sol.t, x1; label=label_low)
        plot!(p, sim1.sol.t, x2; label=label_mid)
        plot!(p, sim1.sol.t, x4; label=label_high)
    end
    savefig(p, "particle_location_resolutions.pdf")
    return p
end

function plot_position_errors(sim1, sim2, sim4, scale=2; invert=false)
    pgfplotsx()
    p = plot(; title="Errors")

    nt = length(sim1.sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    x4 = zeros(nt)
    for i in 1:nt
        x1[i] = sim1.sol.u[i].x[2][2]
        x2[i] = sim2.sol.u[2i - 1].x[2][2]
        x4[i] = sim4.sol.u[4i - 3].x[2][2]
    end
    if sim1.params.sf
        label_low = "sf: low - mid"
        label_mid = "sf: mid - high"
    else
        label_low = "low - mid"
        label_mid = "mid - high"
    end
    if invert
        xlabel!(p, "x")
        ylabel!(p, "t")
        plot!(p, x1 .- x2, sim1.sol.t; label=label_low)
        plot!(p, scale .* (x2 .- x4), sim1.sol.t; label=label_mid)
        # scatter!(p, x4, sim1.sol.t; label=label_high)
    else
        ylabel!(p, "x")
        xlabel!(p, "t")
        plot!(p, sim1.sol.t, x1 .- x2; label=label_low)
        plot!(p, sim1.sol.t, scale .* (x2 .- x4); label=label_mid)
    end
    savefig(p, "particle_location_errors.pdf")
    return p
end

function plot_position_conv_factor(sim1, sim2, sim4, scale=2; invert=false)
    pgfplotsx()
    p = plot(; title=L"Convergence Factor $= 2^p$")

    nt = length(sim1.sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    x4 = zeros(nt)
    for i in 1:nt
        x1[i] = sim1.sol.u[i].x[2][2]
        x2[i] = sim2.sol.u[2i - 1].x[2][2]
        x4[i] = sim4.sol.u[4i - 3].x[2][2]
    end
    err1 = x1 - x2
    err2 = x2 - x4
    co = err1 ./ err2

    xlabel!(p, "t")
    plot!(p, sim1.sol.t, co; label=false)
    savefig(p, "particle_location_convergence_factor.pdf")
    return p
end

function plot_oscillation(sim)
    nt = length(sim.sol.t)
    x1s = zeros(nt)
    x2s = zeros(nt)
    x1 = sim.params.x1
    vmax = sim.params.vmax
    A = sim.params.A
    direction1 = sim.params.direction1
    direction2 = sim.params.direction2
    for i in 1:nt
        x1s[i], _, _ = ParticleMotion.oscillator2(sim.sol.t[i], x1, vmax, A,
                                                  direction1)
    end
    p = plot(sim.sol.t, x1s)
    return p
end

end # end of module
