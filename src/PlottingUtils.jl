module PlottingUtils

using LaTeXStrings
using ..ParticleMotion
using ..Run
using LinearAlgebra
using MakieCore
using Makie
using CairoMakie

const PROJ_PATH = pkgdir(@__MODULE__)
const DEFAULT_FIG_PATH = joinpath(PROJ_PATH, "figures")
const FIG_PATH_REF = Ref(DEFAULT_FIG_PATH)
const PAPER_PATH = "/home/svretina/PhD/mypapers/ScalarSelfForcePaper"
const PAPER_FIG_PATH_BASE = joinpath(PAPER_PATH, "figures")
const EXPERIMENT_TYPE_REF = Ref("")
const RUN_NAME_REF = Ref("")

export mytheme_aps, mymean, simple_name
export get_fig_path, set_experiment_path!, reset_fig_path!
export get_paper_fig_path
export apply_scientific_notation_if_needed!

"""
    apply_scientific_notation_if_needed!(ax)

Checks if the y-axis limits have 3 or more decimal places. If so, applies scientific
notation by scaling the y-ticks to simple integers and adding a scale indicator label
at the top of the y-axis.
"""
function apply_scientific_notation_if_needed!(ax)
    # Most reliable way to get data bounds before the scene is fully rendered
    # is to iterate over the plots and collect their data limits.
    ymin = Inf
    ymax = -Inf

    plots_to_check = ax.scene.plots
    # Some axes have a background or other non-data plots, 
    # but Makie.data_limits usually handles them.
    for p in plots_to_check
        # We only care about data plots (Lines, Scatter, etc.)
        # Skip labels/text if any were added
        if p isa Makie.Text || p isa Makie.Annotations
            continue
        end

        bb = Makie.data_limits(p)
        if !isnothing(bb) && isfinite(bb.origin[2]) && isfinite(bb.widths[2])
            ymin = min(ymin, bb.origin[2])
            ymax = max(ymax, bb.origin[2] + bb.widths[2])
        end
    end

    if !isfinite(ymin) || !isfinite(ymax) || ymin >= ymax
        return
    end

    range_val = ymax - ymin

    # Trigger if range implies multiple decimal places (e.g., range < 0.5)
    if range_val < 0.5 && range_val > 0
        exponent = floor(Int, log10(range_val))
        scale = 10.0^exponent

        scaled_min = ymin / scale
        scaled_max = ymax / scale

        # We want integer ticks in the scaled space
        s_min = floor(Int, scaled_min)
        s_max = ceil(Int, scaled_max)

        tick_vals_scaled = collect(s_min:s_max)

        if length(tick_vals_scaled) > 8
            step = ceil(Int, length(tick_vals_scaled) / 5)
            tick_vals_scaled = tick_vals_scaled[1:step:end]
        elseif length(tick_vals_scaled) < 4
            # Fallback to a few nice points if range is too narrow for many integers
            tick_vals_scaled = unique(round.(Int, range(scaled_min, scaled_max; length=5)))
        end

        tick_values = tick_vals_scaled .* scale
        mask = (tick_values .>= ymin - 1e-10) .& (tick_values .<= ymax + 1e-10)

        if any(mask)
            ax.yticks = (tick_values[mask],
                         [string(round(Int, s)) for s in tick_vals_scaled[mask]])

            # Use the approach suggested by the USER for scientific notation scale label
            scale_text = L"\times 10^{%$exponent}"

            # Placing it in the Top() side of the layout position where the axis is located
            try
                # In Makie 0.19+, blocks have layoutmain.gridcontent
                # we can find the parent GridLayout and the span
                gc = ax.layoutmain.gridcontent[]
                if !isnothing(gc)
                    Label(gc.parent[gc.span.rows, gc.span.cols, Top()], scale_text;
                          halign=:left, tellwidth=false, tellheight=false,
                          padding=(0, 0, 2, 0), fontsize=10, font=:bold)
                else
                    # Fallback to text! if gridcontent is not found
                    text!(ax, 0, 1; text=scale_text, align=(:left, :bottom),
                          space=:relative, offset=(0, 2), fontsize=10, font=:bold,
                          overdraw=true)
                end
            catch
                # Defensive fallback
                text!(ax, 0, 1; text=scale_text, align=(:left, :bottom),
                      space=:relative, offset=(0, 2), fontsize=10, font=:bold,
                      overdraw=true)
            end
        end
    end

    return nothing
end

function format_val(x)
    # Round to 4 digits to remove floating point noise like 1.250000002 -> 1.25
    # string(round(x, digits=4))
    # Or simply:
    val = round(x; digits=4)
    if isinteger(val)
        return string(Int(val))
    else
        return string(val)
    end
end

function get_fig_path()
    return FIG_PATH_REF[]
end

function get_paper_fig_path()
    path = joinpath(PAPER_FIG_PATH_BASE, EXPERIMENT_TYPE_REF[], RUN_NAME_REF[])
    if !isdir(path)
        mkpath(path)
    end
    return path
end

function reset_fig_path!()
    FIG_PATH_REF[] = DEFAULT_FIG_PATH
    EXPERIMENT_TYPE_REF[] = ""
    RUN_NAME_REF[] = ""
    return nothing
end

function set_experiment_path!(experiment_type, q, vmax, Nosc, dx, cfl)
    EXPERIMENT_TYPE_REF[] = experiment_type
    run_name = "run_q=$(format_val(q))_vmax=$(format_val(vmax))_Nosc=$(format_val(Nosc))_dx=$(format_val(dx))_cfl=$(format_val(cfl))"
    RUN_NAME_REF[] = run_name
    new_path = joinpath(DEFAULT_FIG_PATH, experiment_type, run_name)
    if !isdir(new_path)
        mkpath(new_path)
    end
    FIG_PATH_REF[] = new_path
    return nothing
end

function simple_name(base_name, time=nothing)
    if time === nothing
        return "$(base_name).pdf"
    else
        return "$(base_name)_t=$(format_val(time)).pdf"
    end
end

function mymean(arr)
    arr2 = zeros(eltype(arr), div(length(arr), 2))
    for i in eachindex(arr2)
        j = 2i - 1
        arr2[i] = 0.5(arr[j] + arr[j + 1])
    end
    return arr2
end

function mytheme_aps()
    return Theme(
                 # Axis attributes
                 ;
                 Axis=Attributes(; spinewidth=1.1,
                                 xgridvisible=true,
                                 xlabelpadding=-0,
                                 xlabelsize=12,
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
                                 ygridvisible=true,
                                 ylabelpadding=2,
                                 ylabelsize=12,
                                 yminortickalign=1,
                                 yminorticks=IntervalsBetween(5, true),
                                 yminorticksize=3,
                                 yminorticksvisible=true,
                                 yminortickwidth=0.75,
                                 ytickalign=1,
                                 yticklabelsize=10,
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
                                  cycle=Cycle([[:color] => :color],
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

function opposite_charges(sim)
    q1 = sim.params.q1
    q2 = sim.params.q2
    if q1 == q2
        return false
    else
        return true
    end
end

end # end of module
