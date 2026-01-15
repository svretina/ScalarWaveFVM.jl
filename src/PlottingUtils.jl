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
