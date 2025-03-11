module PlottingUtils

using LaTeXStrings
using ..ParticleMotion
using ..Run
using LinearAlgebra
using MakieCore
using Makie
using CairoMakie
CairoMakie.activate!()

const PROJ_PATH = pkgdir(@__MODULE__)
const FIG_PATH = joinpath(PROJ_PATH, "figures")
const PAPER_PATH = "/home/svretina/PhD/mypapers/ScalarSelfForcePaper"
const PAPER_FIG_PATH = joinpath(PAPER_PATH, "figs")

export mytheme_aps, mymean, name_potential, name_interacting, name_forced
export FIG_PATH, PAPER_FIG_PATH

function name_potential(pre, sim, post=nothing)
    q = sim.params.q
    v = sim.params.v0
    dx = sim.params.dx
    λ = sim.params.λ
    cfl = sim.params.cfl
    if post === nothing
        tmp = join([pre,
                    "v=$(v)",
                    "q=$(q)",
                    "dx=$(dx)",
                    "cfl=$(cfl)",
                    "pdf"],
                   "_", ".")
    else
        tmp = join([pre,
                    "v=$(v)_q=$(q)_dx=$(dx)_cfl=$(cfl)",
                    post, "pdf"],
                   "_", ".")
    end
    return tmp
end

function name_interacting(pre, sim, post=nothing)
    q1 = sim.params.q1
    q2 = sim.params.q2

    v1 = sim.params.v1
    v2 = sim.params.v2

    dx = sim.params.dx
    cfl = sim.params.cfl
    if post === nothing
        tmp = join([pre,
                    "q1=$(q1)",
                    "q2=$(q2)",
                    "v1=$(v1)",
                    "v2=$(v2)",
                    "dx=$(dx)",
                    "cfl=$(cfl)",
                    "pdf"],
                   "_", ".")
    else
        tmp = join([pre,
                    "q1=$(q1)",
                    "q2=$(q2)",
                    "v1=$(v1)",
                    "v2=$(v2)",
                    "dx=$(dx)",
                    "cfl=$(cfl)",
                    post, "pdf"],
                   "_", ".")
    end
    return tmp
end

function name_forced(pre, sim, post=nothing)
    q = sim.params.q

    vmax = sim.params.vmax
    Nosc = sim.params.Nosc
    dx = sim.params.dx
    cfl = sim.params.cfl
    if post === nothing
        tmp = join([pre,
                    "q=$(q)",
                    "vmax=$(vmax)",
                    "Nosc=$(Nosc)",
                    "dx=$(dx)",
                    "cfl=$(cfl)",
                    "pdf"],
                   "_", ".")
    else
        tmp = join([pre,
                    "q=$(q)",
                    "vmax=$(vmax)",
                    "Nosc=$(Nosc)",
                    "dx=$(dx)",
                    "cfl=$(cfl)",
                    post, "pdf"],
                   "_", ".")
    end
    return tmp
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

end # end of module
