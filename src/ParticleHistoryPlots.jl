module ParticleHistoryPlots

using CairoMakie
using LaTeXStrings
using ..PlottingUtils: mytheme_aps

export plot_particle_history

"""Plot the evolved dynamical mass, position, and velocity of one particle."""
function plot_particle_history(sim; savepath=nothing)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    t = sim.sol.t
    particle = getindex.(getproperty.(sim.sol.u, :x), 2)
    mass = getindex.(particle, 1)
    position = getindex.(particle, 2)
    velocity = getindex.(particle, 3)

    fig = Figure(; size=(253, 520), figure_padding=(1, 5, 5, 1))
    axm = Axis(fig[1, 1]; ylabel=L"m(t)", title=L"\textrm{Particle Mass}")
    axz = Axis(fig[2, 1]; ylabel=L"z(t)", title=L"\textrm{Particle Position}")
    axv = Axis(fig[3, 1]; xlabel=L"t", ylabel=L"v(t)",
               title=L"\textrm{Particle Velocity}")
    lines!(axm, t, mass)
    lines!(axz, t, position)
    lines!(axv, t, velocity)
    hidexdecorations!(axm; grid=false)
    hidexdecorations!(axz; grid=false)
    linkxaxes!(axm, axz, axv)
    rowgap!(fig.layout, 4)
    colsize!(fig.layout, 1, Relative(0.9))

    savepath === nothing || save(savepath, fig)
    return fig
end

end
