module ParticleHistoryPlots

using CairoMakie
using LaTeXStrings
using ..PlottingUtils: mytheme_aps
using ..Diagnostics

export plot_particle_history, plot_remainder_energy_balance

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

"""Plot the exact discrete remainder-energy budget and its closure residual."""
function plot_remainder_energy_balance(sim; savepath=nothing)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    balance = Diagnostics.remainder_energy_balance(sim)

    fig = Figure(; size=(253, 520), figure_padding=(1, 5, 5, 1))
    ax1 = Axis(fig[1, 1]; ylabel=L"\Delta E_{r,h}",
               title=L"\textrm{Discrete Energy Balance}")
    lines!(ax1, balance.t, balance.energy_change; label=L"\Delta E_{r,h}")
    lines!(ax1, balance.t, balance.predicted_change; linestyle=:dash,
           label=L"\int(P_S-D_h+R_h)\,dt")
    axislegend(ax1; position=:rb, labelsize=8, rowgap=0.1)

    ax2 = Axis(fig[2, 1]; ylabel=L"\textrm{cumulative}")
    lines!(ax2, balance.t, balance.integrated_source; label=L"\int P_S\,dt")
    lines!(ax2, balance.t, -balance.integrated_dissipation;
           label=L"-\int D_h\,dt")
    lines!(ax2, balance.t, balance.integrated_reconstruction;
           label=L"\int R_h\,dt")
    axislegend(ax2; position=:rb, labelsize=8, rowgap=0.1)

    ax3 = Axis(fig[3, 1]; xlabel=L"t", ylabel=L"\mathcal{B}_h(t)",
               title=L"\textrm{Closure Residual}")
    lines!(ax3, balance.t, balance.residual)
    hidexdecorations!(ax1; grid=false)
    hidexdecorations!(ax2; grid=false)
    linkxaxes!(ax1, ax2, ax3)
    rowgap!(fig.layout, 4)
    colsize!(fig.layout, 1, Relative(0.9))

    savepath === nothing || save(savepath, fig)
    return fig
end

end
