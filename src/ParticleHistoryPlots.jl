module ParticleHistoryPlots

using CairoMakie
using LaTeXStrings
using ..PlottingUtils: mytheme_aps
using ..Diagnostics

export plot_particle_history, plot_particle_trajectory_resolutions,
       plot_remainder_energy_balance, plot_self_force_validation

function _turning_points(t, position)
    indices = [i for i in 2:(length(position) - 1)
               if abs(position[i]) >= abs(position[i - 1]) &&
                  abs(position[i]) > abs(position[i + 1])]
    return t[indices], abs.(position[indices])
end

function _position_history(data)
    if hasproperty(data, :sol)
        particle = getindex.(getproperty.(data.sol.u, :x), 2)
        return data.sol.t, getindex.(particle, 2)
    end
    return data.t, data.position
end

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

"""Overlay particle trajectories and damping envelopes from several resolutions."""
function plot_particle_trajectory_resolutions(histories, labels;
                                              reference=nothing,
                                              reference_label=L"\textrm{reduced ODE}",
                                              savepath=nothing)
    length(histories) == length(labels) ||
        throw(ArgumentError("histories and labels must have the same length"))
    isempty(histories) && throw(ArgumentError("at least one history is required"))

    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    linestyles = (:dot, :dash, :solid)

    fig = Figure(; size=(253, 390), figure_padding=(1, 5, 5, 1))
    axz = Axis(fig[1, 1]; ylabel=L"z(t)",
               title=L"\textrm{Trajectory Resolution Comparison}")
    axa = Axis(fig[2, 1]; xlabel=L"t", ylabel=L"|z|_{\mathrm{turn}}",
               title=L"\textrm{Damping Envelope}")

    if reference !== nothing
        peak_time, peak_amplitude = _turning_points(reference.t, reference.position)
        lines!(axz, reference.t, reference.position;
               color=:black, linewidth=1, label=reference_label)
        scatterlines!(axa, peak_time, peak_amplitude;
                      color=:black, linewidth=1, markersize=3,
                      label=reference_label)
    end

    for (i, (history, label)) in enumerate(zip(histories, labels))
        t, position = _position_history(history)
        peak_time, peak_amplitude = _turning_points(t, position)
        linestyle = linestyles[min(i, length(linestyles))]
        line = lines!(axz, t, position; label, linestyle)
        scatterlines!(axa, peak_time, peak_amplitude;
                      color=line.color, linestyle, markersize=3, label)
    end

    axislegend(axz; position=:rb, labelsize=8, rowgap=0.1)
    hidexdecorations!(axz; grid=false)
    linkxaxes!(axz, axa)
    rowgap!(fig.layout, 4)
    colsize!(fig.layout, 1, Relative(0.9))

    savepath === nothing || save(savepath, fig)
    return fig
end

"""Plot the three convergence checks for the complete self-force model."""
function plot_self_force_validation(field_h, field_errors, h, force_error,
                                    particle_errors; savepath=nothing)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    fig = Figure(; size=(253, 520), figure_padding=(1, 5, 5, 1))
    ax_field = Axis(fig[1, 1]; xscale=log10, yscale=log10,
                    ylabel=L"L^1\ \textrm{error}",
                    title=L"\textrm{Retarded Field}")
    ax_force = Axis(fig[2, 1]; xscale=log10, yscale=log10,
                    ylabel=L"\max_t|\Psi_r+v\Pi_r|_p",
                    title=L"\textrm{Spatial Self-Force}")
    ax_particle = Axis(fig[3, 1]; xscale=log10, yscale=log10,
                       xlabel=L"\Delta x", ylabel=L"\max_t|e|",
                       title=L"\textrm{Reduced Particle Dynamics}")

    for (errors, label) in zip(field_errors,
                               (L"\Pi_r", L"\Psi_r", L"\Phi_r"))
        scatterlines!(ax_field, field_h, errors; label, markersize=5)
    end
    scatterlines!(ax_force, h, force_error; markersize=5)
    for (errors, label) in zip(particle_errors, (L"z", L"v", L"m"))
        scatterlines!(ax_particle, h, errors; label, markersize=5)
    end

    axislegend(ax_field; position=:lt, labelsize=8, rowgap=0.1)
    axislegend(ax_particle; position=:lt, labelsize=8, rowgap=0.1)
    hidexdecorations!(ax_field; grid=false)
    hidexdecorations!(ax_force; grid=false)
    linkxaxes!(ax_field, ax_force, ax_particle)
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
