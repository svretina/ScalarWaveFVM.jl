module RetardedPotentialPlots

using ..PlottingUtils
using ..RetardedPotentials
using CairoMakie

function plot_retarded_Φ(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x = (-L / 2):dx:(L / 2)
    y = Φ_ret.(x, Ref(sim.params), Ref(t))

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_{ret}",
              title=L"t=%$(round(t, digits=3))")
    lines!(ax, x, y)
    return fig
end

function plot_retarded_Π(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x_range = (-L / 2, L / 2)
    x = (-L / 2):dx:(L / 2)
    y = Π_ret.(x, Ref(sim.params), Ref(t))

    y2 = RetardedPotentials.Π_ret_autodiff.(x, Ref(sim.params), Ref(t))

    ymin = minimum(y)
    ymax = maximum(y)
    @show ymin
    @show ymax
    padding = round(0.1(ymax - ymin); digits=3)
    @show padding
    y_limits = round.((ymin - padding, ymax + padding), digits=3)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_{ret}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    lines!(ax, x, y)
    # lines!(ax, x, y2)
    return fig
end

function plot_retarded_Ψ(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x_range = (-L / 2, L / 2)
    x = (-L / 2):dx:(L / 2)
    y = Ψ_ret.(x, Ref(sim.params), Ref(t))

    ymin = minimum(y)
    ymax = maximum(y)
    @show ymin
    @show ymax
    padding = round(0.1(ymax - ymin); digits=3)
    @show padding
    y_limits = round.((ymin - padding, ymax + padding), digits=3)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_{ret}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    lines!(ax, x, y)
    return fig
end

function plot_retarded_Π_at() end

end # end of module
