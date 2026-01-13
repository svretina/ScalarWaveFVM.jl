module RetardedPotentialPlots

using ..PlottingUtils
using ..RetardedPotentials
using ..ScalarField
using CairoMakie

function plot_retarded_Φ(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x = (-L / 2):dx:(L / 2)
    q = sim.params.q
    A = sim.params.A
    w = sim.params.ω
    z = A * sin(w * 0)
    v = A * w * cos(w * 0)
    y = Φ_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Φs.(x, z, v, q)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_{ret}",
              title=L"t=%$(round(t, digits=3))")
    lines!(ax, x, y)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Π(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x_range = (-L / 2, L / 2)
    x = (-L / 2):dx:(L / 2)
    # y = Π_ret.(x, Ref(sim.params), Ref(t)) .- 0Π_ret.(x, Ref(sim.params), Ref(0.0))
    # y[div(length(y), 2) + 1] = 0.0
    # y = RetardedPotentials.Π_ret_autodiff.(x, Ref(sim.params), Ref(t)) .-
    #     0RetardedPotentials.Π_ret_autodiff.(x, Ref(sim.params), Ref(0.0))
    q = sim.params.q
    A = sim.params.A
    w = sim.params.ω
    z = A * sin(w * t)
    v = A * w * cos(w * t)
    y = Π_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Πs.(x, z, v, q)#0Ψ_ret2.(x, Ref(sim.params), Ref(0.0))
    # y .+= 0.005
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
              ylabel=L"\Pi_{\textrm{ret}}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    lines!(ax, x, y)
    # lines!(ax, x, y2)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Π2(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x_range = (-L / 2, L / 2)
    x = (-L / 2):dx:(L / 2)
    # y = Π_ret.(x, Ref(sim.params), Ref(t)) .- 0Π_ret.(x, Ref(sim.params), Ref(0.0))

    y = RetardedPotentials.Π_ret2.(x, Ref(sim.params), Ref(t)) .-
        0RetardedPotentials.Π_ret2.(x, Ref(sim.params), Ref(0.0))

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
    #
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Ψ(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x_range = (-L / 2, L / 2)
    x = (-L / 2):dx:(L / 2)

    # y = Ψ_ret.(x, Ref(sim.params), Ref(t)) # 0Ψ_ret.(x, Ref(sim.params), Ref(0.0))

    q = sim.params.q
    A = sim.params.A
    w = sim.params.ω
    z = A * sin(w * 0)
    v = A * w * cos(w * 0)
    y = Ψ_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Ψs.(x, z, v, q)#0Ψ_ret2.(x, Ref(sim.params), Ref(0.0))

    y[(div(length(y), 2) + 1)] = 0.0
    ymin = minimum(y)
    ymax = maximum(y)
    # y .-= ymax * sign.(x)
    # y[(div(length(y), 2) + 1)] = 0.0
    # ymin = minimum(y)
    # ymax = maximum(y)

    @show ymin
    @show ymax
    padding = round(0.1(ymax - ymin); digits=5)
    @show padding
    y_limits = round.((ymin - padding, ymax + padding), digits=5)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_{ret}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x, y)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Ψ2(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    dx = sim.params.dx
    L = sim.params.L
    x_range = (-L / 2, L / 2)
    x = (-L / 2):dx:(L / 2)

    q = sim.params.q
    A = sim.params.A
    w = sim.params.ω
    z = A * sin(w * t)
    v = A * w * cos(w * t)
    y = Ψ_ret2.(x, Ref(sim.params), Ref(t)) .- ScalarField.Ψs.(x, z, v, q)#0Ψ_ret2.(x, Ref(sim.params), Ref(0.0))
    @show y
    ymin = minimum(y)
    ymax = maximum(y[(div(length(y), 2) - 1):(div(length(y), 2) + 2)])
    # y .-= ymax * sign.(x)
    # y[(div(length(y), 2) + 1)] = 0.0
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
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Π_at() end

end # end of module
