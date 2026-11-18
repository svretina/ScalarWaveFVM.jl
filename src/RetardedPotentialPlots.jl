module RetardedPotentialPlots

using ..PlottingUtils
using ..RetardedPotentials
using ..ScalarField
using ..ParticleMotion
using CairoMakie

function plot_retarded_Φ(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    h = sim.params.dx
    L = sim.params.L
    x = (-L / 2):h:(L / 2)
    y = Φ_ret.(x, Ref(sim.params), Ref(t)) .- 0Φ_ret.(x, Ref(sim.params), Ref(0.0))
    y2 = Φ_ret_exact.(x, Ref(sim.params), Ref(t)) .-
         0Φ_ret_exact.(x, Ref(sim.params), Ref(0.0))

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_{ret}",
              title=L"t=%$(round(t, digits=3))")
    lines!(ax, x, y; label=L"\Phi_{ret}(t) - \Phi_{s}(t)")
    lines!(ax, x, y2; linestyle=:dash,
           label=L"\Phi_{ret}(t) - \Phi_{ret}^{exact}(t)")
    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, framevisible=false, colgap=5, rowgap=5, patchsize=(15, 10))
    #xlims!(ax, -1A, 10A)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Π(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    h = sim.params.dx
    L = sim.params.L
    x = (-t):h:t

    z, v, _ = ParticleMotion.oscillator(t, sim.params.z0, sim.params.A, sim.params.ω)
    y = Π_ret.(x, Ref(sim.params), Ref(t))
    y2 = Π_ret_exact.(x, Ref(sim.params), Ref(t))
    y3 = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(sim.params.q))
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_{\textrm{ret}}",
              title=L"t=%$(round(t, digits=3))"
              # ,limits=(x[1], x[end], y_limits[1], y_limits[2])
              )
    lines!(ax, x, y; label=L"\Pi_{ret}(t)")
    lines!(ax, x, y2; linestyle=:dash,
           label=L"\Pi_{ret}^{exact}(t)")
    lines!(ax, x, y3; linestyle=:dashdot, label=L"\Pi_{s}(z(t), t)")
    #lines!(ax, x, y4; linestyle=:dot, label=L"\Pi_{s}^{approx}(z(t), t)")
    #     lines!(ax, x, y2_0; linestyle=:dot, label=L"\Pi20)")
    #     lines!(ax, x, y1_0; linestyle=:dot, label=L"\Pi10")

    xlo, xhi = -0.1, 0.1
    mask = (x .>= xlo) .& (x .<= xhi)
    ywin_min = min(minimum(y[mask]), minimum(y2[mask]), minimum(y3[mask]))
    ywin_max = max(maximum(y[mask]), maximum(y2[mask]), maximum(y3[mask]))
    ypad = 0.1 * (ywin_max - ywin_min + eps(Float64))

    xlims!(ax, xlo, xhi)
    ylims!(ax, ywin_min - ypad, ywin_max + ypad)
    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, framevisible=false, colgap=5, rowgap=5, patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.18))
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Π_errors(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    h = sim.params.dx
    L = sim.params.L
    x = (-L / 2):h:(L / 2)

    y = Π_ret.(x, Ref(sim.params), Ref(t))
    y2 = Π_ret_exact.(x, Ref(sim.params), Ref(t))
    y[(div(length(y), 2) + 1)] = 0.0
    y2[(div(length(y), 2) + 1)] = 0.0
    ymin = minimum(y .- y2)
    ymax = maximum(y .- y2)
    padding = round(0.1(ymax - ymin); digits=3)
    y_limits = round.((ymin - padding, ymax + padding), digits=3)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_{\textrm{ret}}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x[1], x[end], y_limits[1], y_limits[2]))

    lines!(ax, x, y .- y2; label=L"\Pi_{s}(t_{ret}) - \Pi_{ret}^{exact}(t_{ret})")
    # lines!(ax, x, y2; linestyle=:dash,
    #        label=L"\Pi_{ret}^{exact}(t) - \Pi_{ret}^{exact}(0)")

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, framevisible=false, colgap=5, rowgap=5, patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.18))
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Π_at_particle(sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    ts = sim.sol.t
    zp = RetardedPotentials.z.(ts, Ref(sim.params))

    y = Π_ret.(zp, Ref(sim.params), ts)
    y2 = Π_ret_exact.(zp, Ref(sim.params), ts)

    ymin = min(minimum(y), minimum(y2))
    ymax = max(maximum(y), maximum(y2))
    padding = ymin == ymax ? 0.1 : 0.1 * (ymax - ymin)
    y_limits = (ymin - padding, ymax + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Pi_{\textrm{ret}}(z(t), t)",
              title=L"x = z(t)",
              limits=(ts[1], ts[end], y_limits[1], y_limits[2]))

    lines!(ax, ts, y; label=L"\Pi_{s}(z(t), t_{ret}) - \Pi_{s}(z(0), 0)")
    lines!(ax, ts, y2; linestyle=:dash,
           label=L"\Pi_{ret}^{\textrm{exact}}(z(t), t) - \Pi_{ret}^{\textrm{exact}}(z(t), 0)")

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, framevisible=false, colgap=5, rowgap=5, patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.18))
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Ψ_at_particle(sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    ts = sim.sol.t
    zp = RetardedPotentials.z.(ts, Ref(sim.params))

    y = Ψ_ret.(zp, Ref(sim.params), ts) .- Ψ_ret.(zp, Ref(sim.params), zero(first(ts)))
    y2 = Ψ_ret_exact.(zp, Ref(sim.params), ts) .-
         Ψ_ret_exact.(zp, Ref(sim.params), zero(first(ts)))

    ymin = min(minimum(y), minimum(y2))
    ymax = max(maximum(y), maximum(y2))
    padding = ymin == ymax ? 0.1 : 0.1 * (ymax - ymin)
    y_limits = (ymin - padding, ymax + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Psi_{\textrm{ret}}(z(t), t)",
              title=L"x = z(t)",
              limits=(ts[1], ts[end], y_limits[1], y_limits[2]))

    lines!(ax, ts, y; label=L"\Psi_{s}(z(t_{ret}), t_{ret}) - \Psi_{s}(z(0), 0)")
    lines!(ax, ts, y2; linestyle=:dash,
           label=L"\Psi_{ret}^{\textrm{exact}}(z(t), t) - \Psi_{ret}^{\textrm{exact}}(z(t), 0)")

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, framevisible=false, colgap=5, rowgap=5, patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.18))
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_self_force(sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    ts = sim.sol.t
    zp = RetardedPotentials.z.(ts, Ref(sim.params))
    vp = RetardedPotentials.v.(ts, Ref(sim.params))
    Pi = Π_ret.(zp, Ref(sim.params), ts) .- 0Π_ret.(zp, Ref(sim.params), zero(first(ts)))
    Psi = Ψ_ret.(zp, Ref(sim.params), ts) .- 0Ψ_ret.(zp, Ref(sim.params), zero(first(ts)))
    y = vp .* Pi .+ Psi

    Pi2 = Π_ret_exact.(zp, Ref(sim.params), ts) .-
          0 * Π_ret_exact.(zp, Ref(sim.params), zero(first(ts)))
    Psi2 = Ψ_ret_exact.(zp, Ref(sim.params), ts) .-
           0 * Ψ_ret_exact.(zp, Ref(sim.params), zero(first(ts)))
    y2 = vp .* Pi2 .+ Psi2

    ymin = min(minimum(y), minimum(y2))
    ymax = max(maximum(y), maximum(y2))
    padding = ymin == ymax ? 0.1 : 0.1 * (ymax - ymin)
    y_limits = (ymin - padding, ymax + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{Self Force}",
              title=L"\textrm{Self Force}",
              limits=(ts[1], ts[end], y_limits[1], y_limits[2]))

    lines!(ax, ts, y; label=L"v\Pi_{ret}(z(t), t) + \Psi_{ret}(z(t), t)")
    lines!(ax, ts, y2; linestyle=:dash,
           label=L"v\Pi_{ret}^{\textrm{exact}}(z(t), t) + \Psi_{ret}^{\textrm{exact}}(z(t), t)")

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, framevisible=false, colgap=5, rowgap=5, patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.18))
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_retarded_Ψ(t, sim)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    h = sim.params.dx
    L = sim.params.L
    x = (-L / 2):h:(L / 2)

    y = Ψ_ret.(x, Ref(sim.params), Ref(t)) .- 0Ψ_ret.(x, Ref(sim.params), Ref(0.0))
    y2 = Ψ_ret_exact.(x, Ref(sim.params), Ref(t)) .-
         0Ψ_ret_exact.(x, Ref(sim.params), Ref(0.0))
    y[(div(length(y), 2) + 1)] = 0.0
    y2[(div(length(y), 2) + 1)] = 0.0

    ymin = min(minimum(y), minimum(y2))
    ymax = max(maximum(y), maximum(y2))

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
              limits=(x[1], x[end], y_limits[1], y_limits[2]))

    lines!(ax, x, y; label=L"\Psi_{ret}(t) - \Psi_{ret}(0)")
    lines!(ax, x, y2; linestyle=:dash,
           label=L"\Psi_{ret}^{\textrm{exact}}(t) - \Psi_{ret}^{\textrm{exact}}(0)")
    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, framevisible=false, colgap=5, rowgap=5, patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.18))
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

end # end of module
