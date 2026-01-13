module PaperPlots

using ..PlottingUtils

using ..Run
using ..RetardedPotentialPlots
using ..ForcedPlots

using ..ParticleMotion
using ..RetardedPotentials
using ..ScalarField

using LinearAlgebra
using CairoMakie

@inline name(pre, sim, post=nothing) = name_forced(pre, sim, post)

function plot_Π(i, sim; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim.sol.t[i]

    dx = sim.params.dx
    L = sim.params.L
    q = sim.params.q
    A = sim.params.A
    ω = sim.params.ω
    z0 = sim.params.z0

    x_range = (-L / 2, L / 2)
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    x = sim.params.x[mask]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Pis .+ sim.sol.u[i][mask, 1]
    else
        # xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        # Pis0 = ScalarField.Πs.(x, Ref(xp0), Ref(vp0), Ref(q))
        y_sim = sim.sol.u[i][mask, 1] #.- 0Pis0
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret.(x, Ref(sim.params), Ref(t))
    else
        y_lw = Π_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Πs.(x, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x, y_sim; label="Simulation")
    lines!(ax, x, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_Π(i, sim1, sim2; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim1.sol.t[i]

    dx1 = sim1.params.dx
    dx2 = sim2.params.dx
    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    x_range = (-L / 2, L / 2)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])

    x1 = sim1.params.x[mask1]
    x2 = sim2.params.x[mask2]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x1, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Pis .+ sim1.sol.u[i][mask1, 1]
    else
        # xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        # Pis0 = ScalarField.Πs.(x, Ref(xp0), Ref(vp0), Ref(q))
        y_sim1 = sim1.sol.u[i][mask1, 1] #.- 0Pis0
        y_sim2 = sim2.sol.u[2i - 1][mask2, 1]
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret.(x1, Ref(sim1.params), Ref(t))
    else
        y_lw = Π_ret.(x1, Ref(sim1.params), Ref(t)) .- ScalarField.Πs.(x1, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim1)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim1)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_r",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x1, y_sim1; label="dx=$dx1")
    lines!(ax, x2, y_sim2; label="dx=$dx2")
    lines!(ax, x1, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))

    savename = name("forced_motion/pifield_retarded_resolutions", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_Ψ(i, sim; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim.sol.t[i]

    dx = sim.params.dx
    L = sim.params.L
    q = sim.params.q
    A = sim.params.A
    ω = sim.params.ω
    z0 = sim.params.z0

    x_range = (-L / 2, L / 2)
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    x = sim.params.x[mask]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Psis = ScalarField.Ψs.(x, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Psis .+ sim.sol.u[i][mask, 2]
    else
        xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        Psis0 = ScalarField.Ψs.(x, Ref(xp0), Ref(vp0), Ref(q))
        y_sim = Psis .+ sim.sol.u[i][mask, 2] .- Psis0
    end

    ## Lienard-Wiechert

    if total
        y_lw = Ψ_ret.(x, Ref(sim.params), Ref(t))
    else
        y_lw = Ψ_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Ψs.(x, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x, y_sim; label="Simulation")
    lines!(ax, x, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

### Convergence

function plot_Π_convergence(i, sim1, sim2; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    dx1 = sim1.params.dx
    dx2 = sim2.params.dx

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    x_range = (-L / 2, L / 2)

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x = sim1.params.x[mask1]

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    @assert t1 == t2
    t = t1

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))
    if total
        y_sim1 = Pis .+ sim1.sol.u[i][mask1, 1]
        y_sim2 = Pis .+ mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    else
        # xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        # Pis0 = ScalarField.Πs.(x, Ref(xp0), Ref(vp0), Ref(q))
        # y_sim1 = Pis .+ sim1.sol.u[i][mask1, 1] .- Pis0
        # y_sim2 = Pis .+ mymean(sim2.sol.u[2i - 1][:, 1])[mask1] .- Pis0

        y_sim1 = sim1.sol.u[i][mask1, 1]
        y_sim2 = mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret.(x, Ref(sim1.params), Ref(t))
    else
        # xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        y_lw = Π_ret.(x, Ref(sim1.params), Ref(t)) .- ScalarField.Πs.(x, z, v, q)
    end

    diff1 = y_sim1 .- y_lw
    diff2 = y_sim2 .- y_lw

    convfac = diff1 ./ diff2
    @show norm(convfac)
    @show norm(diff1) / norm(diff2)
    @show log2(norm(diff1, 1) / norm(diff2, 1))
    # ymin = minimum(convfac)
    # ymax = maximum(convfac)

    # padding = round(0.1(ymax - ymin); digits=3)

    # y_limits = round.((ymin - padding, ymax + padding), digits=3)
    #y_limits = (0.0, 5.0)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              # yscale=log10,
              # ylabel=L"",
              title=L"\textrm{Differences \Pi for t=%$(round(t, digits=3)) }")
    # ,
    # limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    # lines!(ax, x, convfac)
    lines!(ax, x, diff1)
    lines!(ax, x, diff2)

    lines!(ax, x, y_lw; linestyle=:dash)
    lines!(ax, x, y_sim1)
    lines!(ax, x, y_sim2)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

function plot_Π_convergence(i, Nosc, dx0, cfl, vmax, q; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q)

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    nt = length(sim1.sol.t)
    convfac = zeros(nt - 1)

    x_range = (-L / 2, L / 2)

    hs = [dx0, dx0 / 2, dx0 / 4, dx0 / 8, dx0 / 16]
    errors = zeros(length(hs))

    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q)
    sim8 = Run.forced_motion(Nosc, dx0 / 8, cfl, vmax, q)
    sim16 = Run.forced_motion(Nosc, dx0 / 16, cfl, vmax, q)

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    t8 = sim8.sol.t[8i - 7]
    t16 = sim16.sol.t[16i - 15]

    @assert t1 == t2
    @assert t2 == t4
    @assert t4 == t8
    @assert t8 == t16

    # ### Method 1 (averaging)
    # mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # x1 = sim1.params.x[mask1]

    # z, v, _ = ParticleMotion.oscillator(t1, z0, A, ω)
    # y_lw = Π_ret.(x1, Ref(sim1.params), Ref(t1)) .- ScalarField.Πs.(x1, z, v, q)

    # y_sim1 = sim1.sol.u[i][mask1, 1]
    # y_sim2 = mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    # y_sim4 = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))[mask1]
    # y_sim8 = mymean(mymean(mymean(sim8.sol.u[8i - 7][:, 1])))[mask1]
    # y_sim16 = mymean(mymean(mymean(mymean(sim16.sol.u[16i - 15][:, 1]))))[mask1]

    # mode = 1
    # errors[1] = norm(abs.(y_lw .- y_sim1), mode)
    # errors[2] = norm(abs.(y_lw .- y_sim2), mode)
    # errors[3] = norm(abs.(y_lw .- y_sim4), mode)
    # errors[4] = norm(abs.(y_lw .- y_sim8), mode)
    # errors[5] = norm(abs.(y_lw .- y_sim16), mode)

    ### Method 2
    z, v, _ = ParticleMotion.oscillator(t1, z0, A, ω)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1 = sim1.params.x[mask1]

    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    x2 = sim2.params.x[mask2]

    mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    x4 = sim4.params.x[mask4]

    mask8 = (sim8.params.x .>= x_range[1]) .& (sim8.params.x .<= x_range[2])
    x8 = sim8.params.x[mask8]

    mask16 = (sim16.params.x .>= x_range[1]) .& (sim16.params.x .<= x_range[2])
    x16 = sim16.params.x[mask16]

    y_lw1 = Π_ret.(x1, Ref(sim1.params), Ref(t1)) .- ScalarField.Πs.(x1, z, v, q)
    y_lw2 = Π_ret.(x2, Ref(sim2.params), Ref(t1)) .- ScalarField.Πs.(x2, z, v, q)
    y_lw4 = Π_ret.(x4, Ref(sim4.params), Ref(t1)) .- ScalarField.Πs.(x4, z, v, q)
    y_lw8 = Π_ret.(x8, Ref(sim8.params), Ref(t1)) .- ScalarField.Πs.(x8, z, v, q)
    y_lw16 = Π_ret.(x16, Ref(sim16.params), Ref(t1)) .- ScalarField.Πs.(x16, z, v, q)

    y_sim1 = sim1.sol.u[i][mask1, 1]
    y_sim2 = sim2.sol.u[2i - 1][mask2, 1]
    y_sim4 = sim4.sol.u[4i - 3][mask4, 1]
    y_sim8 = sim8.sol.u[8i - 7][mask8, 1]
    y_sim16 = sim16.sol.u[16i - 15][mask16, 1]

    mode = Inf
    errors[1] = norm(abs.(y_lw1 .- y_sim1), mode)
    errors[2] = norm(abs.(y_lw2 .- y_sim2)[begin:2:end], mode)
    errors[3] = norm(abs.(y_lw4 .- y_sim4)[begin:4:end], mode)
    errors[4] = norm(abs.(y_lw8 .- y_sim8)[begin:8:end], mode)
    errors[5] = norm(abs.(y_lw16 .- y_sim16)[begin:16:end], mode)

    @show hs
    @show errors
    @show t1

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"h",
              # ylabel=L"\Pi",
              yscale=log10,
              xscale=log10,
              xreversed=true,
              title=L"\Vert |\Pi_{\textrm{ret}}-\Pi_\textrm{h}| \Vert_{\infty}")

    lines!(ax, hs, errors)
    scatter!(ax, hs, errors)

    xref = hs[1:3]
    yref1 = hs[1:3] .^ 1
    shift1 = (errors[1] / yref1[1]) * 1.3
    yref1 .*= shift1

    lines!(ax, xref, yref1; linestyle=:dash, label=L"h^1")
    # axislegend(ax; position=:rt)

    # leg = Legend(fig[1, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)
    leg = Legend(fig[1, 1], ax; # Note: Same grid position as the Axis!
                 halign=:right,
                 valign=:top,
                 tellwidth=false,      # Important: Don't shrink the axis to fit the legend
                 tellheight=false,     # Important: Don't shrink the axis to fit the legend
                 margin=(10, 10, 10, 10), # Optional: Add padding from the axis lines
                 framevisible=false,
                 labelsize=12)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))

    savename = name("forced_motion/pifield_retarded_convergence", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_Ψ_convergence(i, sim1, sim2; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim.sol.t[i]

    dx = sim.params.dx
    L = sim.params.L
    q = sim.params.q
    A = sim.params.A
    ω = sim.params.ω
    z0 = sim.params.z0

    x_range = (-L / 2, L / 2)
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    x = sim.params.x[mask]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Psis = ScalarField.Ψs.(x, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Psis .+ sim.sol.u[i][mask, 2]
    else
        xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        Psis0 = ScalarField.Ψs.(x, Ref(xp0), Ref(vp0), Ref(q))
        y_sim = Psis .+ sim.sol.u[i][mask, 2] .- Psis0
    end

    ## Lienard-Wiechert

    if total
        y_lw = Ψ_ret.(x, Ref(sim.params), Ref(t))
    else
        y_lw = Ψ_ret.(x, Ref(sim.params), Ref(t)) .- ScalarField.Ψs.(x, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim)
    ymax = maximum([ymax1, ymax2])

    padding = round(0.1(ymax - ymin); digits=3)

    y_limits = round.((ymin - padding, ymax + padding), digits=3)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x, y_sim; label="Simulation")
    lines!(ax, x, y_lw; label="Retarded", linestyle=:dash)

    vlines!(ax, [z]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    return fig
end

end # end of module
