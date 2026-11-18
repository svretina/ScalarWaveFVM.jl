module ParticlePotentialPlots

using ..ParticleMotion
using ..PlottingUtils
using ..ScalarField
using ..InitialData
using ..Run
using ..ODE
using ..Equations
using CairoMakie
using LinearAlgebra
using DataInterpolations
using Roots
using OrdinaryDiffEq

function trapezoidal(y, h)
    n = length(y) - 1
    integral = zero(eltype(y))
    @inbounds for i in 1:n
        integral += y[i]
    end
    @inbounds integral += 0.5 * (y[1] + y[n + 1])  # Adjust endpoints
    integral *= h
    return integral
end

function deriv(y, h)
    n = length(y)
    dy = zeros(eltype(y), n)
    for i in 2:(n - 1)
        dy[i] = (y[i + 1] - y[i - 1]) / (2h)
    end

    dy[1] = (y[1] - y[2]) / h
    dy[n] = (y[n - 1] - y[n]) / h
    return dy
end

function cumulative_trapezoidal(t, y)
    n = length(t)
    out = zeros(promote_type(eltype(t), eltype(y)), n)
    @inbounds for i in 2:n
        dt = t[i] - t[i - 1]
        out[i] = out[i - 1] + 0.5 * dt * (y[i - 1] + y[i])
    end
    return out
end

function cumulative_trapezoidal(t, dt, y)
    n = length(t)
    out = zeros(promote_type(eltype(t), eltype(y)), n)
    @inbounds for i in 2:n
        out[i] = out[i - 1] + 0.5 * dt * (y[i - 1] + y[i])
    end
    return out
end

function get_velocity(sim)
    if hasproperty(sim.sol.prob.p, :equation)
        return sim.sol.prob.p.equation.velocity
    end
    return 1.0
end

function non_radiative_reference_solution(sim)
    return Run.solve_non_radiative_relativistic_particle(sim.params)
end

function particle_centered_interpolant(ξ, y)
    if length(ξ) >= 5
        return AkimaInterpolation(y, ξ; extrapolation=ExtrapolationType.Extension)
    end
    return LinearInterpolation(y, ξ; extrapolation=ExtrapolationType.Extension)
end

function relativistic_self_force_at_particle(sim, i)
    q = sim.sol.prob.p.q
    v = sim.sol.u[i].x[2][3]
    return q * (1 - v * v) * regular_force_combo_at_particle(sim, i)
end

function regular_force_combo_at_particle(sim, i)
    params = sim.sol.prob.p

    z = sim.sol.u[i].x[2][2]
    v = sim.sol.u[i].x[2][3]
    Q_field = sim.sol.u[i].x[1]
    Πrz, Ψrz = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)

    return Ψrz + v * Πrz
end

function plot_regular_force_combo_convergence(sim1, sim2, sim4, sim8, sim16; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    g1 = zeros(nt)
    g2 = zeros(nt)
    g4 = zeros(nt)
    g8 = zeros(nt)
    g16 = zeros(nt)

    for i in 1:nt
        g1[i] = regular_force_combo_at_particle(sim1, i)
        g2[i] = regular_force_combo_at_particle(sim2, 2i - 1)
        g4[i] = regular_force_combo_at_particle(sim4, 4i - 3)
        g8[i] = regular_force_combo_at_particle(sim8, 8i - 7)
        g16[i] = regular_force_combo_at_particle(sim16, 16i - 15)
    end

    diff12 = max.(abs.(g1 .- g2), eps(Float64))
    diff24 = max.(abs.(g2 .- g4), eps(Float64))
    diff48 = max.(abs.(g4 .- g8), eps(Float64))
    diff816 = max.(abs.(g8 .- g16), eps(Float64))

    palette = mytheme.palette
    color_low = palette.color[][1]
    color_mid = palette.color[][2]
    color_high = palette.color[][3]
    color_finer = palette.color[][4]
    color_finest = palette.color[][5]

    fig = Figure(; size=(253, 330), figure_padding=(1, 5, 5, 1))

    ax1 = Axis(fig[1, 1];
               xlabel=L"t",
               ylabel=L"\Psi_r(z) + v \Pi_r(z)",
               title=L"\textrm{Particle-sampled regular field}")
    lines!(ax1, sim1.sol.t, zeros(nt); color=:black, linestyle=:dot, linewidth=1.0,
           label=L"0")
    lines!(ax1, sim1.sol.t, g1; color=color_low, linewidth=1.5,
           label=L"h")
    lines!(ax1, sim1.sol.t, g2; color=color_mid, linestyle=:dash, linewidth=1.4,
           label=L"h/2")
    lines!(ax1, sim1.sol.t, g4; color=color_high, linestyle=:dot, linewidth=1.4,
           label=L"h/4")
    lines!(ax1, sim1.sol.t, g8; color=color_finer, linestyle=:dashdot, linewidth=1.3,
           label=L"h/8")
    lines!(ax1, sim1.sol.t, g16; color=color_finest, linestyle=:solid, linewidth=1.2,
           label=L"h/16")
    axislegend(ax1; position=:rb, labelsize=9, rowgap=0.1)
    apply_scientific_notation_if_needed!(ax1)

    ax2 = Axis(fig[2, 1];
               xlabel=L"t",
               ylabel=L"\Delta",
               yscale=log10,
               title=L"\textrm{Successive resolution differences}")
    lines!(ax2, sim1.sol.t, diff12; color=color_low, linewidth=1.5,
           label=L"|Q_h - Q_{h/2}|")
    lines!(ax2, sim1.sol.t, diff24; color=color_mid, linestyle=:dash, linewidth=1.4,
           label=L"|Q_{h/2} - Q_{h/4}|")
    lines!(ax2, sim1.sol.t, diff48; color=color_high, linestyle=:dot, linewidth=1.4,
           label=L"|Q_{h/4} - Q_{h/8}|")
    lines!(ax2, sim1.sol.t, diff816; color=color_finer, linestyle=:dashdot, linewidth=1.3,
           label=L"|Q_{h/8} - Q_{h/16}|")
    axislegend(ax2; position=:rb, labelsize=9, rowgap=0.1)

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Auto())
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 4)

    savename = simple_name("psi_vpi_particle_convergence")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_phifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    t1 = sim1.sol.t[i]
    t = t1
    @show t
    q = sim1.params.q
    xp = sim1.sol.u[i].x[2][2]
    vp = sim1.sol.u[i].x[2][3]
    # Create figure and axis

    tf = sim1.sol.t[end]
    x_range = (-tf, tf)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = mask # (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 3] .+ ScalarField.Φs.(x_filtered, xp, vp, q)

    # Compute y-limits with padding
    # y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    # y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi",
              title=L"t=%$(round(t, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\Phi")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = simple_name("total_phifield", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_phifield_atx(ix, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    c = get_velocity(sim1)
    # if iszero(t)
    #     x_range = (-sim1.params.L / 2, sim1.params.L / 2)
    # else
    #     x_range = (-c * t, c * t)
    # end
    # mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    # mask1 = mask # (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[ix]  # ~320 elements if symmetric around 0
    y1_filtered = zeros(length(sim1.sol.t))
    for i in 1:length(sim1.sol.t)
        q = sim1.params.q
        xp = sim1.sol.u[i].x[2][2]
        vp = sim1.sol.u[i].x[2][3]
        y1_filtered[i] = sim1.sol.u[i].x[1][ix, 3] .+
                         ScalarField.Φs.(x1_filtered, xp, vp, q)
    end
    ts = sim1.sol.t
    # Compute y-limits with padding
    # y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    # y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)
    @show max(y1_filtered...)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Phi",
              title=L"x=%$(round(x1_filtered, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, ts, y1_filtered; label=L"\Pi")
    # vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # savename = simple_name("total_phifield", t)
    # save(joinpath(get_fig_path(), savename), fig)
    # if paper
    #     save(joinpath(get_paper_fig_path(), savename), fig)
    # end
    return fig
end

function plot_total_pifield_atx(ix, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    nt = length(sim1.sol.t)
    x1_filtered = sim1.params.x[ix]  # ~320 elements if symmetric around 0
    y1_filtered = zeros(nt)
    q = sim1.params.q
    ps = zeros(nt)
    for i in 1:nt
        xp = sim1.sol.u[i].x[2][2]
        vp = sim1.sol.u[i].x[2][3]
        ps[i] = ScalarField.Πs.(x1_filtered, xp, vp, q)
        y1_filtered[i] = sim1.sol.u[i].x[1][ix, 1]
    end
    ts = sim1.sol.t

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Pi",
              title=L"x=%$(round(x1_filtered, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, ts, y1_filtered; label=L"\Pi_{r}")
    lines!(ax, ts, ps; label=L"\Pi_s")
    lines!(ax, ts, y1_filtered .+ ps; linestyle=:dash, label=L"\Pi")
    # vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # savename = simple_name("total_phifield", t)
    # save(joinpath(get_fig_path(), savename), fig)
    # if paper
    #     save(joinpath(get_paper_fig_path(), savename), fig)
    # end
    return fig
end

function plot_total_psifield_atx(ix, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    nt = length(sim1.sol.t)
    x1_filtered = sim1.params.x[ix]  # ~320 elements if symmetric around 0
    y1_filtered = zeros(nt)
    q = sim1.params.q
    ps = zeros(nt)
    for i in 1:nt
        xp = sim1.sol.u[i].x[2][2]
        vp = sim1.sol.u[i].x[2][3]
        ps[i] = ScalarField.Ψs.(x1_filtered, xp, vp, q)
        y1_filtered[i] = sim1.sol.u[i].x[1][ix, 2]
    end
    ts = sim1.sol.t

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Psi",
              title=L"x=%$(round(x1_filtered, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, ts, y1_filtered; label=L"\Psi_{r}")
    lines!(ax, ts, ps; label=L"\Psi_s")
    lines!(ax, ts, y1_filtered .+ ps; linestyle=:dash, label=L"\Psi")
    # vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # savename = simple_name("total_phifield", t)
    # save(joinpath(get_fig_path(), savename), fig)
    # if paper
    #     save(joinpath(get_paper_fig_path(), savename), fig)
    # end
    return fig
end

function plot_total_phifield_diff(ti, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if ti == -1
        ti = lastindex(sim1.sol.t)
    end

    t1 = sim1.sol.t[ti]
    ti0 = 1
    @show t1
    q = sim1.params.q
    n = length(sim1.params.x)
    idx_zero = div(n - 1, 2)
    ## Particle Position and Velocity at t=ti
    ## at a later time
    xp1 = sim1.sol.u[ti].x[2][2]
    vp1 = sim1.sol.u[ti].x[2][3]

    ## Particle Position and Velocity at t=ti0
    ## at the start of the simulation
    xp10 = sim1.sol.u[ti0].x[2][2]
    vp10 = sim1.sol.u[ti0].x[2][3]

    ## Filter our grid affected from boundary information
    ## propagating inwards
    L = sim1.params.L
    Nosc = sim1.params.Nosc
    k = sim1.params.k
    T = 2π / sqrt(abs(k))
    tf = Nosc * T
    L0 = 2tf
    x_range = (10, L0 / 4)#(-L0 / 4, L0 / 4)
    @show x_range
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]

    ## Fields at start of simulation
    Φ10_filtered = sim1.sol.u[ti0].x[1][mask1, 3] .+
                   ScalarField.Φs.(x1_filtered, xp10, vp10, q)
    ## Fields at t=ti of simulation
    Φ1_filtered = sim1.sol.u[ti].x[1][mask1, 3] .+ ScalarField.Φs.(x1_filtered, xp1, vp1, q)

    ## Π fields
    Π10_filtered = sim1.sol.u[ti0].x[1][mask1, 1] .+
                   ScalarField.Πs.(x1_filtered, xp10, vp10, q)
    ## Fields at t=ti of simulation
    Π1_filtered = sim1.sol.u[ti].x[1][mask1, 1] .+ ScalarField.Πs.(x1_filtered, xp1, vp1, q)

    ## Ψ fields
    Ψ10_filtered = sim1.sol.u[ti0].x[1][mask1, 2] .+
                   ScalarField.Ψs.(x1_filtered, xp10, vp10, q)

    ## Fields at t=ti of simulation
    Ψ1_filtered = sim1.sol.u[ti].x[1][mask1, 2] .+ ScalarField.Ψs.(x1_filtered, xp1, vp1, q)

    ## Figure
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi(ti=%$ti0)-\Phi(ti=%$ti)",
              title=L"t=%$(round(t1, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    # lines!(ax, x1_filtered, y10_filtered; label=L"init")
    # lines!(ax, x1_filtered, y1_filtered; linestyle=:dash, label=L"end")
    @show Φ10_filtered[end] .- Φ1_filtered[end]
    @show Π10_filtered[end] .- Π1_filtered[end]
    @show Ψ10_filtered[end] .- Ψ1_filtered[end]
    lines!(ax, x1_filtered, Φ10_filtered .- Φ1_filtered; label=L"Δ\Phi")
    lines!(ax, x1_filtered, Π10_filtered .- Π1_filtered; linestyle=:dash, label=L"Δ\Pi")
    lines!(ax, x1_filtered, Ψ10_filtered .- Ψ1_filtered; linestyle=:dash,
           label=L"\Delta\Psi")

    vlines!(ax, [xp1]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("particle_potential/total_phifield_diff", sim1, "i=$(ti)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_phifield_diff_res(ti, xi, sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if ti == -1
        ti = lastindex(sim1.sol.t)
    end

    t1 = sim1.sol.t[ti]
    t2 = sim2.sol.t[2ti - 1]
    t4 = sim4.sol.t[4ti - 3]
    ti0 = 1

    @assert t1 == t2
    @assert t2 == t4

    @assert sim1.params.q == sim2.params.q
    @assert sim2.params.q == sim4.params.q

    q = sim1.params.q
    n = length(sim1.params.x)
    idx_zero = div(n - 1, 2)
    ## Particle Position and Velocity at t=ti
    ## at a later time
    xp1 = sim1.sol.u[ti].x[2][2]
    vp1 = sim1.sol.u[ti].x[2][3]

    xp2 = sim2.sol.u[2ti - 1].x[2][2]
    vp2 = sim2.sol.u[2ti - 1].x[2][3]

    xp4 = sim4.sol.u[4ti - 3].x[2][2]
    vp4 = sim4.sol.u[4ti - 3].x[2][3]

    ## Particle Position and Velocity at t=ti0
    ## at the start of the simulation
    xp10 = sim1.sol.u[ti0].x[2][2]
    vp10 = sim1.sol.u[ti0].x[2][3]

    xp20 = sim2.sol.u[2ti0 - 1].x[2][2]
    vp20 = sim2.sol.u[2ti0 - 1].x[2][3]

    xp40 = sim4.sol.u[4ti0 - 3].x[2][2]
    vp40 = sim4.sol.u[4ti0 - 3].x[2][3]

    ## Location where to evaluate the fields
    x1 = sim1.params.x[idx_zero + 1 + xi]
    x2 = sim2.params.x[2idx_zero + 1 + 2xi]
    x4 = sim4.params.x[4idx_zero + 1 + 4xi]

    @assert x1 == x2
    @assert x2 == x4

    ## Filter our grid affected from boundary information
    ## propagating inwards
    F = 2
    L = sim1.params.L
    x_range = (-L / F, L / F)
    @assert x1 > x_range[1] && x1 < x_range[2]

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])

    x1_filtered = sim1.params.x[mask1]
    x2_filtered = sim2.params.x[mask2]
    x4_filtered = sim4.params.x[mask4]

    ## Fields at start of simulation
    Φ10_filtered = sim1.sol.u[ti0].x[1][mask1, 3] .+
                   ScalarField.Φs.(x1_filtered, xp10, vp10, q)
    Φ20_filtered = sim2.sol.u[ti0].x[1][mask2, 3] .+
                   ScalarField.Φs.(x2_filtered, xp20, vp20, q)
    Φ40_filtered = sim4.sol.u[ti0].x[1][mask4, 3] .+
                   ScalarField.Φs.(x4_filtered, xp40, vp40, q)

    ## Fields at t=ti of simulation
    Φ1_filtered = sim1.sol.u[ti].x[1][mask1, 3] .+ ScalarField.Φs.(x1_filtered, xp1, vp1, q)
    Φ2_filtered = sim2.sol.u[2ti - 1].x[1][mask2, 3] .+
                  ScalarField.Φs.(x2_filtered, xp2, vp2, q)
    Φ4_filtered = sim4.sol.u[4ti - 3].x[1][mask4, 3] .+
                  ScalarField.Φs.(x4_filtered, xp4, vp4, q)

    ## Π fields
    Π10_filtered = sim1.sol.u[ti0].x[1][mask1, 1] .+
                   ScalarField.Πs.(x1_filtered, xp10, vp10, q)
    Π20_filtered = sim2.sol.u[ti0].x[1][mask2, 1] .+
                   ScalarField.Πs.(x2_filtered, xp20, vp20, q)
    Π40_filtered = sim4.sol.u[ti0].x[1][mask4, 1] .+
                   ScalarField.Πs.(x4_filtered, xp40, vp40, q)

    ## Fields at t=ti of simulation
    Π1_filtered = sim1.sol.u[ti].x[1][mask1, 1] .+ ScalarField.Πs.(x1_filtered, xp1, vp1, q)
    Π2_filtered = sim2.sol.u[2ti - 1].x[1][mask2, 1] .+
                  ScalarField.Πs.(x2_filtered, xp2, vp2, q)
    Π4_filtered = sim4.sol.u[4ti - 3].x[1][mask4, 1] .+
                  ScalarField.Πs.(x4_filtered, xp4, vp4, q)

    ## Ψ fields
    Ψ10_filtered = sim1.sol.u[ti0].x[1][mask1, 2] .+
                   ScalarField.Ψs.(x1_filtered, xp10, vp10, q)
    Ψ20_filtered = sim2.sol.u[ti0].x[1][mask2, 2] .+
                   ScalarField.Ψs.(x2_filtered, xp20, vp20, q)
    Ψ40_filtered = sim4.sol.u[ti0].x[1][mask4, 2] .+
                   ScalarField.Ψs.(x4_filtered, xp40, vp40, q)

    ## Fields at t=ti of simulation
    Ψ1_filtered = sim1.sol.u[ti].x[1][mask1, 2] .+ ScalarField.Ψs.(x1_filtered, xp1, vp1, q)
    Ψ2_filtered = sim2.sol.u[2ti - 1].x[1][mask2, 2] .+
                  ScalarField.Ψs.(x2_filtered, xp2, vp2, q)
    Ψ4_filtered = sim4.sol.u[4ti - 3].x[1][mask4, 2] .+
                  ScalarField.Ψs.(x4_filtered, xp4, vp4, q)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi(ti=2)-\Phi(ti=%$ti)",
              title=L"t=%$(round(t1, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    # lines!(ax, x1_filtered, y10_filtered; label=L"init")
    # lines!(ax, x1_filtered, y1_filtered; linestyle=:dash, label=L"end")

    lines!(ax, x1_filtered, Φ10_filtered .- Φ1_filtered; label=L"\Phi")
    lines!(ax, x2_filtered, Π10_filtered .- Π1_filtered; linestyle=:dash, label=L"\Pi")
    lines!(ax, x4_filtered, Ψ10_filtered .- Ψ1_filtered; linestyle=:dash, label=L"\Psi")

    vlines!(ax, [xp1]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = name("particle_potential/total_phifield_diff_res", sim1, "i=$(ti)")
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_regular_phifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    ti = sim1.sol.t[i]
    tf = sim1.sol.t[end]

    # Create figure and axis
    x_range = (-tf, tf)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = mask# (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 3]  # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp = sim1.sol.u[i].x[2][2]

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_\mathrm{r}",
              title=L"t=%$(round(ti, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\Phi_r")
    apply_scientific_notation_if_needed!(ax)
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # Save the figure
    savename = simple_name("phifield", i)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_regular_pifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    t = sim1.sol.t[i]
    tf = sim1.sol.t[end]

    # Create figure and axis

    x_range = (-tf, tf)
    # Filter sim1 data (640 points)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]  # ~320 elements if symmetric around 0
    y_filtered = sim1.sol.u[i].x[1][mask, 1]  # Matches x1_filtered length

    xp = sim1.sol.u[i].x[2][2]

    y_min = minimum(y_filtered)
    y_max = maximum(y_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=5))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x_filtered, y_filtered; label=L"\Pi_r")

    apply_scientific_notation_if_needed!(ax)
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1], ax; orientation=:horizontal,
    #              tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # Save the figure
    savename = simple_name("pifield", i)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_pifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    t = sim1.sol.t[i]
    k = sim1.params.k
    tf = sim1.sol.t[end]

    # Create figure and axis

    x_range = (-tf, tf)
    # Filter sim1 data (640 points)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]  # ~320 elements if symmetric around 0
    y_filtered = sim1.sol.u[i].x[1][mask, 1]  # Matches x1_filtered length

    q = sim1.params.q
    h = sim1.params.dx
    xp = sim1.sol.u[i].x[2][2]
    vp = sim1.sol.u[i].x[2][3]

    y_min = minimum(y_filtered)
    y_max = maximum(y_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)
    ps = ScalarField.Πs.(x_filtered, xp, vp, q)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x_filtered, y_filtered; label=L"\Pi_r")
    lines!(ax, x_filtered, y_filtered .+ ps; label=L"\Pi")

    # lines!(ax, x_filtered, ScalarField.Πs.(x_filtered, xp, vp, q); label=L"\Pi_s")
    # lines!(ax, x_filtered, y_filtered .+ ScalarField.Πs.(x_filtered, xp, vp, q);
    #        label=L"\Pi")

    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # Save the figure
    savename = simple_name("total_pifield", i)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_resolutions(i, sim1, sim2, sim4; paper=false)
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1

    # Create figure and axis

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = mask# (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 1])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp = sim1.sol.u[i].x[2][2]

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    # Save the figure
    savename = simple_name("pifield_res", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_regular_psifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    ti = sim1.sol.t[i]
    tf = sim1.sol.t[end]

    # Create figure and axis

    x_range = (-tf, tf)

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 2]  # Matches x1_filtered length

    # Compute y-limits with padding
    # y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    # y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)

    q = sim1.params.q
    xp = sim1.sol.u[i].x[2][2]
    vp = sim1.sol.u[i].x[2][3]

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(ti, digits=3))")
    # ,
    # limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x1_filtered, y1_filtered; label=L"\Psi_r")
    # lines!(ax, x1_filtered, ScalarField.Ψs.(x1_filtered, xp, vp, q); label=L"\Psi_s")
    # lines!(ax, x1_filtered, y1_filtered .+ ScalarField.Ψs.(x1_filtered, xp, vp, q);
    #        linestyle=:dash, label=L"\Psi")

    apply_scientific_notation_if_needed!(ax)
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # Save the figure
    savename = simple_name("psifield", i)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_psifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    t = sim1.sol.t[i]
    k = sim1.params.k
    tf = sim1.sol.t[end]

    # Create figure and axis

    x_range = (-tf, tf)
    # Filter sim1 data (640 points)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]  # ~320 elements if symmetric around 0
    y_filtered = sim1.sol.u[i].x[1][mask, 1]  # Matches x1_filtered length
    q = sim1.params.q
    xp = sim1.sol.u[i].x[2][2]
    vp = sim1.sol.u[i].x[2][3]
    ys = ScalarField.Ψs.(x_filtered, xp, vp, q)

    y_min = minimum(ys)
    y_max = maximum(ys)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x_filtered, y_filtered; label=L"\Psi_r")
    lines!(ax, x_filtered, y_filtered .+ ys; label=L"\Psi")

    lines!(ax, x_filtered, ys; label=L"\Psi_s")
    # lines!(ax, x_filtered, ScalarField.Πs.(x_filtered, xp, vp, q); label=L"\Pi_s")
    # lines!(ax, x_filtered, y_filtered .+ ScalarField.Πs.(x_filtered, xp, vp, q);
    #        label=L"\Pi")

    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # Save the figure
    savename = simple_name("total_psifield", i)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function particle_centered_field_slice(i, sim, component)
    t = sim.sol.t[i]
    c = get_velocity(sim)
    L = sim.params.L

    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end

    x = sim.params.x
    mask = (x .>= x_range[1]) .& (x .<= x_range[2])
    xp = sim.sol.u[i].x[2][2]
    ξ = x[mask] .- xp
    y = sim.sol.u[i].x[1][mask, component]
    return ξ, y, xp, t
end

function plot_field_resolutions_particle_centered(i, sim1, sim2, sim4;
                                                  component, ylabel, basename,
                                                  paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1

    ξ1, y1, _, _ = particle_centered_field_slice(i, sim1, component)
    ξ2, y2, _, _ = particle_centered_field_slice(2i - 1, sim2, component)
    ξ4, y4, _, _ = particle_centered_field_slice(4i - 3, sim4, component)

    y_min = minimum([minimum(y1), minimum(y2), minimum(y4)])
    y_max = maximum([maximum(y1), maximum(y2), maximum(y4)])
    padding = 0.1 * (y_max - y_min + eps(eltype(y1)))
    ξ_min = minimum([minimum(ξ1), minimum(ξ2), minimum(ξ4)])
    ξ_max = maximum([maximum(ξ1), maximum(ξ2), maximum(ξ4)])

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"\xi = x - z_p(t)",
              ylabel=ylabel,
              title=L"t=%$(round(t, digits=3))",
              limits=(ξ_min, ξ_max, y_min - padding, y_max + padding))

    lines!(ax, ξ1, y1; label=L"\textrm{low}")
    lines!(ax, ξ2, y2; label=L"\textrm{mid}")
    lines!(ax, ξ4, y4; label=L"\textrm{high}")
    vlines!(ax, [0.0]; color=:darkgray, linestyle=:dash, linewidth=0.8)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.1))
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name(basename, t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_resolutions_particle_centered(i, sim1, sim2, sim4; paper=false)
    return plot_field_resolutions_particle_centered(i, sim1, sim2, sim4;
                                                    component=1,
                                                    ylabel=L"\Pi_\mathrm{r}",
                                                    basename="pifield_res_particle_centered",
                                                    paper=paper)
end

function plot_psifield_resolutions_particle_centered(i, sim1, sim2, sim4; paper=false)
    return plot_field_resolutions_particle_centered(i, sim1, sim2, sim4;
                                                    component=2,
                                                    ylabel=L"\Psi_\mathrm{r}",
                                                    basename="psifield_res_particle_centered",
                                                    paper=paper)
end

function particle_centered_self_convergence_factor(sim1, sim2, sim4, component)
    nt = length(sim1.sol.t)
    convfac = zeros(nt)

    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4

        ξ1, y1, _, _ = particle_centered_field_slice(i, sim1, component)
        ξ2, y2, _, _ = particle_centered_field_slice(2i - 1, sim2, component)
        ξ4, y4, _, _ = particle_centered_field_slice(4i - 3, sim4, component)
        if min(length(ξ1), length(ξ2), length(ξ4)) < 2
            convfac[i] = NaN
            continue
        end

        ξref = ξ1
        y1ref = y1

        interp2 = particle_centered_interpolant(ξ2, y2)
        interp4 = particle_centered_interpolant(ξ4, y4)
        y2ref = interp2.(ξref)
        y4ref = interp4.(ξref)

        diff1 = y1ref .- y2ref
        diff2 = y2ref .- y4ref
        convfac[i] = norm(diff1, 1) / max(norm(diff2, 1), eps(eltype(convfac)))
    end
    return convfac
end

function plot_field_self_convergence_order_particle_centered(sim1, sim2, sim4;
                                                             component, title, basename,
                                                             paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    convfac = particle_centered_self_convergence_factor(sim1, sim2, sim4, component)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              title=title,
              limits=(0, sim1.sol.t[end], 0, 4))
    mask = isfinite.(convfac[2:end]) .& (convfac[2:end] .> 0)
    lines!(ax, sim1.sol.t[2:end][mask], log2.(convfac[2:end][mask]); linewidth=1)

    rowsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name(basename)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_self_convergence_order_particle_centered(sim1, sim2, sim4;
                                                               paper=false)
    return plot_field_self_convergence_order_particle_centered(sim1, sim2, sim4;
                                                               component=1,
                                                               title=L"\textrm{Self Convergence Order of }\Pi_\mathrm{r}\textrm{ in } \xi",
                                                               basename="pifield_convord_particle_centered_potential",
                                                               paper=paper)
end

function plot_psifield_self_convergence_order_particle_centered(sim1, sim2, sim4;
                                                                paper=false)
    return plot_field_self_convergence_order_particle_centered(sim1, sim2, sim4;
                                                               component=2,
                                                               title=L"\textrm{Self Convergence Order of }\Psi_\mathrm{r}\textrm{ in } \xi",
                                                               basename="psifield_convord_particle_centered_potential",
                                                               paper=paper)
end

function plot_psifield_resolutions(i, sim1, sim2, sim4; paper=false)
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @show t1, t2, t4
    @assert t1 == t2
    @assert t2 == t4
    t = t1

    # Create figure and axis

    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = mask #(sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 2]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 2])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 2]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp = sim1.sol.u[i].x[2][2]

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # Save the figure
    savename = simple_name("psifield_res", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_self_convergence_order(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    nt = length(sim1.sol.t)
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    convfac = zeros(nt)
    c = get_velocity(sim1)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

        if iszero(t)
            x_range = (-L / 2, L / 2)
            continue
        else
            x_range = (-c * t, c * t)
        end
        # x_range = (-c * t, c * t)
        mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])

        y1_raw = sim1.sol.u[i].x[1][:, 1]  # Matches x1_filtered length
        y1_filtered = y1_raw[mask1]

        y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 1])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Self Convergence Order of \Pi_\mathrm{r} field}",
              limits=(0, sim1.sol.t[end], 0, 4))
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = simple_name("pifield_convord_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_self_convergence_order(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    nt = length(sim1.sol.t)
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    convfac = zeros(nt)
    c = get_velocity(sim1)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

        if iszero(t)
            x_range = (-L / 2, L / 2)
            continue
        else
            x_range = (-c * t, c * t)
        end
        mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])

        y1_raw = sim1.sol.u[i].x[1][:, 2]  # Matches x1_filtered length
        y1_filtered = y1_raw[mask1]

        y2_raw = mymean(sim2.sol.u[2i - 1].x[1][:, 2])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 2]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Self Convergence Order of \Psi_\mathrm{r} field}",
              limits=(0, sim1.sol.t[end], 0, 4))
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = simple_name("psifield_convord_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Particle
function plot_particle_position(sim1; invert=true, paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xnr = zeros(nt)
    non_radiative_sol = non_radiative_reference_solution(sim1)

    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][2]
        xnr[i] = non_radiative_sol.u[i][2]
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}",
              xticks=LinearTicks(5),
              xgridvisible=true,
              ygridvisible=true,
              limits=(1.1 * minimum(xs1), 1.1 * maximum(xs1), 0, sim1.sol.t[end]))
    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{radiative}")
        lines!(ax, xnr, sim1.sol.t; linestyle=:dash, label=L"\textrm{non-radiative}")
        ax.xlabel = L"z\,(t\,)"
        ax.ylabel = L"t\,"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{radiative}")
        lines!(ax, sim1.sol.t, xnr; linestyle=:dash, label=L"\textrm{non-radiative}")
        ax.ylabel = L"z\,(t\,)"
        ax.xlabel = L"t\,"
    end
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    #rowsize!(fig.layout, 1, Relative(1.1))  # Let the axis size adjust
    #colsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    #rowgap!(fig.layout, 0)                 # Gap between plot and legend

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    # Save the figure
    savename = simple_name("particle_position")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_velocity(sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    vnr = zeros(nt)
    non_radiative_sol = non_radiative_reference_solution(sim1)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][3]
        vnr[i] = non_radiative_sol.u[i][3]
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1]; title=L"\textrm{Particle Velocity}")
    lines!(ax, sim1.sol.t, xs1; label=L"\textrm{radiative}")
    lines!(ax, sim1.sol.t, vnr; linestyle=:dash, label=L"\textrm{non-radiative}")
    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    # Save the figure
    savename = simple_name("particle_velocity")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_mass(sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    Nosc = sim1.params.Nosc
    k = sim1.params.k
    T = 2π / sqrt(abs(k))
    tf = Nosc * T
    L0 = 2tf
    x_range = (-L0 / 2, L0 / 2)

    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    m = zeros(nt)

    k = sim1.params.k
    q = sim1.params.q
    xs = sim1.params.x
    h = sim1.params.dx
    V(x) = InitialData.Harmonic(x, k)
    m0 = sim1.sol.u[1].x[2][1]
    mm = zeros(nt)
    for i in 1:nt
        z = sim1.sol.u[i].x[2][2]
        Q_field = sim1.sol.u[i].x[1]
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, sim1.sol.prob.p)

        mm[i] = m0 + q * Φrz
        m[i] = sim1.sol.u[i].x[2][1]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    # lines!(ax, sim1.sol.t, m; label=L"\textrm{particle}")
    lines!(ax, sim1.sol.t, m; label=L"\textrm{particle mass}")
    lines!(ax, sim1.sol.t, mm; label=L"\textrm{particle mass}")
    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"
    # leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    #rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    # Save the figure
    savename = simple_name("particle_mass")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_position_resolutions(sim1, sim2, sim4;
                                            invert=true, paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xs2 = zeros(nt)
    xs4 = zeros(nt)
    xnr = zeros(nt)
    non_radiative_sol = non_radiative_reference_solution(sim1)

    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][2]
        xs2[i] = sim2.sol.u[2i - 1].x[2][2]
        xs4[i] = sim4.sol.u[4i - 3].x[2][2]
        xnr[i] = non_radiative_sol.u[i][2]
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Position}")
    #,
    #          limits=(-0.025, 0.025, 0, sim1.sol.t[end]))
    palette = mytheme.palette
    color_ref = palette.color[][1]
    color_low = palette.color[][2]
    color_mid = palette.color[][3]
    color_high = palette.color[][4]
    low_style = :solid
    mid_style = :dash
    high_style = :dot

    # Plot the data
    if invert
        lines!(ax, xnr, sim1.sol.t; color=color_ref, linestyle=:solid,
               linewidth=1.2)
        lines!(ax, xs4, sim1.sol.t; color=color_high, linestyle=high_style, linewidth=1.4)
        lines!(ax, xs2, sim1.sol.t; color=color_mid, linestyle=mid_style, linewidth=1.4)
        lines!(ax, xs1, sim1.sol.t; color=color_low, linestyle=low_style, linewidth=1.6)

        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xnr; color=color_ref, linestyle=:solid, linewidth=1.2)
        lines!(ax, sim1.sol.t, xs4; color=color_high, linestyle=high_style, linewidth=1.4)
        lines!(ax, sim1.sol.t, xs2; color=color_mid, linestyle=mid_style, linewidth=1.4)
        lines!(ax, sim1.sol.t, xs1; color=color_low, linestyle=low_style, linewidth=1.6)
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end

    # # 2. Resolution legend
    res_legend = Legend(leg1[1, 1],
                        [LineElement(; color=color_low, linestyle=low_style),
                         LineElement(; color=color_mid, linestyle=mid_style),
                         LineElement(; color=color_high, linestyle=high_style),
                         LineElement(; color=color_ref, linestyle=:solid)],
                        [L"\textrm{low res}", L"\textrm{mid res}", L"\textrm{high res}",
                         L"\textrm{no radiation}"];
                        halign=:center,
                        valign=:top,
                        nbanks=2,
                        tellwidth=false,
                        tellheight=true, orientation=:horizontal)

    # Place legends below the axis
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Auto()) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Removed colsize! restriction to fix right spacing

    savename = simple_name("particle_position_res_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_velocity_resolutions(sim1, sim2, sim4; invert=true, paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xs2 = zeros(nt)
    xs4 = zeros(nt)
    vnr = zeros(nt)
    non_radiative_sol = non_radiative_reference_solution(sim1)

    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][3]
        xs2[i] = sim2.sol.u[2i - 1].x[2][3]
        xs4[i] = sim4.sol.u[4i - 3].x[2][3]
        vnr[i] = non_radiative_sol.u[i][3]
    end

    fig = Figure(; size=(253, 200))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Velocity}")
    palette = mytheme.palette
    color_ref = palette.color[][1]
    color_low = palette.color[][2]
    color_mid = palette.color[][3]
    color_high = palette.color[][4]
    low_style = :solid
    mid_style = :dash
    high_style = :dot

    # Plot the data
    lines!(ax, sim1.sol.t, vnr; color=color_ref, linestyle=:solid,
           linewidth=1.2)
    lines!(ax, sim1.sol.t, xs4; color=color_high, linestyle=high_style, linewidth=1.4)
    lines!(ax, sim1.sol.t, xs2; color=color_mid, linestyle=mid_style, linewidth=1.4)
    lines!(ax, sim1.sol.t, xs1; color=color_low, linestyle=low_style, linewidth=1.6)

    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"

    # # 2. Resolution legend (Low, Med, High Res)
    res_legend = Legend(leg1[1, 1],
                        [LineElement(; color=color_low, linestyle=low_style),
                         LineElement(; color=color_mid, linestyle=mid_style),
                         LineElement(; color=color_high, linestyle=high_style),
                         LineElement(; color=color_ref, linestyle=:solid)],
                        [L"\textrm{low res}", L"\textrm{mid res}", L"\textrm{high res}",
                         L"\textrm{no radiation}"];
                        halign=:center,
                        valign=:top,
                        tellwidth=false,
                        nbanks=2,
                        tellheight=true, orientation=:horizontal)

    # Place legends below the axis
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Auto()) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Removed colsize! restriction to fix right spacing

    savename = simple_name("particle_velocity_resolutions")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_self_force_resolutions(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    f1 = zeros(nt)
    f2 = zeros(nt)
    f4 = zeros(nt)

    for i in 1:nt
        f1[i] = relativistic_self_force_at_particle(sim1, i)
        f2[i] = relativistic_self_force_at_particle(sim2, 2i - 1)
        f4[i] = relativistic_self_force_at_particle(sim4, 4i - 3)
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Self Force}")
    palette = mytheme.palette
    color_zero = palette.color[][1]
    color_low = palette.color[][2]
    color_mid = palette.color[][3]
    color_high = palette.color[][4]
    low_style = :solid
    mid_style = :dash
    high_style = :dot

    lines!(ax, sim1.sol.t, zeros(nt); color=color_zero, linestyle=:solid, linewidth=1.2)
    lines!(ax, sim1.sol.t, f4; color=color_high, linestyle=high_style, linewidth=1.4)
    lines!(ax, sim1.sol.t, f2; color=color_mid, linestyle=mid_style, linewidth=1.4)
    lines!(ax, sim1.sol.t, f1; color=color_low, linestyle=low_style, linewidth=1.6)

    ax.ylabel = L"F_{\mathrm{sf}}(t)"
    ax.xlabel = L"t"
    apply_scientific_notation_if_needed!(ax)

    Legend(leg1[1, 1],
           [LineElement(; color=color_low, linestyle=low_style),
            LineElement(; color=color_mid, linestyle=mid_style),
            LineElement(; color=color_high, linestyle=high_style),
            LineElement(; color=color_zero, linestyle=:solid)],
           [L"\textrm{low res}", L"\textrm{mid res}", L"\textrm{high res}", L"0"];
           halign=:center,
           valign=:top,
           tellwidth=false,
           nbanks=2,
           tellheight=true,
           orientation=:horizontal)

    rowsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Auto())
    rowgap!(fig.layout, 0)

    savename = simple_name("self_force_resolutions")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_position_diff_scaled(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    nt = length(sim1.sol.t)

    diff1 = zeros(nt)
    diff2_scaled = zeros(nt)
    for i in 1:nt
        y1 = sim1.sol.u[i].x[2][2]
        y2 = sim2.sol.u[2i - 1].x[2][2]
        y4 = sim4.sol.u[4i - 3].x[2][2]

        diff1[i] = abs(y1 - y2)
        diff2_scaled[i] = 2^1 * abs(y2 - y4) # Assuming 2nd order
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    # No colsize restriction
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Delta z",
              title=L"\textrm{Scalled Position Error}")

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]

    lines!(ax, sim1.sol.t, diff1; label=L"|z_{\mathrm{low}} - z_{\mathrm{mid}}|",
           color=color1)
    lines!(ax, sim1.sol.t, diff2_scaled;
           label=L"2 \cdot |z_{\mathrm{mid}} - z_{\mathrm{high}}|", color=color2,
           linestyle=:dash)

    # apply_scientific_notation_if_needed!(ax)
    axislegend(ax; position=:lt, labelsize=10, rowgap=0.2) # Top-left legend
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    savename = simple_name("particle_position_diff_scaled_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_velocity_diff_scaled(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    nt = length(sim1.sol.t)

    diff1 = zeros(nt)
    diff2_scaled = zeros(nt)
    for i in 1:nt
        y1 = sim1.sol.u[i].x[2][3]
        y2 = sim2.sol.u[2i - 1].x[2][3]
        y4 = sim4.sol.u[4i - 3].x[2][3]

        diff1[i] = abs(y1 - y2)
        diff2_scaled[i] = 2^1 * abs(y2 - y4) # Assuming 2nd order
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    # No colsize restriction
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Delta v",
              title=L"\textrm{Scaled Velocity Error}")

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]

    lines!(ax, sim1.sol.t, diff1; label=L"|v_{\mathrm{low}} - v_{\mathrm{mid}}|",
           color=color1)
    lines!(ax, sim1.sol.t, diff2_scaled;
           label=L"2 \cdot |v_{\mathrm{mid}} - v_{\mathrm{high}}|", color=color2,
           linestyle=:dash)

    # apply_scientific_notation_if_needed!(ax)
    axislegend(ax; position=:lt, labelsize=10, rowgap=0.1)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    savename = simple_name("particle_velocity_diff_scaled_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_position_self_convergence_order(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    nt = length(sim1.sol.t)

    convfac = zeros(nt)
    for i in 1:nt
        y1 = sim1.sol.u[i].x[2][2]
        y2 = sim2.sol.u[2i - 1].x[2][2]
        y4 = sim4.sol.u[4i - 3].x[2][2]

        diff1 = abs(y1 - y2)
        diff2 = abs(y2 - y4)
        convfac[i] = diff1 / (diff2 + eps(Float64))
    end
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              title=L"\textrm{Convergence order of } z(t)")
    # limits=(0, sim1.sol.t[end], 0, 4))
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    savename = simple_name("particle_position_convord_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_velocity_self_convergence_order(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    nt = length(sim1.sol.t)

    convfac = zeros(nt)
    for i in 1:nt
        y1 = sim1.sol.u[i].x[2][3]
        y2 = sim2.sol.u[2i - 1].x[2][3]
        y4 = sim4.sol.u[4i - 3].x[2][3]

        diff1 = abs(y1 - y2)
        diff2 = abs(y2 - y4)
        convfac[i] = diff1 / (diff2 + eps(Float64))
    end
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              title=L"\textrm{Convergence order of } v(t)")
    # limits=(0, sim1.sol.t[end], 0, 4))
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    savename = simple_name("particle_velocity_convord_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_self_force_self_convergence_order(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    nt = length(sim1.sol.t)

    convfac = zeros(nt)
    for i in 1:nt
        f1 = relativistic_self_force_at_particle(sim1, i)
        f2 = relativistic_self_force_at_particle(sim2, 2i - 1)
        f4 = relativistic_self_force_at_particle(sim4, 4i - 3)

        diff1 = abs(f1 - f2)
        diff2 = abs(f2 - f4)
        convfac[i] = diff1 / (diff2 + eps(Float64))
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              title=L"\textrm{Convergence order of } F_{\mathrm{sf}}(t)")
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)

    rowsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("self_force_convord_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_self_force_diff_scaled(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    nt = length(sim1.sol.t)

    diff1 = zeros(nt)
    diff2_scaled = zeros(nt)
    for i in 1:nt
        f1 = relativistic_self_force_at_particle(sim1, i)
        f2 = relativistic_self_force_at_particle(sim2, 2i - 1)
        f4 = relativistic_self_force_at_particle(sim4, 4i - 3)

        diff1[i] = abs(f1 - f2)
        diff2_scaled[i] = 2^1 * abs(f2 - f4)
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Delta F_{\mathrm{sf}}",
              title=L"\textrm{Scaled Self-Force Error}")

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]

    lines!(ax, sim1.sol.t, diff1;
           label=L"|F_{\mathrm{sf},\mathrm{low}} - F_{\mathrm{sf},\mathrm{mid}}|",
           color=color1)
    lines!(ax, sim1.sol.t, diff2_scaled;
           label=L"2 \cdot |F_{\mathrm{sf},\mathrm{mid}} - F_{\mathrm{sf},\mathrm{high}}|",
           color=color2, linestyle=:dash)

    axislegend(ax; position=:lt, labelsize=10, rowgap=0.1)
    rowsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("self_force_diff_scaled_potential")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## 
function plot_energy(sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    tf = sim1.sol.t[end]

    mask = (sim1.params.x .>= -tf) .& (sim1.params.x .<= tf)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)

    dt = sim1.params.dt
    q = sim1.params.q
    h = sim1.params.dx
    m0 = 1.0

    k = sim1.params.k
    V(x) = InitialData.Harmonic(x, k)
    ts2plot = sim1.sol.t
    field_energy = zeros(nt)
    potential_energy = zeros(nt)
    kinetic_energy = zeros(nt)
    boundary_flux = zeros(nt)
    particle_cross_term = zeros(nt)
    self_energy = zeros(nt)
    for i in 1:nt
        m = sim1.sol.u[i].x[2][1]
        z = sim1.sol.u[i].x[2][2]
        v = sim1.sol.u[i].x[2][3]
        Q_field = sim1.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, sim1.sol.prob.p)

        Πs = ScalarField.Πs.(x_filtered, z, v, q)
        Ψs = ScalarField.Ψs.(x_filtered, z, v, q)
        Πtot = Πr .+ Πs
        Ψtot = Ψr .+ Ψs

        Hf = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs
        Hfr = @. 0.5 * (Πr^2 + Ψr^2)
        field_energy[i] = trapezoidal(Hf, h) #- q^2 / 8
        potential_energy[i] = V(z)
        # With the repo convention Π = ∂ₜΦ and signature (-,+), the energy flux is
        # T^{01} = -ΠΨ, so the boundary contribution is T^{01}(x_R)-T^{01}(x_L).
        boundary_flux[i] = -Πtot[end] * Ψtot[end] + Πtot[1] * Ψtot[1]

        γ = inv(sqrt(1.0 - v * v))
        kinetic_energy[i] = (m0 - q * Φrz) * (γ - 1.0)
        self_energy[i] = 0#-2q * Φrz * γ
        particle_cross_term[i] = particle_singular_cross_term(sim1, i)
    end
    #cumulative_boundary_flux = cumulative_trapezoidal(ts2plot, dt, boundary_flux)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Energy}")

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]
    color3 = palette.color[][3]
    color4 = palette.color[][4]
    color5 = palette.color[][5]
    color6 = palette.color[][6]

    # Plot the data
    # potential_energy /= max(potential_energy...)
    # kinetic_energy /= max(kinetic_energy...)
    # field_energy /= max(field_energy...)
    total_energy = field_energy .+ potential_energy .+ kinetic_energy .+ self_energy
    # total_energy /= max(total_energy...)
    # cumulative_boundary_flux /= max(cumulative_boundary_flux...)

    # ylims!(ax, ymin - 0.1 * yrange, ymax + 0.1 * yrange)

    lines!(ax, ts2plot, field_energy; label=L"E_{\mathrm{field}}", color=color1,
           linewidth=1.2)
    # lines!(ax, ts2plot, potential_energy; label=L"E_{\mathrm{potential}}", color=color2,
    #        linewidth=1.2)
    # lines!(ax, ts2plot, kinetic_energy; label=L"E_{\mathrm{kinetic}}", color=color3,
    #        linewidth=1.2)
    # # lines!(ax, ts2plot, self_energy; label=L"E_{\mathrm{self}}", color=color4,
    #        linewidth=1.2)
    lines!(ax, ts2plot, total_energy;
           label=L"E_{\mathrm{total}}", color=color4, linewidth=1.2)
    # # lines!(ax, ts2plot, particle_cross_term;
    #        label=L"\Psi_r(z)\bar{\Psi}_s-\Pi_r(z)\bar{\Pi}_s",
    #        color=color5, linewidth=1.2, linestyle=:dash)
    # # lines!(ax, ts2plot[1:(end - 30)], boundary_flux[1:(end - 30)];
    #        label=L"\mathcal{F}", color=color5, linewidth=1.2)
    # lines!(ax, ts2plot, corrected_total_energy;
    #        label=L"E_{\mathrm{total}}+\mathcal{F}", color=color6, linewidth=1.2,
    #        linestyle=:dash)

    # ax.ylabel = L""
    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, nbanks=2,
                 tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.16)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    # Save the figure
    savename = simple_name("energy_conservation")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_raw_energy_balance(sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    ts = sim1.sol.t
    tf = ts[end]
    mask = (sim1.params.x .>= -tf) .& (sim1.params.x .<= tf)
    x_filtered = sim1.params.x[mask]

    nt = length(ts)
    q = sim1.params.q
    h = sim1.params.dx
    k = sim1.params.k
    V(x) = InitialData.Harmonic(x, k)

    particle_energy = zeros(nt)
    field_energy = zeros(nt)
    inside_energy = zeros(nt)
    outward_flux = zeros(nt)

    for i in 1:nt
        m = sim1.sol.u[i].x[2][1]
        z = sim1.sol.u[i].x[2][2]
        v = sim1.sol.u[i].x[2][3]

        Πr = @view sim1.sol.u[i].x[1][mask, 1]
        Ψr = @view sim1.sol.u[i].x[1][mask, 2]
        Πs = ScalarField.Πs.(x_filtered, z, v, q)
        Ψs = ScalarField.Ψs.(x_filtered, z, v, q)

        # The raw finite-domain balance uses the physical total field, not only the
        # regular piece, so we reconstruct Π = Πr + Πs and Ψ = Ψr + Ψs everywhere.
        Πtot = Πr .+ Πs
        Ψtot = Ψr .+ Ψs

        Φrz = ODE.gather_muscl(@view(sim1.sol.u[i].x[1][:, 3]), z, sim1.sol.prob.p)

        m0 = hasproperty(sim1.params, :m0) ? sim1.params.m0 : 1.0
        γ = inv(sqrt(1.0 - v * v))
        particle_energy[i] = (m0 - q * Φrz) * (γ - 1.0) + V(z)

        Hfield = @. 0.5 * (-Πtot^2 + Ψtot^2)
        field_energy[i] = trapezoidal(Hfield, h)
        inside_energy[i] = particle_energy[i] + field_energy[i]

        # With Π = ∂ₜΦ and signature (-,+), T^{01} = -ΠΨ. On a finite interval,
        # energy that leaves through xL/xR must be accounted for explicitly via
        # T^{01}(x_R)-T^{01}(x_L).
        outward_flux[i] = -Πtot[end] * Ψtot[end] + Πtot[1] * Ψtot[1]
    end

    cumulative_flux = cumulative_trapezoidal(ts, outward_flux)

    # This residual measures the discrete interior balance:
    # change of energy inside the domain + accumulated outward flux.
    balance_residual = inside_energy .- inside_energy[1] .+ cumulative_flux
    inside_energy_change = inside_energy .- inside_energy[1]

    fig = Figure(; size=(253, 360), figure_padding=(1, 5, 5, 1))
    ax_energy = Axis(fig[1, 1];
                     title=L"\textrm{Raw Interior Energy}",
                     ylabel=L"E")
    ax_balance = Axis(fig[3, 1];
                      title=L"\textrm{Finite-Domain Balance}",
                      xlabel=L"t",
                      ylabel=L"\Delta E,\ I_{\mathrm{flux}},\ R")

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]
    color3 = palette.color[][3]
    color4 = palette.color[][4]
    color5 = palette.color[][5]
    color6 = palette.color[][6]

    lines!(ax_energy, ts, inside_energy; linestyle=:dot, label=L"E_{\mathrm{inside}}",
           color=color1,
           linewidth=1.2)
    lines!(ax_energy, ts, particle_energy; linestyle=:dot, label=L"E_{\mathrm{particle}}",
           color=color2,
           linewidth=1.2)
    lines!(ax_energy, ts, field_energy; linestyle=:dot, label=L"E_{\mathrm{field}}",
           color=color3,
           linewidth=1.2)

    lines!(ax_balance, ts, inside_energy_change;
           label=L"E_{\mathrm{inside}}(t)-E_{\mathrm{inside}}(0)", color=color4,
           linewidth=1.2)
    lines!(ax_balance, ts, cumulative_flux; label=L"I_{\mathrm{flux}}", color=color5,
           linewidth=1.2)
    lines!(ax_balance, ts, balance_residual; label=L"R_{\mathrm{balance}}", color=color6,
           linewidth=1.2)

    apply_scientific_notation_if_needed!(ax_energy)
    apply_scientific_notation_if_needed!(ax_balance)

    leg_energy = Legend(fig[2, 1], ax_energy; orientation=:horizontal,
                        tellwidth=false, tellheight=true)
    leg_balance = Legend(fig[4, 1], ax_balance; orientation=:horizontal,
                         tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.1))
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Relative(0.1))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("raw_energy_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_fluxes(r, sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    F = 2
    # L0 = sim1.params.L0
    x_range = (-L / F, L / F)
    @show x_range
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)

    flux_l = zeros(nt)
    flux_r = zeros(nt)

    q = sim1.params.q
    sf = sim1.params.sf
    h = sim1.params.dx
    xs = sim1.params.x

    for i in 1:nt
        m = sim1.sol.u[i].x[2][1]
        z = sim1.sol.u[i].x[2][2]
        v = sim1.sol.u[i].x[2][3]

        Πr = @view sim1.sol.u[i].x[1][mask, 1]
        Ψr = @view sim1.sol.u[i].x[1][mask, 2]
        Πr_int = PCHIPInterpolation(Πr, x_filtered)
        Ψr_int = PCHIPInterpolation(Ψr, x_filtered)

        flux_r[i] = -Πr_int(r) * Ψr_int(r)
        flux_l[i] = Πr_int(-r) * Ψr_int(-r)
    end

    fig = Figure(; size=(2 * 253, 2 * 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Fluxes}")
    # Plot the data
    lines!(ax, sim1.sol.t, flux_r; label="right")
    lines!(ax, sim1.sol.t, flux_l; label="left")
    lines!(ax, sim1.sol.t, flux_r + flux_l; label="total")

    ax.ylabel = L"fluxes"
    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    # Save the figure
    savename = simple_name("fluxes")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_work_balance_laws(sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    ts = sim1.sol.t
    tf = ts[end]
    boundary = tf / 3
    mask = (sim1.params.x .>= -boundary) .& (sim1.params.x .<= boundary)
    if iszero(tf)
        mask = trues(length(sim1.params.x))
    end
    x_filtered = sim1.params.x[mask]

    nt = length(ts)
    q = sim1.params.q
    h = sim1.params.dx
    k = sim1.params.k
    dt = sim1.params.dt
    V(x) = InitialData.Harmonic(x, k)

    v_arr = [sim1.sol.u[i].x[2][3] for i in 1:nt]
    v_interp = PCHIPInterpolation(v_arr, ts)
    a_arr = [DataInterpolations.derivative(v_interp, t) for t in ts]

    particle_energy = zeros(nt)
    field_energy = zeros(nt)
    inside_energy = zeros(nt)
    outward_flux = zeros(nt)
    W_src = zeros(nt)
    m0 = hasproperty(sim1.params, :m0) ? sim1.params.m0 : 1.0
    for i in 1:nt
        # m = sim1.sol.u[i].x[2][1]
        z = sim1.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        Πs = ScalarField.Πs.(x_filtered, z, v, q)
        Ψs = ScalarField.Ψs.(x_filtered, z, v, q)
        Πr = @view sim1.sol.u[i].x[1][mask, 1]
        Ψr = @view sim1.sol.u[i].x[1][mask, 2]
        Φr = @view sim1.sol.u[i].x[1][mask, 3]

        Φrz = ODE.gather_muscl(@view(sim1.sol.u[i].x[1][:, 3]), z, sim1.sol.prob.p)

        γ = inv(sqrt(1.0 - v * v))

        particle_energy[i] = (m0 - q * Φrz) * (γ - 1.0) + V(z)

        Hfield = @. 0.5 * (Πr^2 + Ψr^2) #+ Πr * Πs + Ψr * Ψs
        field_energy[i] = trapezoidal(Hfield, h)
        inside_energy[i] = particle_energy[i] + field_energy[i]

        outward_flux[i] = -Πr[end] * Ψr[end] + Πr[1] * Ψr[1]

        W_sum = 0.0

        for (idx, x_i) in enumerate(x_filtered)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            W_sum += Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end
        W_src[i] = W_sum * h
    end

    E_tot_interp = PCHIPInterpolation(inside_energy, ts)
    dE_tot_dt = [DataInterpolations.derivative(E_tot_interp, t) for t in ts]

    rhs_rate = -outward_flux .+ W_src
    rate_residual = dE_tot_dt .- rhs_rate

    cumulative_flux = cumulative_trapezoidal(ts, dt, outward_flux)
    cumulative_work = cumulative_trapezoidal(ts, dt, W_src)

    inside_energy_change = inside_energy .- inside_energy[1]
    rhs_int = -cumulative_flux .+ cumulative_work
    int_residual = inside_energy_change .- rhs_int

    fig = Figure(; size=(253, 550), figure_padding=(1, 5, 5, 1))

    # 1. Plot of W_src
    ax_work = Axis(fig[1, 1];
                   title=L"\textrm{Work Term}",
                   ylabel=L"W_{\mathrm{src},h}")
    lines!(ax_work, ts, W_src; label=L"W_{\mathrm{src},h}",
           color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax_work)

    # 2. Rate Form Balance
    ax_rate = Axis(fig[3, 1];
                   title=L"\textrm{Rate Form Balance}",
                   ylabel=L"\mathrm{Rate}")
    lines!(ax_rate, ts, dE_tot_dt; label=L"\frac{d(E_{\mathrm{part},h} + E_{r,h})}{dt}",
           color=mytheme.palette.color[][1], linewidth=1.5)
    lines!(ax_rate, ts, rhs_rate;
           label=L"-\mathcal{F}_{\mathrm{bdry},h} + W_{\mathrm{src},h}",
           color=mytheme.palette.color[][2], linestyle=:dash, linewidth=1.5)
    lines!(ax_rate, ts, rate_residual; label=L"\mathrm{Residual}",
           color=mytheme.palette.color[][4], linewidth=1.2, linestyle=:dot)
    apply_scientific_notation_if_needed!(ax_rate)

    # 3. Integral Form Balance
    ax_int = Axis(fig[5, 1];
                  title=L"\textrm{Integral Form Balance}",
                  xlabel=L"t",
                  ylabel=L"\mathrm{Energy}")
    lines!(ax_int, ts, inside_energy_change; label=L"\Delta(E_{\mathrm{part},h} + E_{r,h})",
           color=mytheme.palette.color[][1], linewidth=1.5)
    lines!(ax_int, ts, rhs_int; label=L"-E_{\mathrm{flux},h} + \int W_{\mathrm{src},h} dt",
           color=mytheme.palette.color[][2], linestyle=:dash, linewidth=1.5)
    lines!(ax_int, ts, int_residual; label=L"\mathrm{Residual}",
           color=mytheme.palette.color[][4], linewidth=1.2, linestyle=:dot)
    apply_scientific_notation_if_needed!(ax_int)

    # Legends
    Legend(fig[2, 1], ax_work; orientation=:horizontal, tellwidth=false, tellheight=true)
    Legend(fig[4, 1], ax_rate; orientation=:horizontal, tellwidth=false, tellheight=true)
    Legend(fig[6, 1], ax_int; orientation=:horizontal, tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.08))
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Relative(0.08))
    rowsize!(fig.layout, 5, Auto())
    rowsize!(fig.layout, 6, Relative(0.08))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("work_balance_laws")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_field_energy_split_balance(sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    ts = sim1.sol.t
    tf = ts[end]
    x = sim1.params.x
    boundary = tf / 3
    mask = (x .>= -boundary) .& (x .<= boundary)
    if iszero(tf) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    nt = length(ts)
    q = sim1.params.q
    h = sim1.params.dx
    params = sim1.sol.prob.p
    equation = params.equation
    slope_limiter = params.slope_limiter
    N = params.N
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)
    v_arr = [sim1.sol.u[i].x[2][3] for i in 1:nt]
    v_interp = PCHIPInterpolation(v_arr, ts)
    a_arr = [DataInterpolations.derivative(v_interp, t) for t in ts]

    E_r = zeros(nt)
    E_corr = zeros(nt)
    E_tot = zeros(nt)

    outward_flux_r = zeros(nt)
    outward_flux_corr = zeros(nt)
    outward_flux_tot = zeros(nt)
    exchange_rate = zeros(nt)

    for i in 1:nt
        z = sim1.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        Q_field = sim1.sol.u[i].x[1]

        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        Πs = similar(Πr)
        Ψs = similar(Ψr)
        source_density = similar(Πr)

        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end

        e_r = @. 0.5 * (Πr^2 + Ψr^2)
        e_cross = @. Πr * Πs + Ψr * Ψs
        e_s = @. 0.5 * (Πs^2 + Ψs^2)
        e_corr = e_s .+ e_cross

        E_r[i] = h * sum(e_r)
        E_corr[i] = h * sum(e_corr)
        E_tot[i] = E_r[i] + E_corr[i]
        exchange_rate[i] = h * sum(source_density)

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, N, L, R, slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, N, L, R, slope_limiter)

        Πr_left_face = q_left_minus[1]
        Ψr_left_face = q_left_minus[2]
        Πr_right_face = q_right_plus[1]
        Ψr_right_face = q_right_plus[2]

        left_face_x = x[left_idx] - 0.5 * h
        right_face_x = x[right_idx] + 0.5 * h
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x, z, v, q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x, z, v,
                                                                       q,
                                                                       -one(v))

        outward_flux_r[i] = -Πr_right_face * Ψr_right_face + Πr_left_face * Ψr_left_face
        outward_flux_corr[i] = -(Πs_right_face * Ψs_right_face +
                                 Πr_right_face * Ψs_right_face +
                                 Πs_right_face * Ψr_right_face) +
                               (Πs_left_face * Ψs_left_face +
                                Πr_left_face * Ψs_left_face +
                                Πs_left_face * Ψr_left_face)

        Πtot_left_face = Πr_left_face + Πs_left_face
        Ψtot_left_face = Ψr_left_face + Ψs_left_face
        Πtot_right_face = Πr_right_face + Πs_right_face
        Ψtot_right_face = Ψr_right_face + Ψs_right_face
        outward_flux_tot[i] = -Πtot_right_face * Ψtot_right_face +
                              Πtot_left_face * Ψtot_left_face
    end

    E_r_interp = PCHIPInterpolation(E_r, ts)
    E_corr_interp = PCHIPInterpolation(E_corr, ts)
    E_tot_interp = PCHIPInterpolation(E_tot, ts)
    dE_r_dt = [DataInterpolations.derivative(E_r_interp, t) for t in ts]
    dE_corr_dt = [DataInterpolations.derivative(E_corr_interp, t) for t in ts]
    dE_tot_dt = [DataInterpolations.derivative(E_tot_interp, t) for t in ts]

    rhs_regular = -outward_flux_r .+ exchange_rate
    rhs_correction = -outward_flux_corr .- exchange_rate
    rhs_total = -outward_flux_tot

    residual_regular = dE_r_dt .- rhs_regular
    residual_correction = dE_corr_dt .- rhs_correction
    residual_total = dE_tot_dt .- rhs_total

    fig = Figure(; size=(560, 620), figure_padding=(6, 8, 8, 6))

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]
    color3 = palette.color[][3]
    color4 = palette.color[][4]
    color5 = palette.color[][5]

    panel11 = GridLayout(fig[1, 1])
    panel12 = GridLayout(fig[1, 2])
    panel21 = GridLayout(fig[3, 1])
    panel22 = GridLayout(fig[3, 2])

    ax_energy = Axis(panel11[1, 1];
                     title=L"\textrm{Field-Energy Split}",
                     ylabel=L"E")
    lines!(ax_energy, ts, E_r; label=L"E_r",
           color=color1, linewidth=1.5)
    lines!(ax_energy, ts, E_corr; label=L"E_{\mathrm{corr}}",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_energy, ts, E_tot; label=L"E_{\mathrm{tot}}",
           color=color3, linewidth=1.4)

    ax_regular = Axis(panel12[1, 1];
                      title=L"\textrm{Regular Balance}",
                      ylabel=L"\mathrm{Rate}")
    lines!(ax_regular, ts, dE_r_dt;
           label=L"\frac{dE_r}{dt}",
           color=color1, linewidth=1.5)
    lines!(ax_regular, ts, rhs_regular;
           label=L"-\mathcal{F}_{r} + \int_{\Omega}(-\Pi_r R_2 - \Psi_r R_3)\,dx",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_regular, ts, residual_regular;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    ax_correction = Axis(panel21[1, 1];
                         title=L"\textrm{Correction Balance}",
                         xlabel=L"t",
                         ylabel=L"\mathrm{Rate}")
    lines!(ax_correction, ts, dE_corr_dt;
           label=L"\frac{dE_{\mathrm{corr}}}{dt}",
           color=color1, linewidth=1.5)
    lines!(ax_correction, ts, rhs_correction;
           label=L"-\mathcal{F}_{\mathrm{corr}} - \int_{\Omega}(-\Pi_r R_2 - \Psi_r R_3)\,dx",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_correction, ts, residual_correction;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    ax_total = Axis(panel22[1, 1];
                    title=L"\textrm{Total Balance}",
                    xlabel=L"t",
                    ylabel=L"\mathrm{Rate}")
    lines!(ax_total, ts, dE_tot_dt;
           label=L"\frac{dE_{\mathrm{tot}}}{dt}",
           color=color1, linewidth=1.5)
    lines!(ax_total, ts, rhs_total;
           label=L"-\mathcal{F}_{\mathrm{tot}}",
           color=color5, linestyle=:dash, linewidth=1.5)
    lines!(ax_total, ts, residual_total;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    apply_scientific_notation_if_needed!(ax_energy)
    apply_scientific_notation_if_needed!(ax_regular)
    apply_scientific_notation_if_needed!(ax_correction)
    apply_scientific_notation_if_needed!(ax_total)

    Legend(panel11[2, 1], ax_energy; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=3)
    Legend(panel12[2, 1], ax_regular; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=1)
    Legend(panel21[2, 1], ax_correction; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=1)
    Legend(panel22[2, 1], ax_total; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=1)

    rowsize!(panel11, 1, Auto())
    rowsize!(panel11, 2, Relative(0.14))
    rowsize!(panel12, 1, Auto())
    rowsize!(panel12, 2, Relative(0.14))
    rowsize!(panel21, 1, Auto())
    rowsize!(panel21, 2, Relative(0.14))
    rowsize!(panel22, 1, Auto())
    rowsize!(panel22, 2, Relative(0.14))

    colsize!(fig.layout, 1, Relative(0.5))
    colsize!(fig.layout, 2, Relative(0.5))
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Fixed(8))
    rowsize!(fig.layout, 3, Auto())
    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 0)

    savename = simple_name("field_energy_split_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_field_energy_balance(sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    ts = sim1.sol.t
    tf = ts[end]
    x = sim1.params.x
    boundary = tf / 3
    mask = (x .>= -boundary) .& (x .<= boundary)
    if iszero(tf) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    nt = length(ts)
    q = sim1.params.q
    h = sim1.params.dx
    m0 = hasproperty(sim1.params, :m) ? sim1.params.m : 1.0
    k = sim1.params.k
    V(x) = InitialData.Harmonic(x, k)
    params = sim1.sol.prob.p
    equation = params.equation
    slope_limiter = params.slope_limiter
    N = params.N
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)
    v_arr = [sim1.sol.u[i].x[2][3] for i in 1:nt]
    v_interp = PCHIPInterpolation(v_arr, ts)
    a_arr = [DataInterpolations.derivative(v_interp, t) for t in ts]

    E_p = zeros(nt)
    E_r = zeros(nt)
    E_corr = zeros(nt)
    E_pr = zeros(nt)
    E_prcorr = zeros(nt)

    particle_power = zeros(nt)
    outward_flux_r = zeros(nt)
    outward_flux_corr = zeros(nt)
    outward_flux_tot = zeros(nt)
    exchange_rate = zeros(nt)

    for i in 1:nt
        z = sim1.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        Q_field = sim1.sol.u[i].x[1]

        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        _, Ψrz = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)

        E_p[i] = 0.5 * m0 * v * v + V(z)
        particle_power[i] = q * v * Ψrz

        Πs = similar(Πr)
        Ψs = similar(Ψr)
        source_density = similar(Πr)

        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end

        e_r = @. 0.5 * (Πr^2 + Ψr^2)
        e_cross = @. Πr * Πs + Ψr * Ψs
        e_s = @. 0.5 * (Πs^2 + Ψs^2)
        e_corr = e_s .+ e_cross

        E_r[i] = h * sum(e_r)
        E_corr[i] = h * sum(e_corr)
        E_pr[i] = E_p[i] + E_r[i]
        E_prcorr[i] = E_pr[i] + E_corr[i]
        exchange_rate[i] = h * sum(source_density)

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, N, L, R, slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, N, L, R, slope_limiter)

        Πr_left_face = q_left_minus[1]
        Ψr_left_face = q_left_minus[2]
        Πr_right_face = q_right_plus[1]
        Ψr_right_face = q_right_plus[2]

        left_face_x = x[left_idx] - 0.5 * h
        right_face_x = x[right_idx] + 0.5 * h
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x, z, v, q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x, z, v,
                                                                       q,
                                                                       -one(v))

        outward_flux_r[i] = -Πr_right_face * Ψr_right_face + Πr_left_face * Ψr_left_face
        outward_flux_corr[i] = -(Πs_right_face * Ψs_right_face +
                                 Πr_right_face * Ψs_right_face +
                                 Πs_right_face * Ψr_right_face) +
                               (Πs_left_face * Ψs_left_face +
                                Πr_left_face * Ψs_left_face +
                                Πs_left_face * Ψr_left_face)

        Πtot_left_face = Πr_left_face + Πs_left_face
        Ψtot_left_face = Ψr_left_face + Ψs_left_face
        Πtot_right_face = Πr_right_face + Πs_right_face
        Ψtot_right_face = Ψr_right_face + Ψs_right_face
        outward_flux_tot[i] = -Πtot_right_face * Ψtot_right_face +
                              Πtot_left_face * Ψtot_left_face
    end

    E_p_interp = PCHIPInterpolation(E_p, ts)
    E_pr_interp = PCHIPInterpolation(E_pr, ts)
    E_prcorr_interp = PCHIPInterpolation(E_prcorr, ts)
    dE_p_dt = [DataInterpolations.derivative(E_p_interp, t) for t in ts]
    dE_pr_dt = [DataInterpolations.derivative(E_pr_interp, t) for t in ts]
    dE_prcorr_dt = [DataInterpolations.derivative(E_prcorr_interp, t) for t in ts]

    rhs_particle = particle_power
    rhs_pr = particle_power .- outward_flux_r .+ exchange_rate
    lhs_prcorr_remainder = dE_prcorr_dt .+ outward_flux_tot
    rhs_prcorr_remainder = particle_power

    residual_particle = dE_p_dt .- rhs_particle
    residual_pr = dE_pr_dt .- rhs_pr
    residual_prcorr = lhs_prcorr_remainder .- rhs_prcorr_remainder

    fig = Figure(; size=(560, 620), figure_padding=(6, 8, 8, 6))

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]
    color3 = palette.color[][3]
    color4 = palette.color[][4]
    color5 = palette.color[][5]

    panel11 = GridLayout(fig[1, 1])
    panel12 = GridLayout(fig[1, 2])
    panel21 = GridLayout(fig[3, 1])
    panel22 = GridLayout(fig[3, 2])

    ax_energy = Axis(panel11[1, 1];
                     title=L"\textrm{Particle + Field Energies}",
                     ylabel=L"E")
    lines!(ax_energy, ts, E_p; label=L"E_p = \frac12 m_0 v^2 + V(z)",
           color=color1, linewidth=1.5)
    lines!(ax_energy, ts, E_r; label=L"E_r",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_energy, ts, E_corr; label=L"E_{\mathrm{corr}}",
           color=color3, linestyle=:dot, linewidth=1.5)
    lines!(ax_energy, ts, E_prcorr; label=L"E_p + E_r + E_{\mathrm{corr}}",
           color=color5, linewidth=1.4)

    ax_particle = Axis(panel12[1, 1];
                       title=L"\textrm{Particle Balance}",
                       ylabel=L"\mathrm{Rate}")
    lines!(ax_particle, ts, dE_p_dt;
           label=L"\frac{dE_p}{dt}",
           color=color1, linewidth=1.5)
    lines!(ax_particle, ts, rhs_particle;
           label=L"q v \Psi_r(z)",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_particle, ts, residual_particle;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    ax_pr = Axis(panel21[1, 1];
                 title=L"\textrm{Particle + Regular Balance}",
                 xlabel=L"t",
                 ylabel=L"\mathrm{Rate}")
    lines!(ax_pr, ts, dE_pr_dt;
           label=L"\frac{d(E_p+E_r)}{dt}",
           color=color1, linewidth=1.5)
    lines!(ax_pr, ts, rhs_pr;
           label=L"q v \Psi_r(z) - \mathcal{F}_r + \int_{\Omega}(-\Pi_r R_2 - \Psi_r R_3)\,dx",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_pr, ts, residual_pr;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    ax_prcorr = Axis(panel22[1, 1];
                     title=L"\textrm{After Adding } E_{\mathrm{corr}}",
                     xlabel=L"t",
                     ylabel=L"\mathrm{Rate}")
    lines!(ax_prcorr, ts, lhs_prcorr_remainder;
           label=L"\frac{d(E_p+E_r+E_{\mathrm{corr}})}{dt} + \mathcal{F}_{\mathrm{tot}}",
           color=color1, linewidth=1.5)
    lines!(ax_prcorr, ts, rhs_prcorr_remainder;
           label=L"q v \Psi_r(z)",
           color=color3, linestyle=:dash, linewidth=1.5)
    lines!(ax_prcorr, ts, residual_prcorr;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    apply_scientific_notation_if_needed!(ax_energy)
    apply_scientific_notation_if_needed!(ax_particle)
    apply_scientific_notation_if_needed!(ax_pr)
    apply_scientific_notation_if_needed!(ax_prcorr)

    Legend(panel11[2, 1], ax_energy; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=2)
    Legend(panel12[2, 1], ax_particle; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=1)
    Legend(panel21[2, 1], ax_pr; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=1)
    Legend(panel22[2, 1], ax_prcorr; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=1)

    rowsize!(panel11, 1, Auto())
    rowsize!(panel11, 2, Relative(0.14))
    rowsize!(panel12, 1, Auto())
    rowsize!(panel12, 2, Relative(0.14))
    rowsize!(panel21, 1, Auto())
    rowsize!(panel21, 2, Relative(0.14))
    rowsize!(panel22, 1, Auto())
    rowsize!(panel22, 2, Relative(0.14))

    colsize!(fig.layout, 1, Relative(0.5))
    colsize!(fig.layout, 2, Relative(0.5))
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Fixed(8))
    rowsize!(fig.layout, 3, Auto())
    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 0)

    savename = simple_name("particle_field_energy_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_field_balance(sim1; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    ts = sim1.sol.t
    x = sim1.params.x

    nt = length(ts)
    q = sim1.params.q
    h = sim1.params.dx
    boundary = 0.18 / 2
    mask = (x .>= -boundary) .& (x .<= boundary)
    if !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    params = sim1.sol.prob.p
    equation = params.equation
    slope_limiter = params.slope_limiter
    N = params.N
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)

    E_field_tot = zeros(nt)
    particle_energy = zeros(nt)
    E_total = zeros(nt)
    outward_flux = zeros(nt)
    particle_inside_mask = falses(nt)
    m0 = hasproperty(sim1.params, :m0) ? sim1.params.m0 : 1.0
    k = sim1.params.k
    V(x) = InitialData.Harmonic(x, k)

    for i in 1:nt
        z = sim1.sol.u[i].x[2][2]
        v = sim1.sol.u[i].x[2][3]
        Q_field = sim1.sol.u[i].x[1]

        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        Πs = similar(Πr)
        Ψs = similar(Ψr)
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        Πtot = Πr .+ Πs
        Ψtot = Ψr .+ Ψs

        # Use cell-average total fields in the masked control volume.
        Htot = @. 0.5 * (Πtot^2 + Ψtot^2) - 0.5 * (Πs^2 + Ψs^2)
        E_field_tot[i] = h * sum(Htot)
        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, N, L, R, slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, N, L, R, slope_limiter)
        left_face_x = x[left_idx] - 0.5 * h
        right_face_x = x[right_idx] + 0.5 * h
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x, z, v, q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x, z, v,
                                                                       q,
                                                                       -one(v))
        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face

        # Flux through the masked control-volume faces, evaluated at the faces.
        outward_flux[i] = -Πtot_right_face * Ψtot_right_face +
                          Πtot_left_face * Ψtot_left_face

        γ = inv(sqrt(1.0 - v * v))
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        Φsz = ODE.gather_muscl(ScalarField.cell_averages_Φs(params.x, h, z, v, q), z,
                               params)
        Φtot_worldline = Φrz + Φsz

        # In the total-field bookkeeping, the particle Hamiltonian is
        # γ (m0 - q Φ_tot(z)) + V(z).
        if first(x_filtered) <= z <= last(x_filtered)
            particle_inside_mask[i] = true
            particle_energy[i] = γ * (m0 - q * Φtot_worldline) + V(z)
        else
            particle_energy[i] = 0.0
        end
        E_total[i] = E_field_tot[i] + particle_energy[i]
    end

    E_interp = PCHIPInterpolation(E_total, ts)
    dE_dt = [DataInterpolations.derivative(E_interp, t) for t in ts]
    rhs_rate = -outward_flux
    rate_residual = dE_dt .- rhs_rate

    # `cumulative_trapezoidal(ts, y)` uses the actual solver time grid `ts`.
    cumulative_rhs = cumulative_trapezoidal(ts, rhs_rate)
    delta_E = E_total .- E_total[1]
    int_residual = delta_E .- cumulative_rhs

    println("max |rate_residual| = ", maximum(abs.(rate_residual)))
    println("mean(rate_residual) = ", sum(rate_residual) / length(rate_residual))
    println("max |int_residual| = ", maximum(abs.(int_residual)))
    println("particle inside masked control volume for ",
            count(particle_inside_mask), " / ", nt, " timesteps")
    println("field energy range = [", minimum(E_field_tot), ", ", maximum(E_field_tot), "]")
    println("particle energy range = [", minimum(particle_energy), ", ",
            maximum(particle_energy), "]")

    fig = Figure(; size=(253, 520), figure_padding=(1, 5, 5, 1))

    ax_rate = Axis(fig[1, 1];
                   title=L"\textrm{Rate Balance}",
                   ylabel=L"\mathrm{Rate}")
    ax_int = Axis(fig[3, 1];
                  title=L"\textrm{Integral Balance}",
                  xlabel=L"t",
                  ylabel=L"\mathrm{Energy}")

    palette = mytheme.palette
    color1 = palette.color[][1]
    color2 = palette.color[][2]
    color4 = palette.color[][4]

    lines!(ax_rate, ts, dE_dt;
           label=L"\frac{dE_{\mathrm{tot}}}{dt}",
           color=color1, linewidth=1.5)
    lines!(ax_rate, ts, rhs_rate;
           label=L"-\mathcal{F}_{\partial\Omega}",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_rate, ts, rate_residual;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    lines!(ax_int, ts, delta_E;
           label=L"\Delta E_{\mathrm{tot}}",
           color=color1, linewidth=1.5)
    lines!(ax_int, ts, cumulative_rhs;
           label=L"-\int_0^t \mathcal{F}_{\partial\Omega}\,dt",
           color=color2, linestyle=:dash, linewidth=1.5)
    lines!(ax_int, ts, int_residual;
           label=L"\mathrm{Residual}",
           color=color4, linestyle=:dot, linewidth=1.2)

    apply_scientific_notation_if_needed!(ax_rate)
    apply_scientific_notation_if_needed!(ax_int)

    Legend(fig[2, 1], ax_rate; orientation=:horizontal, tellwidth=false, tellheight=true)
    Legend(fig[4, 1], ax_int; orientation=:horizontal, tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.09))
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Relative(0.09))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("total_field_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

end # end of module
