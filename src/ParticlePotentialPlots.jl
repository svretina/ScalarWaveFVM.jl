module ParticlePotentialPlots

using CairoMakie: regularly_spaced_array_to_range
using LinearAlgebra: sqrt_diag!
using DataInterpolations: derivative
using ..PlottingUtils
using ..ScalarField
using ..InitialData
using CairoMakie
# using GLMakie
using LinearAlgebra
using DataInterpolations

@inline name(pre, sim, post=nothing) = name_potential(pre, sim, post)
@inline name_analytic(pre, sim, post=nothing) = name_analytic_potential(pre, sim, post)

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

function plot_total_field(i, sim1)
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

    F = 2
    L = sim1.params.L
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 3] .+ ScalarField.Φs.(x_filtered, q, xp, vp)

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
    lines!(ax, x_filtered, y1_filtered; label=L"\Pi")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # savename = name("particle_potential/psifield", sim1, "i=$(i)")
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_total_field_diff(ti, sim1)
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

    # savename = name("particle_potential/psifield", sim1, "i=$(i)")
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_total_field_diff_res(ti, xi, sim1, sim2, sim4)
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
                   ScalarField.Φs.(x1_filtered, q, xp10, vp10)
    Φ20_filtered = sim2.sol.u[ti0].x[1][mask2, 3] .+
                   ScalarField.Φs.(x2_filtered, q, xp20, vp20)
    Φ40_filtered = sim4.sol.u[ti0].x[1][mask4, 3] .+
                   ScalarField.Φs.(x4_filtered, q, xp40, vp40)

    ## Fields at t=ti of simulation
    Φ1_filtered = sim1.sol.u[ti].x[1][mask1, 3] .+ ScalarField.Φs.(x1_filtered, q, xp1, vp1)
    Φ2_filtered = sim2.sol.u[2ti - 1].x[1][mask2, 3] .+
                  ScalarField.Φs.(x2_filtered, q, xp2, vp2)
    Φ4_filtered = sim4.sol.u[4ti - 3].x[1][mask4, 3] .+
                  ScalarField.Φs.(x4_filtered, q, xp4, vp4)

    ## Π fields
    Π10_filtered = sim1.sol.u[ti0].x[1][mask1, 1] .+
                   ScalarField.Πs.(x1_filtered, q, xp10, vp10)
    Π20_filtered = sim2.sol.u[ti0].x[1][mask2, 1] .+
                   ScalarField.Πs.(x2_filtered, q, xp20, vp20)
    Π40_filtered = sim4.sol.u[ti0].x[1][mask4, 1] .+
                   ScalarField.Πs.(x4_filtered, q, xp40, vp40)

    ## Fields at t=ti of simulation
    Π1_filtered = sim1.sol.u[ti].x[1][mask1, 1] .+ ScalarField.Πs.(x1_filtered, q, xp1, vp1)
    Π2_filtered = sim2.sol.u[2ti - 1].x[1][mask2, 1] .+
                  ScalarField.Πs.(x2_filtered, q, xp2, vp2)
    Π4_filtered = sim4.sol.u[4ti - 3].x[1][mask4, 1] .+
                  ScalarField.Πs.(x4_filtered, q, xp4, vp4)

    ## Ψ fields
    Ψ10_filtered = sim1.sol.u[ti0].x[1][mask1, 2] .+
                   ScalarField.Ψs.(x1_filtered, q, xp10, vp10)
    Ψ20_filtered = sim2.sol.u[ti0].x[1][mask2, 2] .+
                   ScalarField.Ψs.(x2_filtered, q, xp20, vp20)
    Ψ40_filtered = sim4.sol.u[ti0].x[1][mask4, 2] .+
                   ScalarField.Ψs.(x4_filtered, q, xp40, vp40)

    ## Fields at t=ti of simulation
    Ψ1_filtered = sim1.sol.u[ti].x[1][mask1, 2] .+ ScalarField.Ψs.(x1_filtered, q, xp1, vp1)
    Ψ2_filtered = sim2.sol.u[2ti - 1].x[1][mask2, 2] .+
                  ScalarField.Ψs.(x2_filtered, q, xp2, vp2)
    Ψ4_filtered = sim4.sol.u[4ti - 3].x[1][mask4, 2] .+
                  ScalarField.Ψs.(x4_filtered, q, xp4, vp4)

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

    # savename = name("particle_potential/psifield", sim1, "i=$(i)")
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_phifield(i, sim1; F=200)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    t1 = sim1.sol.t[i]
    t = t1

    # Create figure and axis

    x_range = (-F, F)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i].x[1][mask1, 3]  # Matches x1_filtered length

    # Compute y-limits with padding
    # y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    # y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)

    q = sim1.params.q
    dx = sim1.params.dx
    xp = sim1.sol.u[i].x[2][2]
    vp = sim1.sol.u[i].x[2][3]
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\Phi_r")
    lines!(ax, x_filtered, ScalarField.Φs.(x_filtered, q, xp, vp); label=L"\Phi_s")
    lines!(ax, x_filtered, y1_filtered .+ ScalarField.Φs.(x_filtered, q, xp, vp);
           label=L"\Phi")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # savename = name("particle_potential/psifield", sim1, "i=$(i)")
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield(i, sim1; F=2)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    t1 = sim1.sol.t[i]
    t = t1

    # Create figure and axis

    L = sim1.params.L
    x_range = (-F, F)
    @show x_range
    # Filter sim1 data (640 points)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]  # ~320 elements if symmetric around 0
    y_filtered = sim1.sol.u[i].x[1][mask, 1]  # Matches x1_filtered length

    # Compute y-limits with padding
    # y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    # y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)
    q = sim1.params.q
    dx = sim1.params.dx
    xp = sim1.sol.u[i].x[2][2]
    vp = sim1.sol.u[i].x[2][3]
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    Πr = y_filtered
    dxΠr = deriv(Πr, dx)
    Ψr = sim1.sol.u[i].x[1][mask, 2]

    denom1 = 2(1 - vp^2)^3.5
    denom2 = 4(1 - vp^2)^3
    term1 = @. (3q * abs.(x_filtered - xp) * dxΠr * vp) / denom1
    term2 = @. -(3q * abs.(x_filtered - xp) * dxΠr * vp^3) / denom1
    term3 = @. q^2 * sign.(x_filtered - xp)^2 / denom2

    @show trapezoidal((@. term1 + term2), dx)
    #     lines!(ax, x_filtered, @. term1 + term2; label=L"\Pi_r")
    lines!(ax, x_filtered, y_filtered; label=L"\Pi_r")
    lines!(ax, x_filtered, ScalarField.Πs.(x_filtered, q, xp, vp); label=L"\Pi_s")
    lines!(ax, x_filtered, y_filtered .+ ScalarField.Πs.(x_filtered, q, xp, vp);
           label=L"\Pi")

    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # savename = name("particle_potential/psifield", sim1, "i=$(i)")
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield_resolutions(i, sim1, sim2, sim4)
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1

    # Create figure and axis

    F = 2
    L = sim1.params.L
    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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

    fig = Figure(; size=(253, 200))
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
    savename = name("particle_potential/pifield", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield(i, sim1; F=100)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    if i == -1
        i = lastindex(sim1.sol.t)
    end
    t1 = sim1.sol.t[i]
    t = t1

    # Create figure and axis

    x_range = (-F, F)
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

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
              title=L"t=%$(round(t, digits=3))")
    #,
    #limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\Psi_r")
    lines!(ax, x_filtered, ScalarField.Ψs.(x_filtered, q, xp, vp); label=L"\Psi_s")
    lines!(ax, x_filtered, y1_filtered .+ ScalarField.Ψs.(x_filtered, q, xp, vp);
           linestyle=:dash, label=L"\Psi")

    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1], ax; orientation=:horizontal,
                 tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    # savename = name("particle_potential/psifield", sim1, "i=$(i)")
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_resolutions(i, sim1, sim2, sim4)
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

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
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

    fig = Figure(; size=(253, 200))
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

    savename = name("particle_potential/psifield", sim1, "i=$(i)")
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_pifield_convord(sim1, sim2, sim4)
    set_theme!(mytheme_aps())
    nt = length(sim1.sol.t)
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

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
    fig = Figure(; ize=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Order of \Pi_\mathrm{r} field}",
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

    savename = name("particle_potential/pifield_convord", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_psifield_convord(sim1, sim2, sim4)
    set_theme!(mytheme_aps())
    nt = length(sim1.sol.t)
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

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
    fig = Figure(; ize=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Order of \Psi_\mathrm{r} field}",
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

    savename = name("particle_potential/psifield_convord", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_position(sim1; invert=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][2]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{particle 1}")
        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle 1}")
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    # savename = name("particle_potential/particle_position", sim1)
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_velocity(sim1)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][3]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1]; title=L"\textrm{Particle Velocity}")
    # Plot the data
    # lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle 1}")

    lines!(ax, sim1.sol.t, tanh.(xs1); label=L"\textrm{particle 1}")
    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    # savename = name("particle_potential/particle_position", sim1)
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_mass(sim1)
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
    m2 = zeros(nt)
    # for i in 1:nt
    #     m[i] = sim1.sol.u[i].x[2][1]
    # end

    k = sim1.params.k
    q = sim1.params.q
    xs = sim1.params.x
    dx = sim1.params.dx
    V(x) = InitialData.Harmonic(x, k)
    m0 = sim1.sol.u[1].x[2][1]
    for i in 1:nt
        z = sim1.sol.u[i].x[2][2]
        v = sim1.sol.u[i].x[2][3]
        Φr = @view sim1.sol.u[i].x[1][mask, 3]
        Πr = @view sim1.sol.u[i].x[1][mask, 1]
        Ψr = @view sim1.sol.u[i].x[1][mask, 2]
        interpolator_Ψ = PCHIPInterpolation(Ψr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Π = PCHIPInterpolation(Πr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Φ = PCHIPInterpolation(Φr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)

        Πrz = interpolator_Π(z)
        Ψrz = interpolator_Ψ(z)
        Φrz = interpolator_Φ(z)

        a2 = 1.0 - v * v
        γ = 1.0 / sqrt(a2)
        # temp = @. q * sign.(xs .- z) * (3Π * v + Ψ * (1 + 2v^2)) / (2a2^2.5)
        # m2[i] = (m0 + q * Φr) / a2^1.5 - trapezoidal(temp, dx)

        # term1 = @. 0.25q^2 * sign.(x_filtered - z)^2 * γ^3

        #term2 = @. q * sign.(xsm - z) * (3Πr * v + Ψr * (1 + 2v^2)) / (2a2^2.5)
        # term2 = @. 0.5q * sign.(x_filtered - z) * (3Πr * v + Ψr * (1 + 2v^2)) * γ^5
        # # term3 = @. -δ * (q * Φr) * γ^3 + term2 + term1
        # # m = -m0 * γ^3 + trapezoidal(term3, dx)
        # term3 = @. term2 + term1
        # m[i] = -1 * (-(m0 + q * Φrz) * γ^3 + trapezoidal(term3, dx))
        term2 = @. 0.5q * γ^5 * sign.(x_filtered - z) * (3Πr * v + Ψr * (1 + 2 * v^2))
        # term3 = @. -δ * (q * Φr) * γ^3 + term2 + term1
        # m = -m0 * γ^3 + trapezoidal(term3, dx)
        # term3 = @. term2 + 0term1
        m[i] = -(-(m0 + q * Φrz) * γ^3 + trapezoidal(term2, dx))
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    # lines!(ax, sim1.sol.t, m; label=L"\textrm{particle}")
    lines!(ax, sim1.sol.t, m; label=L"\textrm{particle mass}")
    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    # savename = name("particle_potential/particle_position", sim1)
    # save(joinpath(FIG_PATH, savename), fig)
    # save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_positions_resolutions(sim1, sim2, sim4, sim1sf, sim2sf, sim4sf;
                                             invert=true)
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    xs1l = zeros(nt)
    xs1m = zeros(nt)
    xs1h = zeros(nt)

    xs2l = zeros(nt)
    xs2m = zeros(nt)
    xs2h = zeros(nt)

    for i in 1:nt
        xs1l[i] = sim1.sol.u[i].x[2][2]
        xs2l[i] = sim1sf.sol.u[i].x[2][2]

        xs1m[i] = sim2.sol.u[2i - 1].x[2][2]
        xs2m[i] = sim2sf.sol.u[2i - 1].x[2][2]

        xs1h[i] = sim4.sol.u[4i - 3].x[2][2]
        xs2h[i] = sim4sf.sol.u[4i - 3].x[2][2]
    end

    fig = Figure(; size=(253, 200))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()
    leg2 = fig[3, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    if invert
        lines!(ax, xs1l, sim1.sol.t; color=color1, linestyle=res_styles[1])
        lines!(ax, xs2l, sim1.sol.t; color=color2, linestyle=res_styles[1])

        lines!(ax, xs1m, sim1.sol.t; color=color1, linestyle=res_styles[2])
        lines!(ax, xs2m, sim1.sol.t; color=color2, linestyle=res_styles[2])

        lines!(ax, xs1h, sim1.sol.t; color=color1, linestyle=res_styles[3])
        lines!(ax, xs2h, sim1.sol.t; color=color2, linestyle=res_styles[3])

        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1)
        lines!(ax, sim1.sol.t, xs2)
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    color_legend = Legend(leg1[1, 1],
                          [LineElement(; color=color1, linestyle=res_styles[1]),
                           LineElement(; color=color2, linestyle=res_styles[1])],
                          [L"\textrm{Particle no SF}", L"\textrm{Particle SF}"];
                          halign=:center,
                          valign=:top,
                          tellwidth=false,
                          tellheight=true,
                          orientation=:horizontal)

    # # 2. Resolution legend (Low, Med, High Res)
    res_legend = Legend(leg2[1, 1],
                        [LineElement(; color=:black, linestyle=res_styles[1]),
                         LineElement(; color=:black, linestyle=res_styles[2]),
                         LineElement(; color=:black, linestyle=res_styles[3])],
                        [L"\textrm{low res}", L"\textrm{med res}", L"\textrm{high res}"];
                        halign=:center,
                        valign=:top,
                        tellwidth=false,
                        tellheight=true, orientation=:horizontal)

    # Place legends below the axis
    rowgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 2, -10)

    savename = name("particle_potential/particle_positions_resolutions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_particle_position_resolutions(sim1, sim2, sim4; invert=true)
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
    xs1l = zeros(nt)
    xs1m = zeros(nt)
    xs1h = zeros(nt)

    for i in 1:nt
        xs1l[i] = sim1.sol.u[i].x[2][2]
        xs1m[i] = sim2.sol.u[2i - 1].x[2][2]
        xs1h[i] = sim4.sol.u[4i - 3].x[2][2]
    end

    fig = Figure(; size=(253, 200))
    plt = fig[1, 1] = GridLayout()
    leg1 = fig[2, 1] = GridLayout()

    ax = Axis(plt[1, 1]; title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    if invert
        lines!(ax, xs1l, sim1.sol.t; color=color1, linestyle=res_styles[1])
        lines!(ax, xs1m, sim1.sol.t; color=color1, linestyle=res_styles[2])
        lines!(ax, xs1h, sim1.sol.t; color=color1, linestyle=res_styles[3])

        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1)
        lines!(ax, sim1.sol.t, xs2)
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    # # 2. Resolution legend (Low, Med, High Res)
    res_legend = Legend(leg1[1, 1],
                        [LineElement(; color=:black, linestyle=res_styles[1]),
                         LineElement(; color=:black, linestyle=res_styles[2]),
                         LineElement(; color=:black, linestyle=res_styles[3])],
                        [L"\textrm{low res}", L"\textrm{med res}", L"\textrm{high res}"];
                        halign=:center,
                        valign=:top,
                        tellwidth=false,
                        tellheight=true, orientation=:horizontal)

    # Place legends below the axis
    rowgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 2, -10)

    savename = name("particle_potential/particle_position_resolutions", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    save(joinpath(PAPER_FIG_PATH, savename), fig)
    return fig
end

function plot_energy(sim1; F=2)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / F, L / F)
    @show x_range
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    particle_energy = zeros(nt)
    potential = zeros(nt)
    kinetic = zeros(nt)
    self_potential = zeros(nt)
    regular_field_energy = zeros(nt)
    singular_field_energy = zeros(nt)
    mixed_terms = zeros(nt)

    q = sim1.params.q
    sf = sim1.params.sf
    dx = sim1.params.dx
    xs = sim1.params.x

    k = sim1.params.k
    V(x) = InitialData.Harmonic(x, k)
    if sf
        m0 = sim1.sol.u[1].x[2][1]
        for i in 1:nt
            m = sim1.sol.u[i].x[2][1]
            z = sim1.sol.u[i].x[2][2]
            v = sim1.sol.u[i].x[2][3]
            Φr = @view sim1.sol.u[i].x[1][mask, 3]
            Πr = @view sim1.sol.u[i].x[1][mask, 1]
            Ψr = @view sim1.sol.u[i].x[1][mask, 2]

            interpolator_Φ = PCHIPInterpolation(Φr, x_filtered;
                                                extrapolation=ExtrapolationType.Extension)
            interpolator_Π = PCHIPInterpolation(Πr, x_filtered;
                                                extrapolation=ExtrapolationType.Extension)
            interpolator_Ψ = PCHIPInterpolation(Ψr, x_filtered;
                                                extrapolation=ExtrapolationType.Extension)

            Φrz = interpolator_Φ(z)
            Πrz = interpolator_Π(z)
            Ψrz = interpolator_Ψ(z)
            Φs = ScalarField.Φs.(x_filtered, q, z, v)
            Πs = ScalarField.Πs.(x_filtered, q, z, v)
            Ψs = ScalarField.Ψs.(x_filtered, q, z, v)
            ∂vΠs = ScalarField.∂vΠs.(x_filtered, q, z, v)
            ∂vΨs = ScalarField.∂vΨs.(x_filtered, q, z, v)

            γ = 1.0 / sqrt(1.0 - v * v)
            kinetic[i] = γ * m0
            potential[i] = -q * V(z)
            self_potential[i] = -γ * q * Φrz

            particle_energy[i] = kinetic[i] + potential[i] + self_potential[i]
            Π = Πr .+ Πs
            Ψ = Ψr .+ Ψs
            ## This is with Φr and z as a generalized coordinate:
            Hfr = @. Πr * Π + 0.5 * (-Π^2 + Ψ^2) + v * (Π * ∂vΠs - Ψ * ∂vΨs)
            regular_field_energy[i] = trapezoidal(Hfr, dx)
        end
    else
        m0 = sim1.sol.u[1].x[2][1]
        for i in 1:nt
            m = sim1.sol.u[i].x[2][1]
            z = sim1.sol.u[i].x[2][2]
            v = sim1.sol.u[i].x[2][3]
            γ = 1.0 / sqrt(1.0 - v * v)
            kinetic[i] = γ * m0
            potential[i] = -q * V(z)
            particle_energy[i] = kinetic[i] + potential[i] + self_potential[i]
        end
    end

    fig = Figure(; size=(2 * 253, 2 * 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Energy}")
    # Plot the data
    # lines!(ax, sim1.sol.t, regular_field_energy .- regular_field_energy[1];
    #        label=L"\textrm{field}")
    lines!(ax, sim1.sol.t, particle_energy;
           label=L"\textrm{particle}")

    # lines!(ax, sim1.sol.t, singular_field_energy;
    #        label=L"\textrm{with singular}")

    # lines!(ax, sim1.sol.t,
    #        regular_field_energy .- regular_field_energy[1] +
    #        particle_energy .- particle_energy[1]; linestyle=:dash,
    #        label=L"\textrm{total}")
    # lines!(ax, sim1.sol.t, kinetic .- kinetic[1]; label=L"\textrm{kinetic}")
    # lines!(ax, sim1.sol.t, potential .- potential[1]; label=L"\textrm{potential}")
    # lines!(ax, sim1.sol.t, self_potential .- self_potential[1]; label=L"\textrm{self-pot}")

    # lines!(ax, sim1.sol.t, particle_energy .- particle_energy[1]; linestyle=:dash,
    #        label=L"\textrm{particle}")
    ax.ylabel = L"energy"
    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    return fig
end

function plot_work(sim1; F=2)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim1.params.L
    x_range = (-L / F, L / F)
    @show x_range
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    nt = length(sim1.sol.t)
    work = zeros(nt)
    work_potential = zeros(nt)
    dEfrdt_analytic = zeros(nt)
    dEfrdt_numerical = zeros(nt)
    dEpdt = zeros(nt)
    fluxes1 = zeros(nt)
    temp = zeros(nt)
    q = sim1.params.q
    sf = sim1.params.sf
    dx = sim1.params.dx
    xs = sim1.params.x

    function deriv(y, h)
        n = length(y)
        dy = zeros(eltype(y), n)
        for i in 2:(n - 1)
            dy[i] = (y[i - 1] - y[i + 1]) / (2h)
        end

        dy[1] = (y[1] - y[2]) / h
        dy[n] = (y[n - 1] - y[n]) / h
        return dy
    end

    k = sim1.params.k
    V(x) = InitialData.Harmonic(x, k)
    dxV(x) = InitialData.dxHarmonic(x, k)
    m0 = sim1.sol.u[1].x[2][1]
    for i in 1:nt
        m = sim1.sol.u[i].x[2][1]
        z = sim1.sol.u[i].x[2][2]
        v = sim1.sol.u[i].x[2][3]
        Φr = @view sim1.sol.u[i].x[1][mask, 3]
        Πr = @view sim1.sol.u[i].x[1][mask, 1]
        Ψr = @view sim1.sol.u[i].x[1][mask, 2]

        interpolator_Φ = PCHIPInterpolation(Φr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Π = PCHIPInterpolation(Πr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Ψ = PCHIPInterpolation(Ψr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)

        Φrz = interpolator_Φ(z)
        Πrz = interpolator_Π(z)
        Ψrz = interpolator_Ψ(z)
        Φs = ScalarField.Φs.(x_filtered, q, z, v)
        Πs = ScalarField.Πs.(x_filtered, q, z, v)
        Ψs = ScalarField.Ψs.(x_filtered, q, z, v)
        ∂vΠs = ScalarField.∂vΠs.(x_filtered, q, z, v)
        ∂vΨs = ScalarField.∂vΨs.(x_filtered, q, z, v)
        dxΠr = deriv(Πr, dx)
        dxΨr = deriv(Ψr, dx)

        a2 = 1.0 - v * v
        γ = 1.0 / sqrt(a2)

        kinetic[i] = γ * m0
        potential[i] = q * V(z)
        self_potential[i] = -γ * q * Φrz

        particle_energy[i] = kinetic[i] + potential[i] + self_potential[i]
        Π = Πr .+ Πs
        Ψ = Ψr .+ Ψs
        ## This is with Φr and z as a generalized coordinate:
        Hfr = @. Πr * Π + 0.5 * (-Π^2 + Ψ^2) + v * (Π * ∂vΠs - Ψ * ∂vΨs)
        regular_field_energy[i] = trapezoidal(Hfr, dx)
    end

    fig = Figure(; size=(2 * 253, 2 * 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Work}")
    # Plot the data
    lines!(ax, sim1.sol.t, work;
           label=L"\textrm{work}")
    # lines!(ax, sim1.sol.t, work_potential;
    #        label=L"\textrm{work of potential}")

    # lines!(ax, sim1.sol.t, work_potential .- dEfrdt;
    #        label=L"\textrm{work of potential}")

    # lines!(ax, sim1.sol.t, dEfrdt;
    #        label=L"\textrm{dE_fr/dt}")
    # lines!(ax, sim1.sol.t, dEfsdt;
    #        label=L"\textrm{dE_fs/dt}")
    # lines!(ax, sim1.sol.t, dEpdt; linestyle=:dash,
    #        label=L"\textrm{dE_p/dt}")
    # lines!(ax, sim1.sol.t, fluxes1;
    #        label=L"\textrm{flux}")
    # lines!(ax, sim1.sol.t, temp;
    #        label=L"\textrm{dEp/dt}")

    # lines!(ax, sim1.sol.t, fluxes1 .+ temp;
    #        label=L"\textrm{total}", linestyle=:dash)
    ax.ylabel = L"work"
    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    return fig
end
using RegularizationTools

function plot_balance_law(sim1; F=2)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    # F = 6
    L = sim1.params.L
    Nosc = sim1.params.Nosc
    k = sim1.params.k
    T = 2π / sqrt(abs(k))
    tf = Nosc * T
    L0 = 2tf
    x_range2 = (-L0 / 2, L0 / 2)
    x_range = (-L0 / 2, L0 / 2) #(-100, 100)
    @show x_range
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    mask2 = (sim1.params.x .>= x_range2[1]) .& (sim1.params.x .<= x_range2[2])

    x_filtered = sim1.params.x[mask]
    x_filtered2 = sim1.params.x[mask2]

    nt = length(sim1.sol.t)
    particle_energy = zeros(nt)
    potential = zeros(nt)
    kinetic = zeros(nt)
    self_potential = zeros(nt)
    regular_field_energy = zeros(nt)
    # total_energyr = zeros(nt)
    # total_energys = zeros(nt)
    flux = zeros(nt)
    flux_analytic = zeros(nt)
    bl = zeros(nt)
    dxTi0 = zeros(nt, length(x_filtered))
    Hfr = zeros(nt, length(x_filtered))
    Hfr2 = zeros(nt, length(x_filtered))
    dHpdt_analytic = zeros(nt)
    dHfrdt_analytic2 = zeros(nt, length(x_filtered))
    tmp = zeros(nt, length(x_filtered))
    # fluxes2 = zeros(nt)
    # fluxes3 = zeros(nt)
    # fluxes4 = zeros(nt)

    dHfrdt_analytic = zeros(nt)

    q = sim1.params.q
    sf = sim1.params.sf
    dx = sim1.params.dx
    xs = sim1.params.x
    dt = sim1.params.dt

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
    # V(x) = ScalarField.Φs(x, -q, xpot, 0.0)

    V(x) = InitialData.Harmonic(x, k)
    dxV(x) = InitialData.dxHarmonic(x, k)
    m0 = sim1.sol.u[1].x[2][1]

    function delta(x, z, h)
        x_left = x - h / 2
        x_right = x + h / 2
        return (z >= x_left && z < x_right) ? 1.0 : 0.0
    end
    for i in 1:nt
        m = sim1.sol.u[i].x[2][1]
        z = sim1.sol.u[i].x[2][2]
        v = sim1.sol.u[i].x[2][3]
        δ = delta.(x_filtered, z, dx)
        ## This is with Φr and z as a generalized coordinate:
        Φr = @view sim1.sol.u[i].x[1][mask, 3]
        Πr = @view sim1.sol.u[i].x[1][mask, 1]
        Ψr = @view sim1.sol.u[i].x[1][mask, 2]

        interpolator_Φ = AkimaInterpolation(Φr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Π = AkimaInterpolation(Πr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)
        interpolator_Ψ = AkimaInterpolation(Ψr, x_filtered;
                                            extrapolation=ExtrapolationType.Extension)

        Φrz = interpolator_Φ(z)
        Πrz = interpolator_Π(z)
        Ψrz = interpolator_Ψ(z)
        # Φrz = Φr .* δ
        # Πrz = Πr .* δ
        # Ψrz = Ψr .* δ

        Φs = ScalarField.Φs.(x_filtered, z, v, q)
        Πs = ScalarField.Πs.(x_filtered, z, v, q)
        Ψs = ScalarField.Ψs.(x_filtered, z, v, q)

        ∂vΦs = ScalarField.Φs001.(x_filtered, z, v, q)
        ∂vΠs = ScalarField.Πs001.(x_filtered, z, v, q)
        ∂vΨs = ScalarField.Ψs001.(x_filtered, z, v, q)
        ∂zΠs = ScalarField.Πs010.(x_filtered, z, v, q)
        ∂xΨs = ScalarField.Ψs100.(x_filtered, z, v, q)

        dxΠr = deriv(Πr, dx)
        dxΨr = deriv(Ψr, dx)

        a2 = 1.0 - v * v
        γ = 1.0 / sqrt(a2)
        signxz = sign.(x_filtered .- z)

        δ2 = delta.(x_filtered2, z, dx)
        ## This is with Φr and z as a generalized coordinate:
        Φr2 = @view sim1.sol.u[i].x[1][mask2, 3]
        Πr2 = @view sim1.sol.u[i].x[1][mask2, 1]
        Ψr2 = @view sim1.sol.u[i].x[1][mask2, 2]

        interpolator_Φ2 = ConstantInterpolation(Φr2, x_filtered;
                                                extrapolation=ExtrapolationType.Extension)
        interpolator_Π2 = ConstantInterpolation(Πr2, x_filtered;
                                                extrapolation=ExtrapolationType.Extension)
        interpolator_Ψ2 = ConstantInterpolation(Ψr2, x_filtered;
                                                extrapolation=ExtrapolationType.Extension)

        Φrz2 = interpolator_Φ2(z)
        Πrz2 = interpolator_Π2(z)
        Ψrz2 = interpolator_Ψ2(z)
        # Φrz = Φr .* δ
        # Πrz = Πr .* δ
        # Ψrz = Ψr .* δ

        Φs2 = ScalarField.Φs.(x_filtered2, z, v, q)
        Πs2 = ScalarField.Πs.(x_filtered2, z, v, q)
        Ψs2 = ScalarField.Ψs.(x_filtered2, z, v, q)

        ∂vΦs2 = ScalarField.Φs001.(x_filtered2, z, v, q)
        ∂vΠs2 = ScalarField.Πs001.(x_filtered2, z, v, q)
        ∂vΨs2 = ScalarField.Ψs001.(x_filtered2, z, v, q)
        ∂zΠs2 = ScalarField.Πs010.(x_filtered2, z, v, q)
        ∂xΨs2 = ScalarField.Ψs100.(x_filtered2, z, v, q)

        dxΠr2 = deriv(Πr2, dx)
        dxΨr2 = deriv(Ψr2, dx)

        signxz2 = sign.(x_filtered2 .- z)

        ### calculate Acceleration
        term1 = @. -q^2 * signxz2^2 / (4(-1 + v^2)^3)
        term2 = @. q * signxz2 * (3Πr2 * v + Ψr2 * (1 + 2v^2)) / (2a2^2.5)
        # term3 = @. -δ2 * (q * Φr2) * γ^3 + term2 + term1
        # m = -m0 * γ^3 + trapezoidal(term3, dx)
        term3 = @. term2 + 0term1
        m = -(m0 + q * Φrz2) * γ^3 + trapezoidal(term3, dx)

        temp = @. -0.5q * signxz2 * (v * dxΠr2 + dxΨr2) * γ^3
        # numerator = @. q * δ2 * (Ψr2 * v^2 - dxV.(x_filtered2) * a2^1.5 + Πr2 * v) * γ^3 +
        #                temp
        # acc = trapezoidal(numerator, dx) / m
        numerator = q * (Ψrz2 * v^2 + Πrz2 * v) * γ^3 + q * dxV(z) +
                    trapezoidal(temp, dx)
        acc = numerator / m
        ####

        kinetic[i] = γ * m0
        potential[i] = q * V(z) # trapezoidal(q * V.(x_filtered) .* δ, dx)
        self_potential[i] = γ * q * Φrz # trapezoidal(γ * q * Φr .* δ, dx)
        particle_energy[i] = kinetic[i] + potential[i] + self_potential[i]
        Π = Πr .+ Πs
        Ψ = Ψr .+ Ψs
        Φ = Φr .+ Φs

        tmp = @. 0.5((Πr + Πs)^2 + Ψ^2) +
                 v * (Π * ∂vΠs -
                      Ψ * ∂vΨs)
        regular_field_energy[i] = γ * (m0 + q * Φrz) + q * V(z) + trapezoidal(tmp, dx)
        ## This is with Φr and z as a generalized coordinate:
        # Hfr[i, :] .= @. Πr * Π + 0.5 * (-Π^2 + Ψ^2) + v * (Π * ∂vΠs - Ψ * ∂vΨs)
        # Hfr[i, :] .= @. 0.5(Πr^2 - 2Πr * v * ∂vΠs + Ψr * Ψr + 2Ψr * Ψs - Ψr * 2v * ∂vΨs)
        # term1 = @. -q * signxz * (Πr * v - Ψr * (1 - 2v^2)) / (2a2^1.5)
        # term1 = @. -0.5q * (Πr * v + Ψr * (-1 + 2v^2)) * γ^3

        # ##################################################################################
        # term2 = -q^2 * (-1 + 2v^2) / (8 * (-1 + v^2)) ## this term is the missing piece
        # # term2 = -0.5q^2 * v^2  ## this term is the missing piece
        # # D[(-1 + 2 v[t]^2)/(8 (-1 + v[t]^2)), t] // FullSimplifyterm2 = -v * acc / (4(-1 + v^2)^2)
        # ## but I dont remember how I derived IT !!!!
        # ##
        # # @show term2
        # # Hfr[i, :] .= term1
        # Hfr[i, :] .= @. 0.5 * (Πr^2 + Ψr^2) + signxz * term1 #+ term2 # + q^2 / 8# + q * V.(x_filtered) .* δ +
        # # γ * q * Φr .* δ + δ * γ * m0 #+ term2

        # ρ = q * δ
        # L = @. -ρ * V.(x_filtered) + 0.5 * (Π^2 - Ψ^2) - δ * m0 * sqrt(a2) -
        #        ρ * Φ * sqrt(a2)
        # DLdΠr = Π
        # ## maybe m0 should get removed?
        # # DLdz = @. Π * ∂vΠs + γ * ((δ * m0 + ρ * Φ) * v - ρ * a2 * ∂vΦs) - Ψ * ∂vΨs
        # DLdz = @. Π * ∂vΠs + γ * ((δ * m0 + ρ * Φ) * v - ρ * a2 * ∂vΦs) - Ψ * ∂vΨs
        # # DLdz = @. δ * m0 * v * γ + γ * ρ * Φ * v + Π * ∂vΠs - ρ * sqrt(a2) * ∂vΦs -
        # #           Ψ * ∂vΨs -
        # #           γ^3 * V.(x_filtered) *
        # #           (-∂zΠs +
        # #            v * (-a2 * ScalarField.Πs011.(x_filtered, z, v, q) +
        # #                 ScalarField.Ψs100.(x_filtered, z, v, q)) +
        # #            a2 * ScalarField.Ψs101.(x_filtered, z, v, q)) -
        # #           γ^2 * Φ *
        # #           (-ScalarField.Πs010.(x_filtered, z, v, q) +
        # #            v * (-a2 * ScalarField.Πs011.(x_filtered, z, v, q) +
        # #                 ScalarField.Ψs100.(x_filtered, z, v, q)) +
        # #            a2 * ScalarField.Ψs101.(x_filtered, z, v, q))

        # temp1 = @. (1 / 8) * (q^2 * signxz^2)
        # Hfr[i, :] .= @. 0.5 * (Πr^2 + Ψr^2) * γ^2 + signxz * term1 +
        #                 abs.(x_filtered .- z) * δ * q^2 + temp1
        # #Hfr[i, :] .= @. DLdΠr * Πr + DLdz * v - L
        # # Hfr[i, :] .-= Hfr[i, 1]
        # # D[H,t] + D[Ti0,x]
        # #
        # bl[i] = 2q * v * dxV(z) - 0.5q * v * γ^3 * (Πrz * v + Ψrz * (-1 + 2v^2))
        # # D[H,t] only
        # # bl[i] = -0.5 * γ^3 *
        # #         (q * Φrz * v * (-1 + 2v^2) + q * Πrz * (-2 + 3v^2) +
        # #          (-1 + v^2) * (4q * dxV(z) * v * sqrt(a2))) +
        # #         0.5 * γ^3 * a2 *
        # #         trapezoidal(q .* signxz .* dxΠr .+ 2sqrt(a2) .* (Ψr .* dxΠr .+ Πr .* dxΨr),
        # #                     dx)
        # # Hfr[i, :] .= Hfr[i, :] .- Hfr[i, 1]
        # # regular_field_energy[i] = q * V(z) + m0 * γ + q * Φrz * γ +
        # #                           trapezoidal(Hfr[i, :], dx)
        # regular_field_energy[i] = m0 * γ + trapezoidal(Hfr[i, :], dx)
        # # @show regular_field_energy[i]
        # ## After substituting particle and field eq.of.motion
        # # dHfrdt_analytic2[i, :] .= @. Ψr * dxΠr + Πr * dxΨr +
        # #                              q * γ * Πr * δ +
        # #                              0.5q * γ * signxz * dxΠr

        # ## After substituting field eq.of.motion only (no particle eqs)

        # denom = 2a2^2.5
        # # dHfrdt_analytic2[i, :] .= @. (-2q * δ * dxV.(x_filtered) * v * a2^2.5 +
        # #                               +2q * δ * Φr * v * a2 * acc + q * signxz * a2 * dxΠr -
        # #                               2q * signxz * v^2 * a2 * dxΠr -
        # #                               Ψr * (2q * δ * v^5 + q * signxz * v * acc -
        # #                                     2q * v^3 * (δ - signxz * acc) - 2sqrt(a2) * dxΠr +
        # #                                     4v^2 * sqrt(a2) * dxΠr - 2v^4 * sqrt(a2) * dxΠr) -
        # #                               q * signxz * v * a2 * dxΨr +
        # #                               Πr *
        # #                               (2v^4 * sqrt(a2) * dxΨr + 2(q * δ + sqrt(a2)dxΨr) -
        # #                                v^2 *
        # #                                (2q * δ + 3q * signxz * acc + 4sqrt(a2) * dxΨr))) /
        # #                              denom +
        # #                              δ * m0 * v * acc / a2^1.5

        # ######### field part only
        # #########
        # ######### sign^2 term after substituting acceleration
        # tmp = @. q^3 * signxz * v * (v * dxΠr + dxΨr) /
        #          ((4(-1 + v^2)^2 * (2m0 * (-1 + v^2) + 2q * δ * Φr * (-1 + v^2) +
        #                             q * signxz * (Ψr + 3Πr * v + 2Ψr * v^2))))

        # # dHfrdt_analytic2[i, :] .= @. (Ψr * (2sqrt(a2) * dxΠr +
        # #                                     v * (q * δ * (-1 + 3v^2 - 2v^4) -
        # #                                          q * signxz * (1 + 2v^2) * acc -
        # #                                          4v * sqrt(a2) * dxΠr + 2v^3 * sqrt(a2) * dxΠr)) +
        # #                               q * signxz * (-1 + v^2) *
        # #                               ((-1 + 2v^2) * dxΠr + v * dxΨr) +
        # #                               Πr * (2sqrt(a2) * dxΨr +
        # #                                     v^2 *
        # #                                     (q * δ - 3q * signxz * acc - 4sqrt(a2) * dxΨr) +
        # #                                     v^4 * (-q * δ + 2sqrt(a2) * dxΨr))) / denom # -
        # # # q^2 * v * acc / (4a2^3)
        # dHfrdt_analytic2[i, :] .= @. 0.5 *
        #                              (-q * δ * Ψr * v / a2^1.5 + q * δ * Πr * v^2 / a2^1.5) +
        #                              q * δ * Ψr * v^3 / a2^1.5 + tmp + Ψr * dxΠr +
        #                              Πr * dxΨr -
        #                              0.5(q * signxz *
        #                                  (3Πr * v^2 * acc + Ψr * v * (1 + 2v^2) * acc +
        #                                   a2 * ((-1 + 2v^2) * dxΠr + v * dxΨr)))

        # ## D[Hf,t] with substitutions, but not substituting acceleration
        # ## agrees with my pure numerical results
        # ## so not subtracting stuff
        # ## if I put the sign^2 term to zero, it agrees with the expression above. and
        # ## cancels with the fluxes
        # temp = @. Ψr * dxΠr + Πr * dxΨr +
        #           signxz * (-0.5q * γ^5 *
        #                     (3Πr * v^2 * acc + Ψr * v * (1 + 2v^2) * acc +
        #                      a2 * ((-1 + 2v^2) * dxΠr + v * dxΨr))) +
        #           0signxz^2 * (q^2 * v * acc) / (4 * (-1 + v^2)^3)
        # dHfrdt_analytic[i] = 0.5q * v * γ^3 * (Πrz * v + Ψrz * (-1 + 2 * v^2)) +
        #                      trapezoidal(temp, dx)

        # dHpdt_analytic[i] = -(-q * dxV(z) * v * a2^1.5 -
        #                       v * ((m0 + q * Φrz) * acc + q * a2 * Ψrz) + q * a2 * Πrz) /
        #                     a2^1.5

        # dHfrdt_analytic[i] = trapezoidal(dHfrdt_analytic2[i, :], dx)
        Ti0 = -Πr .* Ψ
        #intTi0 = AkimaInterpolation(Ti0, x_filtered)
        #dxti0(x) = derivative(intTi0, x, 1)
        #dxTi0[i, :] .= dxti0.(x_filtered)
        dxTi0[i, :] .= deriv(Ti0, dx)

        tmp2 = @. -Ψ * dxΠr - Πr * (dxΨr + ∂xΨs)
        flux_analytic[i] = Ti0[end] - Ti0[begin] # trapezoidal(tmp2, dx)
        flux[i] = trapezoidal(dxTi0[i, :], dx)
        bl[i] = regular_field_energy[i] + flux[i]
    end

    # Ep = AkimaInterpolation(particle_energy, sim1.sol.t)
    # dEdtp(t) = derivative(Ep, t, 1)

    Er = AkimaInterpolation(regular_field_energy, sim1.sol.t)
    dEdtr(t) = derivative(Er, t, 1)

    function time_derivative(y, h, ti)
        nt, nx = size(y)
        dy = zeros(eltype(y), nx)
        if ti == -1
            ti = nt
        end

        if ti == 1
            dy[:] .= @. (y[1, :] - y[2, :]) / h
        elseif ti == nt
            dy[:] .= @. (y[nt - 1, :] - y[nt, :]) / h
        else
            dy[:] .= @. (y[ti + 1, :] - y[ti - 1, :]) / (2h)
        end
        return dy
    end
    fig = Figure(; size=(2 * 253, 2 * 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Energy}")
    # Plot the data
    # ylims!(ax, -0.0001, 0.0001)
    # lines!(ax, sim1.sol.t, dEdtp.(sim1.sol.t); label="dEpdt")
    # lines!(ax, sim1.sol.t, bl; label="fa")
    ti = 100
    # lines!(ax, x_filtered, time_derivative(Hfr, dt, ti) .+ dxTi0[ti, :];
    #        label="Balance with dHfrdt")

    # lines!(ax, x_filtered, time_derivative(Hfr, dt, ti) .+ dxTi0[ti, :];
    #        label="adf")

    # lines!(ax, x_filtered,
    #        time_derivative(Hfr, dt, ti);
    #        label="dHfrdt")
    # lines!(ax, x_filtered, dHfrdt_analytic2[ti, :]; label="analytic",
    #        linestyle=:dash)

    # lines!(ax, x_filtered, dxTi0[ti, :];
    #        label="dxTi0")
    # lines!(ax, x_filtered,
    #        time_derivative(Hfr, dt, ti) .+ dxTi0[ti, :];
    #        label="dHfrdt")
    # lines!(ax, x_filtered, dHfrdt_analytic2[ti, :] .+ dxTi0[ti, :]; label="analytic",
    #        linestyle=:dash)

    # lines!(ax, x_filtered[(end - 200):end], Hfr[ti, (end - 200):end];
    #        label="Hfr")

    # tmp = time_derivative(Hfr, dt, ti) .- dHfrdt_analytic2[ti, :]
    # @show tmp[1:5]

    # lines!(ax, x_filtered, tmp;
    #        label="num-analytic",
    #        linestyle=:dash)

    # lines!(ax, x_filtered, tmp[ti, :]; label="tmp")

    # lines!(ax, x_filtered, dHfrdt_analytic2[ti, :] .+ dxTi0[ti, :]; label="balance")
    # lines!(ax, x_filtered, dxTi0[ti, :]; label="dxTi0")
    # lines!(ax, x_filtered[1:300], Hfr[ti - 1, 1:300]; label="dHfrdt")
    # lines!(ax, x_filtered[1:300], Hfr[ti + 1, 1:300]; label="dHfrdt")
    # lines!(ax, x_filtered[1:300], Hfr2[ti - 1, 1:300]; label="dxTi0")
    # lines!(ax, x_filtered[1:300], Hfr2[ti + 2, 1:300]; label="dxTi0")
    # @show regular_field_energy[1]
    # lines!(ax, sim1.sol.t, regular_field_energy .- regular_field_energy[1];
    #        label=L"field energy")
    # lines!(ax, sim1.sol.t, particle_energy .- particle_energy[1]; label=L"particle energy")
    # lines!(ax, sim1.sol.t,
    #        regular_field_energy .+ particle_energy .- 0particle_energy[1] .-
    #        0regular_field_energy[1];
    #        label=L"Total energy")
    # function filter(y, e)
    #     n = length(y)
    #     for i in 2:(n - 1)
    #         if abs(y[i]) > e
    #             y[i] = 0.5(y[i - 1] + y[i + 1])
    #         end
    #     end
    #     if abs(y[1]) > e
    #         y[1] = 0.0
    #     end
    #     return y
    # end
    # y = deriv(regular_field_energy .+ 0particle_energy, dt)
    # A = RegularizationSmooth(filter(y, 0.0004), sim1.sol.t, 3; alg=:gcv_svd)
    # yhat = A.û

    # lines!(ax, sim1.sol.t, yhat;
    #        label=L"\textrm{total_energy FD}")
    # lines!(ax, sim1.sol.t, yhat;
    #        label=L"\textrm{total_energy FD}")

    # lines!(ax, sim1.sol.t, dEdtr.(sim1.sol.t);
    #        label=L"\textrm{total_energy int}")
    # lines!(ax, sim1.sol.t, dHfrdt_analytic;
    #        label=L"\textrm{dHfrdt analytic}")

    # lines!(ax, sim1.sol.t, dEdtp.(sim1.sol.t);
    #        label=L"\textrm{particle_energy dt}")
    # lines!(ax, sim1.sol.t, dHpdt_analytic;
    #        label=L"\textrm{particle_energy dt analytic}")
    lines!(ax, sim1.sol.t[10:end], dEdtr.(sim1.sol.t[10:end]);
           label=L"\textrm{dHdt}")
    lines!(ax, sim1.sol.t, regular_field_energy .- regular_field_energy[1];
           label=L"\textrm{H}")

    # lines!(ax, sim1.sol.t, dEdtr.(sim1.sol.t) .+ dEdtp.(sim1.sol.t);

    #        label=L"\textrm{dt total energy num}")
    # lines!(ax, sim1.sol.t, dHfrdt_analytic .+ dHpdt_analytic;
    #        label=L"\textrm{dt total energy analytic}")
    # lines!(ax, sim1.sol.t, bl; linestyle=:dash,
    #        label=L"\textrm{bl}")

    ######
    # lines!(ax, sim1.sol.t, flux_analytic; label=L"\textrm{fluxes}")
    # lines!(ax, sim1.sol.t, flux; label=L"\textrm{flux}")

    # lines!(ax, sim1.sol.t, deriv(particle_energy, dt);
    #        label=L"\textrm{particle deriv}")
    # lines!(ax, sim1.sol.t, dHfrdt_analytic;
    #        label=L"\textrm{dHdt analytic}")
    # lines!(ax, sim1.sol.t, dHfrdt_analytic + flux;
    #        label="total")
    # lines!(ax, sim1.sol.t, deriv(regular_field_energy .+ particle_energy, dt) + flux;
    #        label="total")

    # lines!(ax, sim1.sol.t, flux_analytic;
    #        label="flux analytic")

    # lines!(ax, sim1.sol.t, deriv(regular_field_energy, dt) .+ flux_analytic;
    #        label="balance")

    # lines!(ax, sim1.sol.t, fluxes1; linestyle=:dash,
    #        label=L"\textrm{T_r^{i0}}")
    # lines!(ax, sim1.sol.t, fluxes2; linestyle=:dash,
    #        label=L"\textrm{T_r^{0i}}")

    # lines!(ax, sim1.sol.t, fluxes3; linestyle=:dash,
    #        label=L"\textrm{T_s^{i0}}")
    # lines!(ax, sim1.sol.t, fluxes4; linestyle=:dash,
    #        label=L"\textrm{T_s^{0i}}")

    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    return fig
end

function plot_hamiltonian_time_derivative(sim)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    # Get time points and assert equality

    # Create figure and axis
    L = sim.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    nt = length(sim.sol.t)

    q = sim1.params.q
    sf = sim1.params.sf
    dx = sim1.params.dx
    xs = sim1.params.x

    vs = zeros(nt)
    for i in 1:nt
        vs[i] = sim.sol.u[i].x[2][3]
    end

    vint = PCHIPInterpolation(vs, sim.sol.t; extrapolation=ExtrapolationType.Extension)
    a(t) = derivative(vint, t, 1)
    # V(x) = ScalarField.Φs(x, -q, xpot, 0.0)
    k = sim.params.k
    V(x) = InitialData.Harmonic(0.0, x, k)
    for i in 1:nt
        ti = sim.sol.t[i]
        zi = sim.sol.u[i].x[2][2]
        vi = sim.sol.u[i].x[2][3]
        ai = a(t)
        # Φ = @view sim1.sol.u[i].x[1][:, 3]
        Π = @view sim.sol.u[i].x[1][:, 1]
        Ψ = @view sim.sol.u[i].x[1][:, 2]

        γ = 1.0 / sqrt(1.0 - v * v)

        denom1 = (1.0 - vi * vi)^1.5
        denom2 = (1.0 - vi * vi)^2.5
        term1 = @. (3q * sign(x - zi) * Ψ[mask] * vi * ai) / (2denom2)
        term2 = @. (-3q * sign(x - zi) * Π[mask] * vi * vi * ai) / (2denom2)
        term3 = @. (-3q * sign(x - zi) * Ψ[mask] * vi * vi * vi * ai) / denom2
        term4 = @. (Q * sign(x - zi) * Π[mask] * ai) / (2denom1)
        term5 = @. (-2q * sign(x - zi) * Ψ[mask] * vi * ai) / denom1
        term6 = @. Π[mask] * dtΠ[mask] # how can I calculate the time deriv of Π?
        term7 = @. (-q * sign(x - zi) * vi * dtΠ[mask]) / (2denom1)
        term8 = @. Ψ[mask] * dtΨ[mask]
        term9 = @. (q * sign(x - zi) * dtΨ[mask]) / (2denom1)
        term10 = @. (-q * sign(x - zi) * vi * vi * dtΨ[mask]) / denom1
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Energy}")
    # Plot the data
    lines!(ax, sim1.sol.t, regular_field_energy1;
           label=L"\textrm{no singular}")

    # lines!(ax, sim1.sol.t,
    #        regular_field_energy .- regular_field_energy[1] +
    #        particle_energy .- particle_energy[1]; linestyle=:dash,
    #        label=L"\textrm{total}")

    lines!(ax, sim1.sol.t, particle_energy .- particle_energy[1]; linestyle=:dash,
           label=L"\textrm{particle}")
    ax.ylabel = L"energy"
    ax.xlabel = L"t"

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend
    return fig
end

### Particle rhs only
###
###

function plot_particle_position_analytic(sim1; invert=true, paper=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    analytic(t) = x0 * sin(sqrt(k) * t + 0.5pi)
    # Create figure and axis
    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    # vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i].x[2][2]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{particle}")
        lines!(ax, analytic.(sim1.sol.t), sim1.sol.t; label=L"\textrm{analytic}")
        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle}")
        lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
               label=L"\textrm{analytic}")
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_position_analytic_res(sim1, sim2, sim4; invert=true, paper=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    # Create figure and axis
    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xs2 = zeros(nt)
    xs4 = zeros(nt)

    # vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i][2]
        xs2[i] = sim2.sol.u[2i - 1][2]
        xs4[i] = sim4.sol.u[4i - 3][2]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    if invert
        lines!(ax, xs1, sim1.sol.t; label=L"\textrm{particle 1}",
               color=color1, linestyle=res_styles[1])
        lines!(ax, xs2, sim1.sol.t; color=color2, linestyle=res_styles[2])
        lines!(ax, xs4, sim1.sol.t; linestyle=res_styles[3])
        # lines!(ax, analytic.(sim1.sol.t), sim1.sol.t; label=L"\textrm{analytic}")
        ax.xlabel = L"z(t)"
        ax.ylabel = L"t"
    else
        lines!(ax, sim1.sol.t, xs1; label=L"\textrm{particle 1}",
               color=color1, linestyle=res_styles[1])
        lines!(ax, sim1.sol.t, xs2; color=color1, linestyle=res_styles[2])
        lines!(ax, sim1.sol.t, xs4; color=color1, linestyle=res_styles[3])
        lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
               label=L"\textrm{analytic}")
        ax.ylabel = L"z(t)"
        ax.xlabel = L"t"
    end
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position_res", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_position_analytic_diff(sim1, sim2, sim4; invert=true, paper=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)

    err1 = zeros(nt)
    err2 = zeros(nt)
    for i in 1:nt
        err1[i] = sim1.sol.u[i].x[2][2] - sim2.sol.u[2i - 1].x[2][2]
        err2[i] = sim2.sol.u[2i - 1].x[2][2] - sim4.sol.u[4i - 3].x[2][2]
    end
    # convfac = @. (xs1 - xs2) / (xs2 - xs4)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    lines!(ax, sim1.sol.t, err1; label=L"\textrm{low-mid}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, 2^4 * err2; label=L"\textrm{2^4\cdot(mid-high)}",
           color=color2, linestyle=res_styles[2])

    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position_diff", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_position_analytic_convfac(sim1, sim2, sim4; invert=true, paper=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    nt = length(sim1.sol.t)
    xs1 = zeros(nt)
    xs2 = zeros(nt)
    xs4 = zeros(nt)
    convfac = zeros(nt)
    # vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        xs1[i] = sim1.sol.u[i][2]
        xs2[i] = sim2.sol.u[2i - 1][2]
        xs4[i] = sim4.sol.u[4i - 3][2]
        err1 = xs1[i] - xs2[i]
        err2 = xs2[i] - xs4[i]
        if err1 == err2
            convfac[i] = 1.0
        else
            convfac[i] = err1 / err2
        end
    end
    # convfac = @. (xs1 - xs2) / (xs2 - xs4)
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Position Convergence Factor}")
    # palette = mytheme.palette
    # color1 = palette.color[][1]  # Blue-ish for Particle 1
    # color2 = palette.color[][2]  # Orange-ish for Particle 2
    # res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    lines!(ax, sim1.sol.t, convfac)
    ylims!(ax, 10, 40)
    ax.xlabel = L"t"
    # leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    #rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    #rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_position_convfac", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_velocity_analytic(sim1; paper=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    vs1 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        vs1[i] = sim1.sol.u[i].x[2][3]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity}")
    # Plot the data
    lines!(ax, sim1.sol.t, vs1; label=L"\textrm{particle 1}")
    lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
           label=L"\textrm{analytic}")
    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_velocity_analytic_res(sim1, sim2, sim4; paper=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    vs1 = zeros(nt)
    vs2 = zeros(nt)
    vs4 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        vs1[i] = sim1.sol.u[i][3]
        vs2[i] = sim2.sol.u[2i - 1][3]
        vs4[i] = sim4.sol.u[4i - 3][3]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, vs1; label=L"\textrm{particle 1}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, vs2; color=color1, linestyle=res_styles[2])
    lines!(ax, sim1.sol.t, vs4; color=color1, linestyle=res_styles[3])
    lines!(ax, sim1.sol.t, analytic.(sim1.sol.t); linestyle=:dash,
           label=L"\textrm{analytic}")
    ax.ylabel = L"v(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity_res", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_velocity_analytic_diff(sim1, sim2, sim4; paper=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    err1 = zeros(nt)
    err2 = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        err1[i] = sim1.sol.u[i][3] - sim2.sol.u[2i - 1][3]
        err2[i] = sim2.sol.u[2i - 1][3] - sim4.sol.u[4i - 3][3]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, err1; label=L"\textrm{low-mid}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, 2^4 * err2; label=L"\textrm{2^4 \cdot (mid-high)}",
           color=color2, linestyle=res_styles[2])
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity_diff", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_velocity_analytic_convfac(sim1, sim2, sim4; paper=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # analytic(t) = x0 * sin(sqrt(k / m) * t + 0.5pi)
    analytic(t) = -x0 * sqrt(k / m) * sin(sqrt(k / m) * t)
    # Create figure and axis
    nt = length(sim1.sol.t)
    # xs1 = zeros(nt)
    convfac = zeros(nt)
    # ms1 = zeros(nt)
    for i in 1:nt
        err1 = sim1.sol.u[i][3] - sim2.sol.u[2i - 1][3]
        err2 = sim2.sol.u[2i - 1][3] - sim4.sol.u[4i - 3][3]
        if err1 == err2
            convfac[i] = 1.0
        else
            convfac[i] = err1 / err2
        end
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Velocity Convergence Factor}")
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot
    ylims!(ax, 0, 40)
    # Plot the data
    lines!(ax, sim1.sol.t, convfac; color=color1, linestyle=res_styles[1])
    ax.xlabel = L"t"
    #leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    # rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_velocity_convfac", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_mass_analytic(sim1; paper=true)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x0 = sim1.params.x0
    k = sim1.params.k
    m = sim1.params.m
    # Create figure and axis
    nt = length(sim1.sol.t)
    ms1 = zeros(nt)
    for i in 1:nt
        ms1[i] = sim1.sol.u[i].x[2][1]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    lines!(ax, sim1.sol.t, ms1; label=L"\textrm{particle 1}")
    # lines!(ax, sim1.sol.t, analytic_mass)
    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_mass", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end

    return fig
end

function plot_particle_mass_analytic_res(sim1, sim2, sim4; paper=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    # Create figure and axis
    nt = length(sim1.sol.t)
    ms1 = zeros(nt)
    ms2 = zeros(nt)
    ms4 = zeros(nt)
    for i in 1:nt
        ms1[i] = sim1.sol.u[i][1]
        ms2[i] = sim2.sol.u[2i - 1][1]
        ms4[i] = sim4.sol.u[4i - 3][1]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, ms1; label=L"\textrm{particle 1}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, ms2; color=color1, linestyle=res_styles[2])
    lines!(ax, sim1.sol.t, ms4; color=color1, linestyle=res_styles[3])

    ax.ylabel = L"m(t)"
    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_mass_res", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

function plot_particle_mass_analytic_diff(sim1, sim2, sim4; paper=true)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    err1 = zeros(nt)
    err2 = zeros(nt)
    for i in 1:nt
        err1[i] = sim1.sol.u[i][1] - sim2.sol.u[2i - 1][1]
        err2[i] = sim2.sol.u[2i - 1][1] - sim4.sol.u[4i - 3][1]
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Mass}")
    # Plot the data
    palette = mytheme.palette
    color1 = palette.color[][1]  # Blue-ish for Particle 1
    color2 = palette.color[][2]  # Orange-ish for Particle 2
    res_styles = [palette.linestyle[][1], palette.linestyle[][2], palette.linestyle[][3]]  # Solid, Dash, Dot

    # Plot the data
    lines!(ax, sim1.sol.t, err1; label=L"\textrm{low-mid}",
           color=color1, linestyle=res_styles[1])
    lines!(ax, sim1.sol.t, 2^4 * err2; label=L"\textrm{2^4 \cdot (mid-high)}",
           color=color2, linestyle=res_styles[2])

    ax.xlabel = L"t"
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    colsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1)) # Space for the legend (adjustable)
    rowgap!(fig.layout, 0)                 # Gap between plot and legend

    # Save the figure
    savename = name_analytic("particle_analytic_potential/particle_mass_diff", sim1)
    save(joinpath(FIG_PATH, savename), fig)
    if paper
        save(joinpath(PAPER_FIG_PATH, savename), fig)
    end
    return fig
end

end
