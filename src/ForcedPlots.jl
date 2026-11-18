module ForcedPlots

using ..PlottingUtils
import ..PlottingUtils: get_velocity
using ..ParticleMotion
using ..ScalarField
using CairoMakie
using LinearAlgebra
using DataInterpolations
using OrdinaryDiffEq
using Printf
using QuadGK: quadgk

using ..ODE
using ..Equations
using ..RetardedPotentials
using ..Run

function plain_log_ticklabels(values)
    return map(values) do value
        if value >= 1
            return string(round(value; sigdigits=4))
        end
        return string(round(value; digits=6))
    end
end

function convergence_xticks(hs; max_labels=4)
    n = length(hs)
    if n <= max_labels
        return (hs, plain_log_ticklabels(hs))
    end

    idxs = unique(round.(Int, range(1, n; length=max_labels)))
    tickvals = hs[idxs]
    return (tickvals, plain_log_ticklabels(tickvals))
end

function _snapped_spacetime_tick_values(minval, maxval; n=5)
    tick_values = collect(range(minval, maxval; length=n))
    span = maxval - minval
    iszero(span) && return tick_values

    raw_step = span / max(n - 1, 1)
    magnitude = 10.0^floor(log10(abs(raw_step)))
    candidates = magnitude .* (1.0, 2.5, 5.0, 10.0)
    step = candidates[argmin(abs.(candidates .- raw_step))]

    start = ceil(minval / step) * step
    stop = floor(maxval / step) * step
    start > stop && return tick_values

    snapped = collect(start:step:stop)
    snapped = map(snapped) do value
        abs(value) < 1e-10 ? 0.0 : value
    end
    return snapped
end

function spacetime_colorbar_ticks(minval, maxval; n=5)
    tick_values = _snapped_spacetime_tick_values(minval, maxval; n=n)
    tick_labels = [string(round(abs(val) < 1e-10 ? 0.0 : val; digits=2))
                   for val in tick_values]
    return (tick_values, tick_labels)
end

function spacetime_colorbar_ticks_trimmed(minval, maxval; n=5)
    tick_values = _snapped_spacetime_tick_values(minval, maxval; n=n)
    tick_labels = [replace(replace(@sprintf("%.3f", abs(val) < 1e-10 ? 0.0 : val),
                                   r"0+$" => ""),
                           r"\.$" => "")
                   for val in tick_values]
    return (tick_values, tick_labels)
end

function spacetime_xticklabels(values)
    return [@sprintf("%.3f", abs(value) < 1e-10 ? 0.0 : value) for value in values]
end

function spacetime_integer_xticklabels(values)
    return [string(round(Int, abs(value) < 1e-10 ? 0.0 : value)) for value in values]
end

function padded_limits(ymin, ymax; fraction=0.1)
    ymin_f = Float64(ymin)
    ymax_f = Float64(ymax)
    span = ymax_f - ymin_f
    if iszero(span)
        scale = max(max(abs(ymin_f), abs(ymax_f)), 1.0)
        padding = fraction * scale
    else
        padding = fraction * span
    end
    return (ymin_f - padding, ymax_f + padding)
end

function plot_particle_position(sim; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim.sol.t)
    # Create figure and axis

    L = sim.params.L
    c = get_velocity(sim)
    if iszero(sim.sol.t[end])
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * sim.sol.t[end], c * sim.sol.t[end])
    end
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    x_filtered = sim.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    h = sim.params.dx
    dt = h * sim.params.cfl
    pos = zeros(nt)
    for i in 1:nt
        pos[i], _, _ = ParticleMotion.oscillator(sim.sol.t[i], sim.params.z0, sim.params.A,
                                                 sim.params.ω)
    end

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"t",
              title="Particle Position",)
    # Plot the data

    lines!(ax, pos, sim.sol.t)

    # vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("particle_position")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_spacetime(sim; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x = sim.params.x
    t = sim.sol.t
    ztraj = [ParticleMotion.oscillator(ti, sim.params.z0, sim.params.A, sim.params.ω)[1]
             for ti in t]
    L = sim.params.L
    mask = (x .>= -L / 2) .& (x .<= L / 2)
    x = x[mask]
    field = Matrix{eltype(sim.sol.u[1])}(undef, length(x), length(sim.sol.u))
    maxval = -Inf
    minval = Inf

    z0 = sim.params.z0
    v0 = sim.params.v0
    A = sim.params.A
    ω = sim.params.ω
    q = sim.params.q

    Pis0 = ScalarField.Πs.(x, Ref(z0), Ref(v0), Ref(q))

    function g(x, z, lval, rval)
        if x < z
            return lval
        elseif x > z
            return rval
        end
    end

    if total
        for (i, u) in enumerate(sim.sol.u)
            z, v, _ = ParticleMotion.oscillator(t[i], z0, A, ω)
            Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))
            column = u[mask, 1] .+ Pis .- 0g.(x, z, Pis0[1], Pis0[end])
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    else
        for (i, u) in enumerate(sim.sol.u)
            column = u[mask, 1]
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    end

    figtitle = total ? L"\Pi" : L"\Pi_\mathrm{r}"

    fig = Figure(; size=(256, 150), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"t",
              title=figtitle,
              xticklabelpad=2,
              xtickformat=spacetime_integer_xticklabels,
              yticklabelpad=4,
              aspect=Float64(2last(x)) / Float64(last(t)))
    hm = heatmap!(ax, x, t, field; rasterize=5, interpolate=true, colormap=Reverse(:RdBu),
                  colorrange=(minval, maxval))
    lines!(ax, ztraj, t; color=:black, linewidth=1.4)
    Colorbar(fig[2, 1], hm; vertical=false,
             ticks=spacetime_colorbar_ticks_trimmed(minval, maxval))

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.03))
    rowgap!(fig.layout, 1, 2)
    colsize!(fig.layout, 1, Relative(0.9))
    ax.xlabelpadding = 0.5

    savename = simple_name("pifield_spacetime")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_spacetime(sim; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    x = sim.params.x
    t = sim.sol.t
    ztraj = [ParticleMotion.oscillator(ti, sim.params.z0, sim.params.A, sim.params.ω)[1]
             for ti in t]
    L = sim.params.L
    mask = (x .>= -L / 2) .& (x .<= L / 2)
    x = x[mask]
    field = Matrix{eltype(sim.sol.u[1])}(undef, length(x), length(sim.sol.u))
    maxval = -Inf
    minval = Inf
    z0 = sim.params.z0
    v0 = sim.params.v0
    A = sim.params.A
    ω = sim.params.ω
    q = sim.params.q

    Psis0 = ScalarField.Ψs.(x, Ref(z0), Ref(v0), Ref(q))

    function g(x, z, lval, rval)
        if x < z
            return lval
        elseif x > z
            return rval
        end
    end

    if total
        for (i, u) in enumerate(sim.sol.u)
            z, v, _ = ParticleMotion.oscillator(t[i], z0, A, ω)
            Psis = ScalarField.Ψs.(x, Ref(z), Ref(v), Ref(q))
            column = u[mask, 2] .+ Psis .- 0g.(x, z, Psis0[1], Psis0[end])
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    else
        for (i, u) in enumerate(sim.sol.u)
            column = u[mask, 2]
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    end

    figtitle = total ? L"\Psi" : L"\Psi_\mathrm{r}"

    fig = Figure(; size=(256, 150), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"t",
              title=figtitle,
              xticklabelpad=2,
              xtickformat=spacetime_integer_xticklabels,
              yticklabelpad=4,
              aspect=Float64(2last(x)) / Float64(last(t)))
    hm = heatmap!(ax, x, t, field; rasterize=5, interpolate=true, colormap=Reverse(:RdBu),
                  colorrange=(minval, maxval))
    lines!(ax, ztraj, t; color=:black, linewidth=1.4)
    Colorbar(fig[2, 1], hm; vertical=false,
             ticks=spacetime_colorbar_ticks_trimmed(minval, maxval))

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.03))
    rowgap!(fig.layout, 1, 2)
    colsize!(fig.layout, 1, Relative(0.9))
    ax.xlabelpadding = 0.5

    savename = simple_name("psifield_spacetime")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_spacetime_near_particle(sim; total=false, step=1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    @assert step >= 1

    x = sim.params.x
    t = sim.sol.t[1:step:end]
    ztraj = [ParticleMotion.oscillator(ti, sim.params.z0, sim.params.A, sim.params.ω)[1]
             for ti in t]

    num = 1.5
    x_range = (-num * sim.params.A, num * sim.params.A)
    mask = (x .>= 1.1 * x_range[1]) .& (x .<= 1.1 * x_range[2])
    x = x[mask]
    u_samples = sim.sol.u[1:step:end]
    field = Matrix{eltype(sim.sol.u[1])}(undef, length(x), length(u_samples))
    maxval = -Inf
    minval = Inf
    z0 = sim.params.z0
    v0 = sim.params.v0
    A = sim.params.A
    ω = sim.params.ω
    q = sim.params.q

    Pis0 = ScalarField.Πs.(x, Ref(z0), Ref(v0), Ref(q))

    function g(x, z, lval, rval)
        if x < z
            return lval
        elseif x > z
            return rval
        end
    end

    if total
        for (i, u) in enumerate(u_samples)
            z, v, _ = ParticleMotion.oscillator(t[i], z0, A, ω)
            Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))
            column = u[mask, 1] .+ Pis .- g.(x, z, Pis0[1], Pis0[end])
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    else
        for (i, u) in enumerate(u_samples)
            column = u[mask, 1]
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    end

    figtitle = total ? L"\Pi" : L"\Pi_\mathrm{r}"

    fig = Figure(; size=(256, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"t",
              title=figtitle,
              xticklabelpad=2,
              yticklabelpad=4,
              xtickformat=spacetime_xticklabels,
              aspect=1.0)
    #aspect=Float64(last(x)) / Float64(last(t)))
    hm = heatmap!(ax, x, t, field; rasterize=5, interpolate=true, colormap=Reverse(:RdBu),
                  colorrange=(minval, maxval))
    lines!(ax, ztraj, t; color=:black, linewidth=1.4)
    xlims!(ax, x_range[1], x_range[2])
    Colorbar(fig[1, 2], hm; vertical=true,
             ticks=spacetime_colorbar_ticks_trimmed(minval, maxval))

    rowsize!(fig.layout, 1, Auto())
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))
    colsize!(fig.layout, 2, Relative(0.045))
    colgap!(fig.layout, 1, 2)
    ax.xlabelpadding = 0.5

    savename = simple_name("pifield_spacetime_near_particle")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_spacetime_near_particle(sim; total=false, step=1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    @assert step >= 1

    x = sim.params.x
    t = sim.sol.t[1:step:end]
    ztraj = [ParticleMotion.oscillator(ti, sim.params.z0, sim.params.A, sim.params.ω)[1]
             for ti in t]
    x_range = (-1.5 * sim.params.A, 1.5 * sim.params.A)
    mask = (x .>= 1.1 * x_range[1]) .& (x .<= 1.1 * x_range[2])
    x = x[mask]
    u_samples = sim.sol.u[1:step:end]
    field = Matrix{eltype(sim.sol.u[1])}(undef, length(x), length(u_samples))
    maxval = -Inf
    minval = Inf
    z0 = sim.params.z0
    v0 = sim.params.v0
    A = sim.params.A
    ω = sim.params.ω
    q = sim.params.q

    Psis0 = ScalarField.Ψs.(x, Ref(z0), Ref(v0), Ref(q))

    function g(x, z, lval, rval)
        if x < z
            return lval
        elseif x > z
            return rval
        end
    end

    if total
        for (i, u) in enumerate(u_samples)
            z, v, _ = ParticleMotion.oscillator(t[i], z0, A, ω)
            Psis = ScalarField.Ψs.(x, Ref(z), Ref(v), Ref(q))
            column = u[mask, 2] .+ Psis .- g.(x, z, Psis0[1], Psis0[end])
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    else
        for (i, u) in enumerate(u_samples)
            column = u[mask, 2]
            field[:, i] = column
            maxval = max(maxval, maximum(column))
            minval = min(minval, minimum(column))
        end
    end

    figtitle = total ? L"\Psi" : L"\Psi_\mathrm{r}"

    fig = Figure(; size=(256, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"t",
              title=figtitle,
              xticklabelpad=2,
              yticklabelpad=4,
              xtickformat=spacetime_xticklabels,
              aspect=1.0)
    #aspect=Float64(last(x)) / Float64(last(t)))
    hm = heatmap!(ax, x, t, field; rasterize=5, interpolate=true, colormap=Reverse(:RdBu),
                  colorrange=(minval, maxval))
    lines!(ax, ztraj, t; color=:black, linewidth=1.4)
    xlims!(ax, x_range[1], x_range[2])
    Colorbar(fig[1, 2], hm; vertical=true,
             ticks=spacetime_colorbar_ticks_trimmed(minval, maxval))

    rowsize!(fig.layout, 1, Auto())
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))
    colsize!(fig.layout, 2, Relative(0.045))
    colgap!(fig.layout, 1, 2)
    ax.xlabelpadding = 0.5

    savename = simple_name("psifield_spacetime_near_particle")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_regular_fields(sim; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim.sol.t)
    Πrz = zeros(nt)
    Ψrz = zeros(nt)
    vΠrz = zeros(nt)

    for i in 1:nt
        z, v, _ = ParticleMotion.oscillator(sim.sol.t[i], sim.params.z0, sim.params.A,
                                            sim.params.ω)
        Πrz[i] = AkimaInterpolation(sim.sol.u[i][:, 1], sim.params.x;
                                    extrapolation=ExtrapolationType.Extension)(z)
        Ψrz[i] = AkimaInterpolation(sim.sol.u[i][:, 2], sim.params.x;
                                    extrapolation=ExtrapolationType.Extension)(z)
        vΠrz[i] = v * Πrz[i]
    end

    y_min = minimum([minimum(Ψrz), minimum(Πrz), minimum(vΠrz)])
    y_max = maximum([maximum(Ψrz), maximum(Πrz), maximum(vΠrz)])
    padding = 0.1 * (y_max - y_min + eps(Float64))

    palette = mytheme.palette
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{field value}",
              title=L"\textrm{Particle-sampled regular fields}",
              limits=(0, sim.sol.t[end], y_min - padding, y_max + padding))

    lines!(ax, sim.sol.t, Ψrz; color=palette.color[][2], linewidth=1.6, label=L"\Psi_r(z)")
    lines!(ax, sim.sol.t, Πrz; color=palette.color[][3], linestyle=:dash, linewidth=1.5,
           label=L"\Pi_r(z)")
    lines!(ax, sim.sol.t, vΠrz; color=palette.color[][4], linestyle=:dot, linewidth=1.5,
           label=L"v \Pi_r(z)")
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("forced_motion_psipi_particle")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_regularized_trace_fields(sim; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    ts = sim.sol.t
    nt = length(ts)
    ΦR = zeros(nt)
    ΠR = zeros(nt)
    ΨR = zeros(nt)

    params = sim.sol.prob.p
    x = sim.params.x
    h = sim.params.dx
    N = params.N
    slope_limiter = params.slope_limiter
    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    function reconstruct_scalar_cell(q, idx)
        qi = ODE.get_scalar_periodic(q, idx, N)
        qim1 = ODE.get_scalar_periodic(q, idx - 1, N)
        qip1 = ODE.get_scalar_periodic(q, idx + 1, N)
        σi = slope_limiter(qi - qim1, qip1 - qi)
        return qi - 0.5 * σi, qi + 0.5 * σi
    end

    for i in 1:nt
        z, _, _ = ParticleMotion.oscillator(ts[i], sim.params.z0, sim.params.A, sim.params.ω)
        Q = sim.sol.u[i]
        cell_idx = ODE.cell_index(z, x, h, N)
        cell_center = x[cell_idx]
        cell_left_face = cell_center - 0.5 * h
        frac_left = (z - cell_left_face) / h
        frac_right = one(frac_left) - frac_left

        q_particle_minus, q_particle_plus = ODE.reconstruct_cell(Q, cell_idx, N, L, R,
                                                                 params.slope_limiter)
        Πtrace = frac_right * q_particle_minus[1] + frac_left * q_particle_plus[1]
        Ψtrace = frac_right * q_particle_minus[2] + frac_left * q_particle_plus[2]

        ϕ_minus, ϕ_plus = reconstruct_scalar_cell(@view(Q[:, 3]), cell_idx)
        Φtrace = frac_right * ϕ_minus + frac_left * ϕ_plus

        ΦR[i] = 0.5 * (Φtrace + Φtrace)
        ΠR[i] = 0.5 * (Πtrace + Πtrace)
        ΨR[i] = 0.5 * (Ψtrace + Ψtrace)
    end

    y_min = minimum([minimum(ΦR), minimum(ΠR), minimum(ΨR)])
    y_max = maximum([maximum(ΦR), maximum(ΠR), maximum(ΨR)])
    padding = 0.1 * (y_max - y_min + eps(Float64))

    palette = mytheme.palette
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{field value}",
              title=L"\textrm{Particle regularized trace values}",
              limits=(ts[1], ts[end], y_min - padding, y_max + padding))

    lines!(ax, ts, ΦR; color=palette.color[][1], linewidth=1.6,
           label=L"\Phi_R(t,z)")
    lines!(ax, ts, ΠR; color=palette.color[][2], linewidth=1.6, linestyle=:dash,
           label=L"\Pi_R(t,z)")
    lines!(ax, ts, ΨR; color=palette.color[][3], linewidth=1.6, linestyle=:dot,
           label=L"\Psi_R(t,z)")
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("forced_motion_regularized_trace_fields_particle")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_exact_regular_combo(sim; paper=false, step=1)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)
    @assert step >= 1

    ts = sim.sol.t
    nt = length(ts)
    combo = zeros(nt)
    combo_sampled = fill(NaN, nt)
    q = sim.params.q
    h = sim.params.dx
    x = sim.params.x
    params = sim.sol.prob.p

    for i in 1:nt
        t = ts[i]
        z, v, _ = ParticleMotion.oscillator(t, sim.params.z0, sim.params.A, sim.params.ω)
        Ψret = Ψ_ret_exact(z, sim.params, t)
        Πret = Π_ret_exact(z, sim.params, t)
        Ψs = ScalarField.Ψs(z, z, v, q)
        Πs = ScalarField.Πs(z, z, v, q)
        combo[i] = (Ψret - Ψs) + v * (Πret - Πs)

        if mod(i - 1, step) == 0
            sampled_combo = (Ψ_ret_exact.(x, Ref(sim.params), Ref(t)) .-
                             ScalarField.Ψs.(x, Ref(z), Ref(v), Ref(q))) .+
                            v .* (Π_ret_exact.(x, Ref(sim.params), Ref(t)) .-
                                  ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q)))

            integrand(xx) = (Ψ_ret_exact(xx, sim.params, t) - ScalarField.Ψs(xx, z, v, q)) +
                            v * (Π_ret_exact(xx, sim.params, t) - ScalarField.Πs(xx, z, v, q))

            cell_idx = ODE.cell_index(z, x, h, params.N)
            for offset in -2:2
                idx = mod(cell_idx + offset - 1, params.N) + 1
                xi = x[idx]
                left_face = xi - 0.5 * h
                right_face = xi + 0.5 * h

                if idx == cell_idx
                    I1, _ = quadgk(integrand, left_face, z; rtol=1e-10, atol=1e-12)
                    I2, _ = quadgk(integrand, z, right_face; rtol=1e-10, atol=1e-12)
                    sampled_combo[idx] = (I1 + I2) / h
                else
                    I, _ = quadgk(integrand, left_face, right_face; rtol=1e-10, atol=1e-12)
                    sampled_combo[idx] = I / h
                end
            end

            combo_sampled[i] = ODE.gather_muscl(sampled_combo, z, params)
        end
    end

    valid_sampled = filter(!isnan, combo_sampled)
    y_min = isempty(valid_sampled) ? minimum(combo) :
            minimum([minimum(combo), minimum(valid_sampled)])
    y_max = isempty(valid_sampled) ? maximum(combo) :
            maximum([maximum(combo), maximum(valid_sampled)])
    padding = 0.1 * (y_max - y_min + eps(Float64))

    palette = mytheme.palette
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\Psi_\mathrm{reg}^{\mathrm{exact}}(z) + v \Pi_\mathrm{reg}^{\mathrm{exact}}(z)",
              title=L"\textrm{Exact particle-sampled regular combo}",
              limits=(ts[1], ts[end], y_min - padding, y_max + padding))

    lines!(ax, ts, combo; color=palette.color[][2], linewidth=1.6,
           label=L"(\Psi_{ret}^{\mathrm{exact}}-\Psi_s)(z) + v(\Pi_{ret}^{\mathrm{exact}}-\Pi_s)(z)")
    if !isempty(valid_sampled)
        lines!(ax, ts, combo_sampled; color=palette.color[][3], linestyle=:dash, linewidth=1.5,
               label=L"\mathrm{MUSCL}\,\mathrm{gather}\left[(\Psi_{ret}^{\mathrm{exact}}-\Psi_s) + v(\Pi_{ret}^{\mathrm{exact}}-\Pi_s)\right]")
    end
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.16))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("forced_motion_exact_regular_combo_particle")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_mass(sim; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    ts = sim.sol.t
    nt = length(ts)
    q = sim.params.q
    m0 = sim.params.m

    meff = zeros(nt)
    mfield = zeros(nt)
    mref = zeros(nt)
    dmdt_from_mass = zeros(nt)
    dmdt_from_field = zeros(nt)
    dmdt_ref = zeros(nt)

    for i in 1:nt
        z, v, _ = ParticleMotion.oscillator(ts[i], sim.params.z0, sim.params.A,
                                            sim.params.ω)
        Πrz = AkimaInterpolation(sim.sol.u[i][:, 1], sim.params.x;
                                 extrapolation=ExtrapolationType.Extension)(z)
        Ψrz = AkimaInterpolation(sim.sol.u[i][:, 2], sim.params.x;
                                 extrapolation=ExtrapolationType.Extension)(z)
        Φrz = AkimaInterpolation(sim.sol.u[i][:, 3], sim.params.x;
                                 extrapolation=ExtrapolationType.Extension)(z)
        γ = inv(sqrt(1.0 - v * v))

        meff[i] = m0 - q * Φrz
        dmdt_from_field[i] = -q * (v * Ψrz + Πrz)
        dmdt_ref[i] = -0.5 * q * q / γ
    end

    function rhs_field(dU, U, p, t)
        z, v, _ = ParticleMotion.oscillator(t, sim.params.z0, sim.params.A, sim.params.ω)
        Q = sim.sol(t)
        Πrz = AkimaInterpolation(Q[:, 1], sim.params.x;
                                 extrapolation=ExtrapolationType.Extension)(z)
        Ψrz = AkimaInterpolation(Q[:, 2], sim.params.x;
                                 extrapolation=ExtrapolationType.Extension)(z)
        dU[1] = -q * (v * Ψrz + Πrz)
        return nothing
    end

    function rhs_ref(dU, U, p, t)
        _, v, _ = ParticleMotion.oscillator(t, sim.params.z0, sim.params.A, sim.params.ω)
        γ = inv(sqrt(1.0 - v * v))
        dU[1] = -0.5 * q * q / γ
        return nothing
    end

    prob_field = ODEProblem(rhs_field, [meff[1]], (ts[1], ts[end]))
    sol_field = solve(prob_field, RK4(); adaptive=false, dt=sim.params.dt, saveat=ts)
    mfield = Float64[u[1] for u in sol_field.u]

    prob_ref = ODEProblem(rhs_ref, [meff[1]], (ts[1], ts[end]))
    sol_ref = solve(prob_ref, RK4(); adaptive=false, dt=sim.params.dt, saveat=ts)
    mref = Float64[u[1] for u in sol_ref.u]

    dmdt_from_mass[1] = (meff[2] - meff[1]) / sim.params.dt
    for i in 2:(nt - 1)
        dmdt_from_mass[i] = (meff[i + 1] - meff[i - 1]) / (2sim.params.dt)
    end
    dmdt_from_mass[end] = (meff[end] - meff[end - 1]) / sim.params.dt

    mass_min = minimum([minimum(mref), minimum(mfield)])
    mass_max = maximum([maximum(mref), maximum(mfield)])
    mass_pad = 0.1 * (mass_max - mass_min + eps(Float64))

    rate_min = minimum([minimum(dmdt_from_mass), minimum(dmdt_from_field),
                        minimum(dmdt_ref)])
    rate_max = maximum([maximum(dmdt_from_mass), maximum(dmdt_from_field),
                        maximum(dmdt_ref)])
    rate_pad = 0.1 * (rate_max - rate_min + eps(Float64))

    palette = mytheme.palette
    fig = Figure(; size=(360, 253), figure_padding=(1, 5, 5, 1))
    ax_mass = Axis(fig[1, 1];
                   title=L"\textrm{Effective Mass}",
                   ylabel=L"m_\mathrm{eff}",
                   limits=(ts[1], ts[end], mass_min - mass_pad, mass_max + mass_pad))
    # ax_rate = Axis(fig[3, 1];
    #                title=L"\textrm{Effective Mass Rate}",
    #                xlabel=L"t",
    #                ylabel=L"\mathrm{d}m/\mathrm{d}t",
    #                limits=(ts[1], ts[end], rate_min - rate_pad, rate_max + rate_pad))

    #lines!(ax_mass, ts, meff; color=palette.color[][1], label=L"m_0 - q \Phi_r(z)")
    lines!(ax_mass, ts, mfield; color=palette.color[][2], linestyle=:dash,
           label=L"m_{\mathrm{field}}")
    lines!(ax_mass, ts, mref; color=palette.color[][3], linestyle=:dot,
           label=L"m_{q^2}")

    # lines!(ax_rate, ts, dmdt_from_mass; color=palette.color[][2],
    #        label=L"\partial_t(m_0 - q \Phi_r(z))")
    # lines!(ax_rate, ts, dmdt_from_field; color=palette.color[][3], linestyle=:dash,
    #        label=L"-q(v \Psi_r(z) + \Pi_r(z))")
    # lines!(ax_rate, ts, dmdt_ref; color=palette.color[][4], linestyle=:dot,
    #        label=L"-q^2/(2\gamma)")

    apply_scientific_notation_if_needed!(ax_mass)
    # apply_scientific_notation_if_needed!(ax_rate)

    Legend(fig[2, 1], ax_mass; orientation=:horizontal, tellwidth=false, tellheight=true)
    # Legend(fig[4, 1], ax_rate; orientation=:horizontal, nbanks=2,
    #        tellwidth=false, tellheight=true)

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.1))
    # rowsize!(fig.layout, 3, Auto())
    # rowsize!(fig.layout, 4, Relative(0.1))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("effective_mass")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_particle_regular_fields_resolutions(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    nt = length(sim1.sol.t)
    Ψ1 = zeros(nt)
    Ψ2 = zeros(nt)
    Ψ4 = zeros(nt)
    Π1 = zeros(nt)
    Π2 = zeros(nt)
    Π4 = zeros(nt)

    for i in 1:nt
        z1, _, _ = ParticleMotion.oscillator(sim1.sol.t[i], sim1.params.z0, sim1.params.A,
                                             sim1.params.ω)
        z2, _, _ = ParticleMotion.oscillator(sim2.sol.t[2i - 1], sim2.params.z0,
                                             sim2.params.A,
                                             sim2.params.ω)
        z4, _, _ = ParticleMotion.oscillator(sim4.sol.t[4i - 3], sim4.params.z0,
                                             sim4.params.A,
                                             sim4.params.ω)

        Π1[i] = AkimaInterpolation(sim1.sol.u[i][:, 1], sim1.params.x;
                                   extrapolation=ExtrapolationType.Extension)(z1)
        Ψ1[i] = AkimaInterpolation(sim1.sol.u[i][:, 2], sim1.params.x;
                                   extrapolation=ExtrapolationType.Extension)(z1)
        Π2[i] = AkimaInterpolation(sim2.sol.u[2i - 1][:, 1], sim2.params.x;
                                   extrapolation=ExtrapolationType.Extension)(z2)
        Ψ2[i] = AkimaInterpolation(sim2.sol.u[2i - 1][:, 2], sim2.params.x;
                                   extrapolation=ExtrapolationType.Extension)(z2)
        Π4[i] = AkimaInterpolation(sim4.sol.u[4i - 3][:, 1], sim4.params.x;
                                   extrapolation=ExtrapolationType.Extension)(z4)
        Ψ4[i] = AkimaInterpolation(sim4.sol.u[4i - 3][:, 2], sim4.params.x;
                                   extrapolation=ExtrapolationType.Extension)(z4)
    end

    y_min = minimum([minimum(Ψ1), minimum(Ψ2), minimum(Ψ4),
                     minimum(Π1), minimum(Π2), minimum(Π4)])
    y_max = maximum([maximum(Ψ1), maximum(Ψ2), maximum(Ψ4),
                     maximum(Π1), maximum(Π2), maximum(Π4)])
    padding = 0.1 * (y_max - y_min + eps(Float64))

    palette = mytheme.palette
    fig = Figure(; size=(253, 245), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{field value}",
              title=L"\textrm{Particle-sampled regular fields by resolution}",
              limits=(0, sim1.sol.t[end], y_min - padding, y_max + padding))

    lines!(ax, sim1.sol.t, Ψ1; color=palette.color[][2], linewidth=1.6,
           label=L"\textrm{low}, \Psi_r(z)")
    lines!(ax, sim1.sol.t, Ψ2; color=palette.color[][3], linewidth=1.4,
           label=L"\textrm{mid}, \Psi_r(z)")
    lines!(ax, sim1.sol.t, Ψ4; color=palette.color[][4], linewidth=1.4,
           label=L"\textrm{high}, \Psi_r(z)")
    lines!(ax, sim1.sol.t, Π1; color=palette.color[][2], linestyle=:dash, linewidth=1.6,
           label=L"\textrm{low}, \Pi_r(z)")
    lines!(ax, sim1.sol.t, Π2; color=palette.color[][3], linestyle=:dash, linewidth=1.4,
           label=L"\textrm{mid}, \Pi_r(z)")
    lines!(ax, sim1.sol.t, Π4; color=palette.color[][4], linestyle=:dash, linewidth=1.4,
           label=L"\textrm{high}, \Pi_r(z)")
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, nbanks=2, tellwidth=false,
           tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.18))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("forced_motion_psipi_particle_resolutions")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Fields
function plot_pifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    t = sim1.sol.t[i]
    c = get_velocity(sim1)
    L = sim1.params.L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    y1_filtered = sim1.sol.u[i][mask, 1]

    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("pifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
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
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),
              xticks=LinearTicks(5),
              yticks=LinearTicks(5),
              xgridvisible=true,
              ygridvisible=true)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("psifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_phifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t1 = sim1.sol.t[i]
    t = t1
    q = sim1.params.q
    # Create figure and axis
    xp, vp, ap = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A,
                                           sim1.params.ω)

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
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]
    y1_filtered = sim1.sol.u[i][mask1, 3] .- 0ScalarField.Πs.(x_filtered, xp, vp, q)

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered;
           label=L"\textrm{low}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("phifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Total field
function plot_total_pifield(i, sim1; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

    # Create figure and axis

    L = sim1.params.L
    c = get_velocity(sim1)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Pis = ScalarField.Πs.(x_filtered, xp, vp, sim1.params.q)
    xp0, vp0, _ = ParticleMotion.oscillator(0.0, sim1.params.z0, sim1.params.A,
                                            sim1.params.ω)
    Pis0 = ScalarField.Πs.(x_filtered, xp0, vp0, sim1.params.q)
    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = Pis .+ sim1.sol.u[i][mask1, 1] .- Pis0 # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = round(0.1 * (y_max - y_min); digits=4)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=3)

    @show padding
    @show y_limits
    @show y_max
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered)
    # lines!(ax, x_filtered, sim1.sol.u[i][mask1, 1])
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("pifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_total_psifield(i, sim1; F=2, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality
    t = sim1.sol.t[i]

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

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Psis = ScalarField.Ψs.(x_filtered, xp, vp, sim1.params.q)
    xp0, vp0, _ = ParticleMotion.oscillator(0.0, sim1.params.z0, sim1.params.A,
                                            sim1.params.ω)
    Psis0 = ScalarField.Ψs.(x_filtered, xp0, vp0, sim1.params.q)

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = Psis .+ sim1.sol.u[i][mask1, 2] .- Psis0 # Matches x1_filtered length

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = round(0.1 * (y_max - y_min); digits=4)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=4)

    @show padding
    @show y_limits
    @show y_min, y_max
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    # lines!(ax, x_filtered, y1_filtered)

    lines!(ax, x_filtered, y1_filtered)

    # lines!(ax, x_filtered, Psis0)
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("psifield_forced", t)
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
    t = sim1.sol.t[i]

    # Create figure and axis
    c = get_velocity(sim1)
    L = sim1.params.L
    @show L
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Phis = ScalarField.Φs.(x_filtered, Ref(xp), Ref(vp), Ref(sim1.params.q))

    xp0, vp0, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Phis0 = ScalarField.Φs.(x_filtered, Ref(xp0), Ref(vp0), Ref(sim1.params.q))

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]
    y1_filtered = Phis .+ sim1.sol.u[i][mask1, 3] .- Phis0

    # Compute y-limits with padding
    y_min = minimum(y1_filtered)
    y_max = maximum(y1_filtered)
    padding = round(0.1 * (y_max - y_min); digits=4)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=4)

    @show padding
    @show y_limits
    @show y_max

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Phi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    # Plot the data

    lines!(ax, x_filtered, y1_filtered)

    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    # leg = Legend(fig[2, 1],
    #              ax;
    #              orientation=:horizontal,
    #              tellwidth=false,
    #              tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    #rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("phifield_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Resolutions  
function plot_total_pifield_resolutions(i, sim1, sim2, sim4; F=8, paper=false)
    CairoMakie.activate!()
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
        x_range = (-L / F, L / F)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    xp, vp, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)
    Pis = ScalarField.Πs.(x_filtered, xp, vp, sim1.params.q)
    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = Pis .+ sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
    y2_filtered = Pis .+ y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
    y4_filtered = Pis .+ y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = round(0.1 * (y_max - y_min); digits=3)  # Add 10% padding
    y_limits = round.((y_min - padding, y_max + padding); digits=3)

    @show padding
    @show y_limits
    @show y_max
    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("pifield_res_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_resolutions(i, sim1, sim2, sim4; F=8, paper=false)
    CairoMakie.activate!()
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
        x_range = (-L / F, L / F)
    else
        x_range = (-c * t, c * t)
    end
    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Psi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("psifield_res_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_resolutions(i, sim1, sim2, sim4; F=8, paper=false)
    CairoMakie.activate!()
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
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
    y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
    y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

    # Compute y-limits with padding
    y_min = minimum([minimum(y1_filtered), minimum(y2_filtered), minimum(y4_filtered)])
    y_max = maximum([maximum(y1_filtered), maximum(y2_filtered), maximum(y4_filtered)])
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\Pi_\mathrm{r}",
              title=L"t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),)
    # Plot the data
    lines!(ax, x_filtered, y1_filtered; label=L"\textrm{low}")
    lines!(ax, x_filtered, y2_filtered; label=L"\textrm{mid}")
    lines!(ax, x_filtered, y4_filtered; label=L"\textrm{high}")
    vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line
    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("pifield_res_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Differences
function plot_pifield_diff(i, sim1, sim2, sim4; order=2, paper=false)
    CairoMakie.activate!()
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
    x_range = (-c * t, c * t)
    x_pad = 0.1 * (x_range[2] - x_range[1])
    x_range = (x_range[1] - x_pad, x_range[2] + x_pad)

    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    # mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # x1_filtered = sim1.params.x[mask1]  # ~320 elements if symmetric around 0
    y1_filtered = sim1.sol.u[i][mask, 1]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
    y2_filtered = y2_raw[mask]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
    y4_filtered = y4_raw[mask]   # Should reduce to ~320 elements

    diff1 = y1_filtered .- y2_filtered
    diff2 = y2_filtered .- y4_filtered
    # Compute y-limits with padding
    y_min = minimum([minimum(diff1), minimum(diff2)])
    y_max = maximum([maximum(diff1), maximum(diff2)])
    padding = 0.2 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    # xp, _, _ = ParticleMotion.oscillator(t, sim1.params.z0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\textrm{Error}",
              title=L"\textrm{Scaled Errors at } t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),
              xticks=LinearTicks(5),
              yticks=LinearTicks(5),
              xgridvisible=true,
              ygridvisible=true)

    # Plot the data
    lines!(ax, x_filtered, diff1; label=L"\textrm{low-mid}")
    lines!(ax, x_filtered, 2^order * diff2; label=L"2^{%$order}\cdot(\textrm{mid-high})",
           linestyle=:dot)

    # vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    apply_scientific_notation_if_needed!(ax)

    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,
                 colgap=2.5)

    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 10% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("pifield_diff_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_diff(i, sim1, sim2, sim4; order=2, paper=false)
    CairoMakie.activate!()
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
    x_range = (-c * t, c * t)
    x_pad = 0.1 * (x_range[2] - x_range[1])
    x_range = (x_range[1] - x_pad, x_range[2] + x_pad)

    mask = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    y1_filtered = sim1.sol.u[i][mask, 2]  # Matches x1_filtered length

    y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
    y2_filtered = y2_raw[mask]           # Should reduce to ~320 elements

    y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
    y4_filtered = y4_raw[mask]   # Should reduce to ~320 elements

    diff1 = y1_filtered .- y2_filtered
    diff2 = y2_filtered .- y4_filtered
    # Compute y-limits with padding
    y_min = minimum([minimum(diff1), minimum(diff2)])
    y_max = maximum([maximum(diff1), maximum(diff2)])
    padding = 0.3 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)

    # xp, _, _ = ParticleMotion.oscillator(t, sim1.params.x0, sim1.params.A, sim1.params.ω)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=L"\textrm{Errors}",
              title=L"\textrm{Scaled Errors at } t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]),
              xticks=LinearTicks(5),
              yticks=LinearTicks(5),
              xgridvisible=true,
              ygridvisible=true)
    # Plot the data
    lines!(ax, x_filtered, diff1; label=L"\textrm{low-mid}")
    lines!(ax, x_filtered, 2^order * diff2; label=L"2^{%$order}\cdot(\textrm{mid-high})",
           linestyle=:dash)
    # vlines!(ax, [xp]; color=:darkgray, linestyle=:dash, linewidth=0.8)  # Add vertical line

    apply_scientific_notation_if_needed!(ax)

    # Add legend at the bottom
    leg = Legend(fig[2, 1],
                 ax;
                 orientation=:horizontal,
                 tellwidth=false,
                 tellheight=true,
                 colgap=2.5)

    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure
    savename = simple_name("psifield_diff_forced", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Convergence Factor
function plot_pifield_convfac(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    # Create figure and axis
    F = 2
    L = sim1.params.L

    x_range = (-L / F, L / F)
    mask = (sim1.params.x .>= -L / F) .& (sim1.params.x .<= L / F)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # Filter sim2 data (1280 points)
    #mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    # mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 1:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements
        #
        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1) / norm(diff2)
        @show convfac[i]
    end
    # Compute y-limits with padding
    y_min = minimum(convfac)
    y_max = maximum(convfac)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Factor of \Pi_\mathrm{r} field}",
              limits=(sim1.sol.t[2], sim1.sol.t[end], 0, 4),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], convfac[2:end]; linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = simple_name("pifield_convfac_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_convfac(sim1, sim2, sim4; paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    # Filter sim2 data (1280 points)
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 1:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]

        @assert t1 == t2
        @assert t2 == t4
        t = t1

        y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements
        #
        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1) / norm(diff2)
    end
    # Compute y-limits with padding
    y_min = minimum(convfac)
    y_max = maximum(convfac)
    padding = 0.1 * (y_max - y_min)  # Add 10% padding
    y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              # ylabel=L"",
              title=L"\textrm{Convergence Factor of \Psi_\mathrm{r} field}",
              limits=(sim1.sol.t[2], sim1.sol.t[end], 1, 3),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], convfac[2:end]; linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = simple_name("psifield_convfac_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

## Convergence Order
function plot_pifield_self_convergence_order(sim1, sim2, sim4; step=1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    nt = length(sim1.sol.t)
    # Create figure and axis

    L = sim1.params.L
    x_range = (-L / 2, L / 2)
    mask = (sim1.params.x .>= -L / 2) .& (sim1.params.x .<= L / 2)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)

    # Filter sim2 data (1280 points)
    # mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    # mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:step:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2 * i - 1]
        t4 = sim4.sol.t[4 * i - 3]

        @assert t1 == t2
        @assert t2 == t4
        t = t1

        x_range = (-t, t)
        mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
        y1_filtered = sim1.sol.u[i][mask1, 1]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 1])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 1]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    # Compute y-limits with padding
    # y_min = minimum(convfac)
    # y_max = maximum(convfac)
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{Convergence Order}",
              title=L"\textrm{Self Convergence Test for \Pi_\mathrm{r}}",
              limits=(0, sim1.sol.t[end], 1.5, 2.5),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = simple_name("pifield_convord_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_self_convergence_order(sim1, sim2, sim4; step=1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    # Get time points and assert equality

    vmax = sim1.params.vmax
    A = sim1.params.A

    nt = length(sim1.sol.t)
    # Create figure and axis
    tf = sim1.sol.t[end]
    x_range = (-tf, tf)
    mask = (sim1.params.x .>= -tf) .& (sim1.params.x .<= tf)
    x_filtered = sim1.params.x[mask]

    # Filter sim1 data (640 points)

    # Filter sim2 data (1280 points)
    # mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])
    # Filter sim4 data (2560 points)
    # mask4 = (sim4.params.x .>= x_range[1]) .& (sim4.params.x .<= x_range[2])
    convfac = zeros(nt)
    for i in 2:step:nt
        t1 = sim1.sol.t[i]
        t2 = sim2.sol.t[2i - 1]
        t4 = sim4.sol.t[4i - 3]
        @assert t1 == t2
        @assert t2 == t4
        t = t1
        x_range = (-t, t)
        mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
        y1_filtered = sim1.sol.u[i][mask1, 2]  # Matches x1_filtered length

        y2_raw = mymean(sim2.sol.u[2i - 1][:, 2])  # ~640 elements
        y2_filtered = y2_raw[mask1]           # Should reduce to ~320 elements

        y4_raw = mymean(mymean(sim4.sol.u[4i - 3][:, 2]))  # ~1280 elements
        y4_filtered = y4_raw[mask1]   # Should reduce to ~320 elements

        diff1 = y1_filtered .- y2_filtered
        diff2 = y2_filtered .- y4_filtered
        convfac[i] = norm(diff1, 1) / norm(diff2, 1)
    end
    # Compute y-limits with padding
    # y_min = minimum(convfac)
    # y_max = maximum(convfac)
    # padding = 0.1 * (y_max - y_min)  # Add 10% padding
    # y_limits = (y_min - padding, y_max + padding)
    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{Convergence Order}",
              title=L"\textrm{Self Convergence Test for \Psi_\mathrm{r}}",
              limits=(0, sim1.sol.t[end], 0, 3),)
    # Plot the data
    lines!(ax, sim1.sol.t[2:end], log2.(convfac[2:end]); linewidth=1)
    # scatter!(ax, sim1.sol.t[2:end], convfac[2:end]; markersize=5)
    # Adjust layout to ensure proper spacing
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.2))  # Allocate 15% of height for legend
    # rowgap!(fig.layout, 2)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    # Save the figure

    savename = simple_name("psifield_convord_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_retarded_convergence_order(sim_coarse, sim_fine; norm_mode=1,
                                                 step=1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    @assert step >= 1
    ratio = sim_coarse.params.dx / sim_fine.params.dx
    ratio_int = round(Int, ratio)
    @assert isapprox(ratio, ratio_int; atol=1e-12, rtol=1e-12)
    @assert ratio_int >= 2

    L = sim_coarse.params.L
    q = sim_coarse.params.q
    A = sim_coarse.params.A
    ω = sim_coarse.params.ω
    z0 = sim_coarse.params.z0

    nt = length(sim_coarse.sol.t)
    convfac = fill(NaN, nt)
    c = get_velocity(sim_coarse)

    for i in 2:step:nt
        t1 = sim_coarse.sol.t[i]
        t2 = sim_fine.sol.t[ratio_int * (i - 1) + 1]

        @assert t1 == t2
        t = t1

        x_range = (-c * t, c * t)

        z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)

        mask1 = (sim_coarse.params.x .>= x_range[1]) .& (sim_coarse.params.x .<= x_range[2])
        x1_filtered = sim_coarse.params.x[mask1]
        yret1 = Π_ret_exact.(x1_filtered, Ref(sim_coarse.params), Ref(t)) .-
                ScalarField.Πs.(x1_filtered, z, v, q)
        y1_filtered = sim_coarse.sol.u[i][mask1, 1]
        errvec1 = Float64[(yret1[j] - y1_filtered[j]) for j in eachindex(y1_filtered)]
        err1 = norm(errvec1, norm_mode)

        @assert ispow2(ratio_int)
        y2_coarsened = sim_fine.sol.u[ratio_int * (i - 1) + 1][:, 1]
        for _ in 1:trailing_zeros(ratio_int)
            y2_coarsened = mymean(y2_coarsened)
        end
        y2_coarse = y2_coarsened[mask1]
        @assert length(y2_coarse) == length(y1_filtered)
        errvec2 = Float64[(yret1[j] - y2_coarse[j]) for j in eachindex(y1_filtered)]
        err2 = norm(errvec2, norm_mode)

        convfac[i] = log(err1 / err2) / log(ratio_int)
    end

    plot_indices = collect(2:step:nt)

    fig = Figure(; size=(260, 220), figure_padding=(1, 1, 1, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{Convergence order}",
              title=L"\textrm{Comparison to analytic solution for } \Pi_\mathrm{r}",
              limits=(0, sim_coarse.sol.t[end], 0, 3),)
    lines!(ax, sim_coarse.sol.t[plot_indices], convfac[plot_indices]; linewidth=1)
    rowsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("pifield_retarded_convord_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_retarded_convergence_order(sim_coarse, sim_fine; norm_mode=1,
                                                  step=1, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())
    @assert step >= 1
    nt = length(sim_coarse.sol.t)
    L = sim_coarse.params.L
    q = sim_coarse.params.q
    A = sim_coarse.params.A
    ω = sim_coarse.params.ω
    z0 = sim_coarse.params.z0
    c = get_velocity(sim_coarse)
    ratio = sim_coarse.params.dx / sim_fine.params.dx
    ratio_int = round(Int, ratio)
    @assert isapprox(ratio, ratio_int; atol=1e-12, rtol=1e-12)
    @assert ratio_int >= 2

    convfac = fill(NaN, nt)

    tf = sim_coarse.sol.t[end]
    x_range = (-c * tf, c * tf)
    #@show x_range
    #x_range = (-1, 1)
    #mask1 = (sim_coarse.params.x .>= x_range[1]) .& (sim_coarse.params.x .<= x_range[2])
    mask1 = (sim_coarse.params.x .> x_range[1]) .& (sim_coarse.params.x .< x_range[2])
    x1_filtered = sim_coarse.params.x[mask1]

    for i in 2:step:nt
        t1 = sim_coarse.sol.t[i]
        t2 = sim_fine.sol.t[ratio_int * (i - 1) + 1]
        @assert isapprox(t1, t2; atol=1e-12, rtol=1e-12)
        t = t1
        buffer = 0.025 * 8 #5 * sim_coarse.params.dx
        x_range = (-c * t + buffer, c * t - buffer)
        #x_range = (-c * t, c * t)
        mask1 = (sim_coarse.params.x .> x_range[1]) .& (sim_coarse.params.x .< x_range[2])
        x1_filtered = sim_coarse.params.x[mask1]

        z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)

        yret1 = Ψ_ret.(x1_filtered, Ref(sim_coarse.params), Ref(t)) .-
                ScalarField.Ψs.(x1_filtered, z, v, q)

        y1_filtered = sim_coarse.sol.u[i][mask1, 2]
        errvec1 = Float64.(yret1 .- y1_filtered)
        #@show errvec1
        h1 = sim_coarse.params.dx
        err1 = h1 * norm(errvec1, norm_mode)
        # i == 20 && throw(DomainError)

        @assert ispow2(ratio_int)
        y2_coarsened = sim_fine.sol.u[ratio_int * (i - 1) + 1][:, 2]

        for _ in 1:trailing_zeros(ratio_int)
            y2_coarsened = mymean(y2_coarsened)
        end
        y2_coarse = y2_coarsened[mask1]
        @assert length(y2_coarse) == length(y1_filtered)

        errvec2 = Float64[(yret1[j] - y2_coarse[j]) for j in eachindex(y1_filtered)]
        err2 = h1 * norm(errvec2, norm_mode)

        convfac[i] = log2(err1 / err2)
    end
    plot_indices = collect(2:step:nt)
    @show convfac[plot_indices]
    fig = Figure(; size=(260, 220), figure_padding=(1, 1, 1, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{Convergence order}",
              title=L"\textrm{Comparison to analytic solution for } \Psi_\mathrm{r}",
              limits=(0, sim_coarse.sol.t[end], 0, 3),)
    lines!(ax, sim_coarse.sol.t[plot_indices], convfac[plot_indices]; linewidth=1)
    rowsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("psifield_retarded_convord_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_convergence_panel(sim1, sim2, sim4; paper=false)
    nt = length(sim1.sol.t)
    fig1 = plot_pifield_resolutions(nt, sim1, sim2, sim4; paper=paper)  # Your existing function
    fig2 = plot_psifield_resolutions(nt, sim1, sim2, sim4; paper=paper)  # Replace with your second function
    fig3 = plot_pifield_self_convergence_order(sim1, sim2, sim4; paper=paper)  # Replace with your third function
    fig4 = plot_psifield_self_convergence_order(sim1, sim2, sim4; paper=paper)  # Replace with your fourth function

    fig = Figure(; size=(600, 500))

    # Add supertitle
    supertitle = Label(fig[1, 1:2],
                       L"v_{max}=%$(sim1.params.vmax)";
                       fontsize=20,
                       font=:bold,
                       halign=:center,)

    # Function to recreate Axis content
    function recreate_axis!(fig, src_fig, row, col)
        src_ax = contents(src_fig[1, 1])[1]  # Original Axis

        # Create new Axis with same properties
        new_ax = Axis(fig[row, col];
                      xlabel=src_ax.xlabel[],
                      ylabel=src_ax.ylabel[],
                      title=src_ax.title[],
                      limits=src_ax.limits[],)

        # Copy plot objects
        for plot in src_ax.scene.plots
            if plot isa MakieCore.Scatter
                label = haskey(plot.attributes, :label) ? plot.label[] : nothing
                scatter!(new_ax,
                         plot[1][],
                         plot[2][];
                         label=label,
                         markersize=plot.markersize[],)
            elseif plot isa MakieCore.Lines
                data = plot[1][]  # Single Observable
                label = haskey(plot.attributes, :label) ? plot.label[] : nothing
                if isa(data, Vector{Point{2,Float64}})  # Regular Lines
                    x = [p[1] for p in data]
                    y = [p[2] for p in data]
                    lines!(new_ax, x, y; label=label, linewidth=plot.linewidth[])
                elseif isa(data, Vector{Float64})  # VLines misclassified as Lines
                    @show "ASDFA"
                    vlines!(new_ax,
                            data;
                            color=plot.color[],
                            linestyle=plot.linestyle[],
                            linewidth=plot.linewidth[],)
                end
            end
        end

        # Add Legend if it exists, place in the next row under the column
        leg_content = contents(src_fig[2, 1])
        if !isempty(leg_content)
            Legend(fig[row + 1, col],
                   new_ax;
                   orientation=:horizontal,
                   tellwidth=false,
                   tellheight=true,)
        end
    end

    # Recreate each plot in the grid
    recreate_axis!(fig, fig1, 2, 1)  # Top-left: Π field resolutions (Axis in row 2, Legend in row 3)
    recreate_axis!(fig, fig2, 2, 2)  # Top-right: Ψ field resolutions (Axis in row 2, Legend in row 3)
    recreate_axis!(fig, fig3, 4, 1)  # Bottom-left: Π field convergence (Axis in row 4)
    recreate_axis!(fig, fig4, 4, 2)  # Bottom-right: Ψ field convergence (Axis in row 4)

    # Adjust layout with padding on the right
    rowsize!(fig.layout, 1, Relative(0.1)) # Supertitle row
    rowsize!(fig.layout, 2, Auto())       # Top row of plots (fig1, fig2)
    rowsize!(fig.layout, 3, Relative(0.05)) # Legend row for fig1 and fig2
    rowsize!(fig.layout, 4, Auto())       # Bottom row of plots (fig3, fig4)
    # rowsize!(fig.layout, 5, Relative(0.0)) # Empty row
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.5)) # Left column (50% of width)
    colsize!(fig.layout, 2, Relative(0.45)) # Right column (45% of width, leaves 5% padding)
    # Save the combined figure

    savename = simple_name("convergence_pane_forced")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

### Comparison with Retarded
function plot_pifield_vs_retarded(i, sim; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim.sol.t[i]

    h = sim.params.dx
    L = sim.params.L
    q = sim.params.q
    A = sim.params.A
    ω = sim.params.ω
    z0 = sim.params.z0
    z0 = sim.params.z0

    # c = get_velocity(sim)
    # if iszero(t)
    #     x_range = (-L / 2, L / 2)
    # else
    #     x_range = (-c * t, c * t)
    # end
    x_range = (-L / 2, L / 2)
    mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
    x = sim.params.x[mask]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))

    if total
        y_sim = Pis .+ sim.sol.u[i][mask, 1]
    else
        y_sim = sim.sol.u[i][mask, 1] #.- 0Pis0
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret_exact.(x, Ref(sim.params), Ref(t))
    else
        y_lw = Π_ret_exact.(x, Ref(sim.params), Ref(t)) .- ScalarField.Πs.(x, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim)
    ymax = maximum([ymax1, ymax2])

    y_limits = padded_limits(ymin, ymax)

    fig = Figure(; size=(253, 200))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=total ? L"\Pi" : L"\Pi_r",
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

function plot_psifield_vs_retarded(i, sim; total=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim.sol.t[i]

    h = sim.params.dx
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
        y_sim = sim.sol.u[i][mask, 2]
    end

    ## Lienard-Wiechert

    if total
        y_lw = Ψ_ret_exact.(x, Ref(sim.params), Ref(t))
    else
        y_lw = Ψ_ret_exact.(x, Ref(sim.params), Ref(t)) .- ScalarField.Ψs.(x, z, v, q)
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim)
    ymax = maximum([ymax1, ymax2])

    y_limits = padded_limits(ymin, ymax)

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

function plot_pifield_resolutions_vs_retarded(i, sim1, sim2; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim1.sol.t[i]

    h1 = sim1.params.dx
    h2 = sim2.params.dx
    n1 = round(Int, inv(h1))
    n2 = round(Int, inv(h2))
    label1 = isapprox(inv(h1), n1; atol=1e-12, rtol=1e-12) ? L"dx=1/%$n1" : L"dx=%$h1"
    label2 = isapprox(inv(h2), n2; atol=1e-12, rtol=1e-12) ? L"dx=1/%$n2" : L"dx=%$h2"
    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    c = get_velocity(sim1)
    #  x_range = (-c * t, c * t)
    x_range = (-L / 2, L / 2)
    x_pad = 0.1 * (x_range[2] - x_range[1])
    x_range = (x_range[1] - x_pad, x_range[2] + x_pad)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])

    x1 = sim1.params.x[mask1]
    x2 = sim2.params.x[mask2]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    z0, v0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
    Pis = ScalarField.Πs.(x1, Ref(z), Ref(v), Ref(q))

    function g(x, z, lval, rval)
        if x < z
            return lval
        elseif x > z
            return rval
        end
    end

    if total
        Pis2 = ScalarField.Πs.(x2, Ref(z), Ref(v), Ref(q))
        Pis10 = ScalarField.Πs.(x1, Ref(z0), Ref(v0), Ref(q))
        Pis20 = ScalarField.Πs.(x2, Ref(z0), Ref(v0), Ref(q))
        y_sim1 = sim1.sol.u[i][mask1, 1] .+ Pis .- 0g.(x1, z, Pis10[1], Pis10[end])
        y_sim2 = sim2.sol.u[2i - 1][mask2, 1] .+ Pis2 .- 0g.(x2, z, Pis20[1], Pis20[end])
    else
        # xp0, vp0, _ = ParticleMotion.oscillator(0.0, z0, A, ω)
        # Pis0 = ScalarField.Πs.(x, Ref(xp0), Ref(vp0), Ref(q))
        y_sim1 = sim1.sol.u[i][mask1, 1] #.- 0Pis0
        y_sim2 = sim2.sol.u[2i - 1][mask2, 1]
    end

    ## Lienard-Wiechert

    if total
        vals = Π_ret_exact.(x1, Ref(sim1.params), Ref(0.0))
        y_lw = Π_ret_exact.(x1, Ref(sim1.params), Ref(t)) .- 0g.(x1, z, vals[1], vals[end])
    else
        y_lw = Π_ret_exact.(x1, Ref(sim1.params), Ref(t)) .- Pis
    end
    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim1)

    ymin = minimum([ymin1, ymin2])
    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim1)
    ymax = maximum([ymax1, ymax2])
    y_limits = padded_limits(ymin, ymax)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=total ? L"\Pi" : L"\Pi_\mathrm{r}",
              title=L"\textrm{Solution at } t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # xlims!(ax, -0.02, 0.0)
    @show y_lw[(length(y_lw) ÷ 2 - 10):(length(y_lw) ÷ 2 + 10)]
    @show x1[(length(x1) ÷ 2 - 10):(length(x1) ÷ 2 + 10)]
    @show z, v

    lines!(ax, x1, y_lw; label="Analytic")
    lines!(ax, x1, y_sim1; label=label1)
    lines!(ax, x2, y_sim2; label=label2, linestyle=:dash)

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

    savename = simple_name("pifield_retarded_resolutions", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

function plot_psifield_resolutions_vs_retarded(i, sim1, sim2; total=false, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    t = sim1.sol.t[i]

    h1 = sim1.params.dx
    h2 = sim2.params.dx
    n1 = round(Int, inv(h1))
    n2 = round(Int, inv(h2))
    label1 = isapprox(inv(h1), n1; atol=1e-12, rtol=1e-12) ? L"dx=1/%$n1" : L"dx=%$h1"
    label2 = isapprox(inv(h2), n2; atol=1e-12, rtol=1e-12) ? L"dx=1/%$n2" : L"dx=%$h2"
    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    c = get_velocity(sim1)
    x_range = (-c * t, c * t)
    x_pad = 0.1 * (x_range[2] - x_range[1])
    x_range = (x_range[1] - x_pad, x_range[2] + x_pad)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    mask2 = (sim2.params.x .>= x_range[1]) .& (sim2.params.x .<= x_range[2])

    x1 = sim1.params.x[mask1]
    x2 = sim2.params.x[mask2]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Psis = ScalarField.Ψs.(x1, Ref(z), Ref(v), Ref(q))

    if total
        Psis2 = ScalarField.Ψs.(x2, Ref(z), Ref(v), Ref(q))
        y_sim1 = sim1.sol.u[i][mask1, 2] .+ Psis
        y_sim2 = sim2.sol.u[2i - 1][mask2, 2] .+ Psis2
    else
        y_sim1 = sim1.sol.u[i][mask1, 2]
        y_sim2 = sim2.sol.u[2i - 1][mask2, 2]
    end
    ## Lienard-Wiechert

    if total
        y_lw = Ψ_ret_exact.(x1, Ref(sim1.params), Ref(t))
    else
        y_lw = Ψ_ret_exact.(x1, Ref(sim1.params), Ref(t)) .- Psis
    end

    ymin1 = minimum(y_lw)
    ymin2 = minimum(y_sim1)

    ymin = minimum([ymin1, ymin2])

    ymax1 = maximum(y_lw)
    ymax2 = maximum(y_sim1)
    ymax = maximum([ymax1, ymax2])

    y_limits = padded_limits(ymin, ymax)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"x",
              ylabel=total ? L"\Psi" : L"\Psi_\mathrm{r}",
              title=L"\textrm{Solution at } t=%$(round(t, digits=3))",
              limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    lines!(ax, x1, y_lw; label="Analytic")
    lines!(ax, x1, y_sim1; label=label1)
    lines!(ax, x2, y_sim2; label=label2, linestyle=:dash)
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

    savename = simple_name("psifield_retarded_resolutions", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

function plot_pifield_diff_vs_retarded(i, sim1, sim2; total=false, order=0, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    h1 = sim1.params.dx
    h2 = sim2.params.dx

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    @assert t1 == t2
    t = t1

    c = get_velocity(sim1)
    # if iszero(t)
    #     x_range = (-L / 2, L / 2)
    # else
    #     x_range = (-c * t, c * t)
    # end
    x_range = (-1.1c * t, 1.1c * t)
    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x = sim1.params.x[mask1]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Pis = ScalarField.Πs.(x, Ref(z), Ref(v), Ref(q))
    if total
        y_sim1 = Pis .+ sim1.sol.u[i][mask1, 1]
        y_sim2 = Pis .+ mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    else
        y_sim1 = sim1.sol.u[i][mask1, 1]
        y_sim2 = mymean(sim2.sol.u[2i - 1][:, 1])[mask1]
    end

    ## Lienard-Wiechert

    if total
        y_lw = Π_ret.(x, Ref(sim1.params), Ref(t))
    else
        y_lw = Π_ret.(x, Ref(sim1.params), Ref(t)) .- Pis
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
              ylabel=L"\textrm{\Pi_r - \Pi_{ret}}",
              title=L"\textrm{Scaled Errors at t=%$(round(t, digits=3)) }")
    # ,
    # limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))
    # lines!(ax, x, convfac)
    lines!(ax, x, diff1; label=L"\textrm{low resolution} - \textrm{analytic}")
    lines!(ax, x, 2^(order) * diff2;
           label=L"2^{%$(order)} (\textrm{high resolution} - \textrm{analytic})",
           linestyle=:dash)

    # lines!(ax, x, y_sim1)
    # lines!(ax, x, y_sim2)
    # lines!(ax, x, y_lw; linestyle=:dash, label="Retarded")

    ylims!(ax, -0.0005, 0.0004)
    xlims!(ax, -2.35, -1.95)

    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false,
                 tellheight=true, nbanks=2, framevisible=false, colgap=5, rowgap=5,
                 patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    savename = simple_name("force_pifield_diff_ret", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_psifield_diff_vs_retarded(i, sim1, sim2; total=false, order=0, paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    h1 = sim1.params.dx
    h2 = sim2.params.dx

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    @assert t1 == t2
    t = t1

    c = get_velocity(sim1)
    if iszero(t)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t, c * t)
    end

    mask1 = (sim1.params.x .>= x_range[1]) .& (sim1.params.x .<= x_range[2])
    x = sim1.params.x[mask1]

    z, v, _ = ParticleMotion.oscillator(t, z0, A, ω)
    Psis = ScalarField.Ψs.(x, Ref(z), Ref(v), Ref(q))
    if total
        y_sim1 = Psis .+ sim1.sol.u[i][mask1, 2]
        y_sim2 = Psis .+ mymean(sim2.sol.u[2i - 1][:, 2])[mask1]
    else
        y_sim1 = sim1.sol.u[i][mask1, 2]
        y_sim2 = mymean(sim2.sol.u[2i - 1][:, 2])[mask1]
    end

    ## Lienard-Wiechert

    if total
        y_lw = Ψ_ret_exact.(x, Ref(sim1.params), Ref(t))
    else
        y_lw = Ψ_ret.(x, Ref(sim1.params), Ref(t)) .- Psis
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
              ylabel=L"\textrm{\Psi_r - \Psi_{ret}}",
              title=L"\textrm{Scaled Errors at t=%$(round(t, digits=3)) }")
    # ,
    # limits=(x_range[1], x_range[2], y_limits[1], y_limits[2]))

    # lines!(ax, x, convfac)
    lines!(ax, x, diff1; label=L"\textrm{low resolution} - \textrm{analytic}")
    lines!(ax, x, 2^(order) * diff2;
           label=L"2^{%$order} (\mathrm{medium} \; \mathrm{resolution} - \mathrm{analytic})",
           linestyle=:dash)

    # lines!(ax, x, y_sim1)
    # lines!(ax, x, y_sim2)
    # lines!(ax, x, y_lw; linestyle=:dash, label="Retarded")
    leg = Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false,
                 tellheight=true, nbanks=2, framevisible=false, colgap=5, rowgap=5,
                 patchsize=(15, 10))
    rowsize!(fig.layout, 1, Auto())  # Let the axis size adjust
    # rowsize!(fig.layout, 2, Relative(0.1))  # Allocate 15% of height for legend
    rowgap!(fig.layout, 0)  # Small gap between plot and legend
    colsize!(fig.layout, 1, Relative(0.9))
    savename = simple_name("force_psifield_diff_ret", t)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_pifield_error_convergence(t, sims...; exponent=2.0, norm_mode=Inf,
                                        total=false,
                                        paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    @assert !isempty(sims)
    sim1 = first(sims)
    i1 = argmin(abs.(sim1.sol.t .- t))
    t1 = sim1.sol.t[i1]

    indices = Int[i1]
    for sim in Iterators.drop(sims, 1)
        ratio = sim1.params.dx / sim.params.dx
        ratio_int = round(Int, ratio)
        @assert isapprox(ratio, ratio_int; atol=1e-12, rtol=1e-12)
        idx = ratio_int * (i1 - 1) + 1
        @assert 1 <= idx <= length(sim.sol.t)
        @assert isapprox(sim.sol.t[idx], t1; atol=1e-12, rtol=1e-12)
        push!(indices, idx)
    end

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    hs = [sim.params.dx for sim in sims]
    errors = zeros(length(sims))

    c = get_velocity(sim1)
    if iszero(t1)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t1, c * t1)
    end

    z, v, _ = ParticleMotion.oscillator(t1, z0, A, ω)
    for (k, (sim, idx)) in enumerate(zip(sims, indices))
        mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
        x = sim.params.x[mask]
        y_ret = Π_ret.(x, Ref(sim.params), Ref(t1)) .- ScalarField.Πs.(x, z, v, q)
        y_sim = sim.sol.u[idx][mask, 1]
        ratio = sim1.params.dx / sim.params.dx
        ratio_int = round(Int, ratio)
        @assert isapprox(ratio, ratio_int; atol=1e-12, rtol=1e-12)
        errors[k] = norm(abs.(y_ret .- y_sim)[begin:ratio_int:end], norm_mode)
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    xticks = convergence_xticks(hs)
    if norm_mode == 2
        title_str = L"\Vert \Pi_{\textrm{ret}}-\Pi_\textrm{h} \Vert_2"
    elseif norm_mode == 1
        title_str = L"\Vert \Pi_{\textrm{ret}}-\Pi_\textrm{h} \Vert_1"
    else
        title_str = L"\Vert \Pi_{\textrm{ret}}-\Pi_\textrm{h} \Vert_{\infty}"
    end

    ax = Axis(fig[1, 1];
              xlabel=L"dx",
              yscale=log10,
              xscale=log10,
              xticks=xticks,
              ytickformat=plain_log_ticklabels,
              xreversed=true,
              title=title_str)

    lines!(ax, hs, errors)
    scatter!(ax, hs, errors)

    nref = length(hs) #min(3, length(hs))
    xref = hs[1:nref]
    yref1 = xref .^ exponent
    shift1 = (errors[1] / yref1[1]) * 2.5
    yref1 .*= shift1

    lines!(ax, xref, yref1; linestyle=:dash, label=L"h^{%$exponent}")

    Legend(fig[1, 1], ax;
           halign=:right,
           valign=:top,
           tellwidth=false,
           tellheight=false,
           margin=(10, 10, 10, 10),
           framevisible=false,
           labelsize=12)

    ax.ylabel = L"\mathrm{Error}"

    rowsize!(fig.layout, 1, Auto())
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("pifield_retarded_convergence", t1)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

function plot_pifield_error_convergence(i, Nosc, dx0, cfl, vmax, q;
                                        exponent=2.0,
                                        norm_mode=Inf,
                                        total=false,
                                        paper=false)
    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q)
    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q)
    sim8 = Run.forced_motion(Nosc, dx0 / 8, cfl, vmax, q)
    sim16 = Run.forced_motion(Nosc, dx0 / 16, cfl, vmax, q)
    sim32 = Run.forced_motion(Nosc, dx0 / 32, cfl, vmax, q)
    t = sim1.sol.t[i]
    return plot_pifield_error_convergence(t, sim1, sim2, sim4, sim8, sim16, sim32;
                                          exponent=exponent, norm_mode=norm_mode,
                                          total=total, paper=paper)
end

function plot_psifield_error_convergence(t, sims...; exponent=1.0, norm_mode=Inf,
                                         total=false,
                                         paper=false)
    CairoMakie.activate!()
    set_theme!(mytheme_aps())

    @assert !isempty(sims)
    sim1 = first(sims)
    i1 = argmin(abs.(sim1.sol.t .- t))
    t1 = sim1.sol.t[i1]

    indices = Int[i1]
    for sim in Iterators.drop(sims, 1)
        ratio = sim1.params.dx / sim.params.dx
        ratio_int = round(Int, ratio)
        @assert isapprox(ratio, ratio_int; atol=1e-12, rtol=1e-12)
        idx = ratio_int * (i1 - 1) + 1
        @assert 1 <= idx <= length(sim.sol.t)
        @assert isapprox(sim.sol.t[idx], t1; atol=1e-12, rtol=1e-12)
        push!(indices, idx)
    end

    L = sim1.params.L
    q = sim1.params.q
    A = sim1.params.A
    ω = sim1.params.ω
    z0 = sim1.params.z0

    hs = [sim.params.dx for sim in sims]
    errors = zeros(length(sims))

    c = get_velocity(sim1)
    if iszero(t1)
        x_range = (-L / 2, L / 2)
    else
        x_range = (-c * t1, c * t1)
    end

    z, v, _ = ParticleMotion.oscillator(t1, z0, A, ω)
    for (k, (sim, idx)) in enumerate(zip(sims, indices))
        mask = (sim.params.x .>= x_range[1]) .& (sim.params.x .<= x_range[2])
        x = sim.params.x[mask]
        y_ret = Ψ_ret.(x, Ref(sim.params), Ref(t1)) .- ScalarField.Ψs.(x, z, v, q)
        y_sim = sim.sol.u[idx][mask, 2]
        ratio = sim1.params.dx / sim.params.dx
        ratio_int = round(Int, ratio)
        @assert isapprox(ratio, ratio_int; atol=1e-12, rtol=1e-12)
        errors[k] = norm((y_ret .- y_sim)[begin:ratio_int:end], norm_mode)
    end

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    xticks = convergence_xticks(hs)
    if norm_mode == 2
        title_str = L"\Vert \Psi_{\textrm{ret}}-\Psi_\textrm{h} \Vert_2"
    elseif norm_mode == 1
        title_str = L"\Vert \Psi_{\textrm{ret}}-\Psi_\textrm{h} \Vert_1"
    else
        title_str = L"\Vert \Psi_{\textrm{ret}}-\Psi_\textrm{h} \Vert_{\infty}"
    end

    ax = Axis(fig[1, 1];
              xlabel=L"dx",
              yscale=log10,
              xscale=log10,
              xticks=xticks,
              ytickformat=plain_log_ticklabels,
              xreversed=true,
              title=title_str)

    lines!(ax, hs, errors)
    scatter!(ax, hs, errors)

    nref = length(hs) #min(3, length(hs))
    xref = hs[1:nref]
    yref1 = xref .^ exponent
    shift1 = (errors[1] / yref1[1]) * 2.5
    yref1 .*= shift1

    lines!(ax, xref, yref1; linestyle=:dash, label=L"h^{%$exponent}")

    Legend(fig[1, 1], ax;
           halign=:right,
           valign=:top,
           tellwidth=false,
           tellheight=false,
           margin=(10, 10, 10, 10),
           framevisible=false,
           labelsize=12)

    ax.ylabel = L"\mathrm{Error}"

    rowsize!(fig.layout, 1, Auto())
    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 1, Relative(0.9))

    savename = simple_name("psifield_retarded_convergence", t1)
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end

    return fig
end

function plot_psifield_error_convergence(i, Nosc, dx0, cfl, vmax, q; norm_mode=Inf,
                                         total=false,
                                         paper=false)
    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q)
    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q)
    sim8 = Run.forced_motion(Nosc, dx0 / 8, cfl, vmax, q)
    sim16 = Run.forced_motion(Nosc, dx0 / 16, cfl, vmax, q)
    sim32 = Run.forced_motion(Nosc, dx0 / 32, cfl, vmax, q)
    t = sim1.sol.t[i]
    return plot_psifield_error_convergence(t, sim1, sim2, sim4, sim8, sim16, sim32;
                                           norm_mode=norm_mode, total=total, paper=paper)
end

end
