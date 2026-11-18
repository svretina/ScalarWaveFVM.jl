module EnergyBalanceIngredientPlots

using ..PlottingUtils
using ..ParticlePotentialPlots
using ..InitialData
using ..ODE
using ..Equations
using ..ScalarField
using CairoMakie
using DataInterpolations
using LinearAlgebra

const GLMAKIE_PKGID = Base.PkgId(Base.UUID("e9467ef8-e4e7-5192-8a1a-b1aee30e663a"),
                                 "GLMakie")

function glmakie_module()
    if isdefined(@__MODULE__, :GLMakie)
        return getfield(@__MODULE__, :GLMakie)
    end

    try
        return Base.require(GLMAKIE_PKGID)
    catch err
        throw(ArgumentError("GLMakie failed to load for the interactive energy-balance plots. Ensure GLMakie is installed and that an OpenGL display is available in this Julia session. Original error: $(err)"))
    end
end

function activate_glmakie!(theme=nothing)
    mod = glmakie_module()
    Base.invokelatest(getproperty(mod, :activate!))
    if !isnothing(theme)
        set_theme!(theme)
    else
        set_theme!()
    end
    return nothing
end

# Only math/plotting helpers live at module scope. Every plot function below
# computes the quantities it needs directly, with the bookkeeping kept local to
# the plot that uses it.

function integrate_over_domain(values, h)
    # Cell-centered FV data are interpreted as piecewise-constant/averaged values,
    # and we integrate them with the trapezoidal helper already used elsewhere.
    return ParticlePotentialPlots.trapezoidal(values, h)
end

function time_derivative(values, ts)
    # Reconstruct a smooth time signal and differentiate it at the saved timesteps.
    interpolant = PCHIPInterpolation(values, ts)
    return [DataInterpolations.derivative(interpolant, t) for t in ts]
end

function _windowed_quadratic_derivative(values, ts, t0; half_window=3)
    n = length(ts)
    n == length(values) ||
        throw(DimensionMismatch("values and ts must have the same length"))

    if n == 1
        return zero(promote_type(eltype(values), eltype(ts)))
    elseif n == 2
        return (values[2] - values[1]) / (ts[2] - ts[1])
    end

    center = clamp(searchsortedfirst(ts, t0), 2, n - 1)
    target_len = min(n, 2 * half_window + 1)
    start = max(1, center - half_window)
    stop = min(n, start + target_len - 1)
    start = max(1, stop - target_len + 1)

    radius = max(t0 - ts[start], ts[stop] - t0)
    if radius <= zero(radius)
        return (values[stop] - values[start]) / (ts[stop] - ts[start])
    end

    T = promote_type(eltype(values), eltype(ts))
    s0 = zero(T)
    s1 = zero(T)
    s2 = zero(T)
    s3 = zero(T)
    s4 = zero(T)
    sy = zero(T)
    sdy = zero(T)
    sd2y = zero(T)

    @inbounds for j in start:stop
        dt = ts[j] - t0
        u = abs(dt) / radius
        w = u < one(u) ? (one(u) - u^3)^3 : zero(u)
        dt2 = dt * dt
        s0 += w
        s1 += w * dt
        s2 += w * dt2
        s3 += w * dt2 * dt
        s4 += w * dt2 * dt2
        sy += w * values[j]
        sdy += w * dt * values[j]
        sd2y += w * dt2 * values[j]
    end

    M = T[s0 s1 s2;
          s1 s2 s3;
          s2 s3 s4]
    rhs = T[sy, sdy, sd2y]

    try
        coeffs = M \ rhs
        return coeffs[2]
    catch
        return (values[stop] - values[start]) / (ts[stop] - ts[start])
    end
end

function time_derivative(values, ts, sample_ts)
    # On custom grids such as midpoint samples, differentiate a local quadratic fit
    # instead of the piecewise interpolant itself to reduce spike-like artifacts.
    return [_windowed_quadratic_derivative(values, ts, t) for t in sample_ts]
end

function value_at_particle(values, x, z, interpolation_method)
    # Sample a field quantity at the particle worldline x = z(t).
    interpolant = interpolation_method(values, x;
                                       extrapolation=ExtrapolationType.Extension)
    return interpolant(z)
end

function _require_coupled_field_particle_sim(sim, fname)
    u0 = first(sim.sol.u)
    if !(hasproperty(u0, :x) && length(u0.x) >= 2)
        throw(ArgumentError("$(fname) expects a coupled field+particle simulation with state entries `(field_state, particle_state)`, but received states of type $(typeof(u0))."))
    end
    return nothing
end

# Plot the regular-field energy E_r(t) = ∫_Ω 1/2 (Π_r^2 + Ψ_r^2) dx.
function plot_regular_field_energy(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx

    # Build the finite control volume Ω used in the balance-law plots.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end

    regular_field_energy = zeros(nt)

    # At each saved time, integrate the regular-field density
    #   e_r = 1/2 (Π_r^2 + Ψ_r^2)
    # over the masked FV cells.
    for i in eachindex(ts)
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy[i] = integrate_over_domain(regular_density, h)
    end

    # Plot the energy history directly in this function so the full workflow is local.
    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Regular Field Energy}",
              xlabel=L"t",
              ylabel=L"E_r")
    lines!(ax, ts, regular_field_energy;
           label=L"E_r", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("regular_field_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the cross term E_×(t) = ∫_Ω (Π_r Π_s + Ψ_r Ψ_s) dx.
function plot_cross_field_energy(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q

    # Build the control volume Ω and the corresponding cell centers inside it.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    cross_field_energy = zeros(nt)

    # For each time, reconstruct the singular cell averages on Ω and integrate
    # the cross density
    #   e_× = Π_r Π_s + Ψ_r Ψ_s.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        cross_density = @. Πr * Πs + Ψr * Ψs
        cross_field_energy[i] = integrate_over_domain(cross_density, h)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Cross Field Energy}",
              xlabel=L"t",
              ylabel=L"E_\times")
    lines!(ax, ts, cross_field_energy;
           label=L"E_\times", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("cross_field_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the singular-field energy E_s(t) = ∫_Ω 1/2 (Π_s^2 + Ψ_s^2) dx.
function plot_singular_field_energy(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q

    # Build the control volume Ω and the cell centers on which puncture averages
    # are evaluated.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    singular_field_energy = zeros(nt)

    # At each time, compute the puncture cell averages and integrate
    #   e_s = 1/2 (Π_s^2 + Ψ_s^2).
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        singular_density = @. 0.5 * (Πs^2 + Ψs^2)
        singular_field_energy[i] = integrate_over_domain(singular_density, h)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Singular Field Energy}",
              xlabel=L"t",
              ylabel=L"E_s")
    lines!(ax, ts, singular_field_energy;
           label=L"E_s", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("singular_field_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the literal total-field energy E_r + E_× + E_s on the chosen control volume.
function plot_total_field_energy(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q

    # Build the finite control volume Ω used for the balance-law diagnostics.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    total_field_energy = zeros(nt)

    # Combine the regular, cross, and singular densities
    #   e_tot = 1/2(Π_r^2 + Ψ_r^2) + Π_rΠ_s + Ψ_rΨ_s + 1/2(Π_s^2 + Ψ_s^2)
    # and integrate them over Ω.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        total_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs + 0.5 * (Πs^2 + Ψs^2)
        total_field_energy[i] = integrate_over_domain(total_density, h)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Total Field Energy}",
              xlabel=L"t",
              ylabel=L"E_{\mathrm{tot,field}}")
    lines!(ax, ts, total_field_energy;
           label=L"E_{\mathrm{tot,field}}", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("total_field_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Decompose the full field energy into regular, cross, and puncture pieces and
# compare their changes over time:
#   ΔE_field^full = ΔE_r + ΔE_cross + ΔE_punct.
# The puncture piece uses the closed-form boundary expression
#   E_punct = q^2/8 (1 + v^2) [(z-a)/(1+v)^2 + (b-z)/(1-v)^2]
# on the current control volume [a, b].
function plot_full_field_energy_decomposition(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q

    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    a = x[left_idx] - 0.5 * h
    b = x[right_idx] + 0.5 * h

    regular_field_energy = zeros(nt)
    cross_field_energy = zeros(nt)
    puncture_field_energy = zeros(nt)
    total_field_energy = zeros(nt)

    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        cross_density = @. Πr * Πs + Ψr * Ψs

        regular_field_energy[i] = integrate_over_domain(regular_density, h)
        cross_field_energy[i] = integrate_over_domain(cross_density, h)
        puncture_field_energy[i] = 0.125 * q^2 * (1.0 + v * v) *
                                   ((z - a) / (1.0 + v)^2 + (b - z) / (1.0 - v)^2)
        total_field_energy[i] = regular_field_energy[i] + cross_field_energy[i] +
                                puncture_field_energy[i]
    end

    delta_regular = regular_field_energy .- regular_field_energy[1]
    delta_cross = cross_field_energy .- cross_field_energy[1]
    delta_puncture = puncture_field_energy .- puncture_field_energy[1]
    delta_total = total_field_energy .- total_field_energy[1]
    decomposition_sum = delta_regular .+ delta_cross .+ delta_puncture
    decomposition_residual = delta_total .- decomposition_sum

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(360, 600), figure_padding=(1, 5, 5, 1))
    colors = mytheme.palette.color[]

    ax_energy = Axis(fig[1, 1];
                     title=L"\textrm{Field Energy Pieces}",
                     xlabel=L"t",
                     ylabel=L"E")
    lines!(ax_energy, ts, regular_field_energy;
           label=L"E_r",
           color=colors[1], linewidth=1.5)
    lines!(ax_energy, ts, cross_field_energy;
           label=L"E_{\times}",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax_energy, ts, puncture_field_energy;
           label=L"E_{\mathrm{punct}}",
           color=colors[3], linestyle=:dot, linewidth=1.5)
    lines!(ax_energy, ts, total_field_energy;
           label=L"E_{\mathrm{field}}^{\mathrm{full}}",
           color=colors[4], linestyle=:dashdot, linewidth=1.5)
    # The puncture term can be very large and nearly constant, which makes the
    # generic scientific-notation helper generate pathological tick ranges.
    ax_energy.yminorticksvisible = false

    ax_delta = Axis(fig[3, 1];
                    title=L"\textrm{Field-Energy Change Decomposition}",
                    xlabel=L"t",
                    ylabel=L"\Delta E")
    lines!(ax_delta, ts, delta_total;
           label=L"\Delta E_{\mathrm{field}}^{\mathrm{full}}",
           color=colors[4], linewidth=1.5)
    lines!(ax_delta, ts, delta_regular;
           label=L"\Delta E_r",
           color=colors[1], linestyle=:dash, linewidth=1.5)
    lines!(ax_delta, ts, delta_cross;
           label=L"\Delta E_{\times}",
           color=colors[2], linestyle=:dot, linewidth=1.5)
    lines!(ax_delta, ts, delta_puncture;
           label=L"\Delta E_{\mathrm{punct}}",
           color=colors[3], linestyle=:dashdot, linewidth=1.5)
    lines!(ax_delta, ts, decomposition_sum;
           label=L"\Delta E_r + \Delta E_{\times} + \Delta E_{\mathrm{punct}}",
           color=colors[5], linewidth=1.5)
    ax_delta.yminorticksvisible = false

    ax_residual = Axis(fig[5, 1];
                       title=L"\textrm{Decomposition Residual}",
                       xlabel=L"t",
                       ylabel=L"\Delta E")
    lines!(ax_residual, ts, decomposition_residual;
           label=L"\Delta E_{\mathrm{field}}^{\mathrm{full}} - (\Delta E_r + \Delta E_{\times} + \Delta E_{\mathrm{punct}})",
           color=:black, linewidth=1.3)
    ax_residual.yminorticksvisible = false

    Legend(fig[2, 1], ax_energy; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    Legend(fig[4, 1], ax_delta; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    Legend(fig[6, 1], ax_residual; orientation=:horizontal, tellwidth=false,
           tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.08))
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Relative(0.08))
    rowsize!(fig.layout, 5, Auto())
    rowsize!(fig.layout, 6, Relative(0.08))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("full_field_energy_decomposition")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the particle-energy change alongside the regular, cross, and puncture
# field-energy changes, together with the total change of their sum:
#   ΔE_p, ΔE_r, ΔE_×, ΔE_punct, Δ(E_p + E_r + E_× + E_punct).
function plot_particle_full_field_energy_change_decomposition(sim; paper=false,
                                                              boundary=nothing)
    _require_coupled_field_particle_sim(sim,
                                        "plot_particle_full_field_energy_change_decomposition")

    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the same finite control volume Ω used in the field-energy diagnostics.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    a = x[left_idx] - 0.5 * h
    b = x[right_idx] + 0.5 * h

    particle_energy = zeros(nt)
    regular_field_energy = zeros(nt)
    cross_field_energy = zeros(nt)
    puncture_field_energy = zeros(nt)
    total_field_energy = zeros(nt)

    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        # Reconstruct the singular/particle field on Ω at the current particle state.
        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        # Compute the three field-energy pieces E_r, E_×, and E_punct at this time.
        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        cross_density = @. Πr * Πs + Ψr * Ψs
        regular_field_energy[i] = integrate_over_domain(regular_density, h)
        cross_field_energy[i] = integrate_over_domain(cross_density, h)
        puncture_field_energy[i] = 0.125 * q^2 * (1.0 + v * v) *
                                   ((z - a) / (1.0 + v)^2 + (b - z) / (1.0 - v)^2)
        total_field_energy[i] = regular_field_energy[i] + cross_field_energy[i] +
                                puncture_field_energy[i]

        # Compute the particle-energy candidate E_p = (m0 - q Φr(z))(γ - 1) + V(z).
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        particle_energy[i] = (m0 - q * Φrz) * (γ - 1.0) + V(z)
    end

    delta_particle = particle_energy .- particle_energy[1]
    delta_regular = regular_field_energy .- regular_field_energy[1]
    delta_cross = cross_field_energy .- cross_field_energy[1]
    delta_puncture = puncture_field_energy .- puncture_field_energy[1]

    total_candidate = particle_energy .+ total_field_energy
    delta_total_candidate = total_candidate .- total_candidate[1]

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(380, 250), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Energy-Change Decomposition}",
              xlabel=L"t",
              ylabel=L"\Delta E")

    colors = mytheme.palette.color[]
    lines!(ax, ts, delta_particle;
           label=L"\Delta E_p",
           color=colors[1], linewidth=1.5)
    lines!(ax, ts, delta_regular;
           label=L"\Delta E_r",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax, ts, delta_cross;
           label=L"\Delta E_{\times}",
           color=colors[3], linestyle=:dot, linewidth=1.5)
    lines!(ax, ts, delta_puncture;
           label=L"\Delta E_{\mathrm{punct}}",
           color=colors[4], linestyle=:dashdot, linewidth=1.5)
    lines!(ax, ts, delta_total_candidate;
           label=L"\Delta(E_p + E_r + E_{\times} + E_{\mathrm{punct}})",
           color=:black, linewidth=1.3)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.2))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("particle_full_field_energy_change_decomposition")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the regular-plus-cross field energy E_r + E_× used in the finite conserved candidate.
function plot_regular_cross_field_energy(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q

    # Restrict the FV bookkeeping to the same control volume Ω used elsewhere.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    regular_cross_field_energy = zeros(nt)

    # Integrate
    #   e_{r+×} = 1/2(Π_r^2 + Ψ_r^2) + Π_rΠ_s + Ψ_rΨ_s
    # on the masked domain.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        regular_cross_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs
        regular_cross_field_energy[i] = integrate_over_domain(regular_cross_density, h)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Regular + Cross Field Energy}",
              xlabel=L"t",
              ylabel=L"E_r + E_\times")
    lines!(ax, ts, regular_cross_field_energy;
           label=L"E_r + E_\times", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("regular_cross_field_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the external potential energy V(z(t)) of the particle.
function plot_potential_energy(sim; paper=false, boundary=nothing)
    _ = boundary
    ts = sim.sol.t
    nt = length(ts)
    k = sim.params.k
    V(x) = InitialData.Harmonic(x, k)

    potential_energy = zeros(nt)

    # The external potential depends only on the particle position z(t).
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        potential_energy[i] = V(z)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Potential Energy}",
              xlabel=L"t",
              ylabel=L"V(z)")
    lines!(ax, ts, potential_energy;
           label=L"V(z)", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("potential_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the dressed kinetic part (m_0 - qΦ_r(z))(γ - 1).
function plot_dressed_mass_energy(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    params = sim.sol.prob.p

    # Restrict the regular field to the balance-law control volume before
    # sampling Φ_r on the particle worldline.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    dressed_mass_energy = zeros(nt)

    # Evaluate
    #   E_dress = (m_0 - q Φ_r(z))(γ - 1)
    # at each saved time.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Φr = @view Q_field[mask, 3]
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        dressed_mass_energy[i] = (m0 - q * Φrz) * (γ - 1.0)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Dressed Mass Energy}",
              xlabel=L"t",
              ylabel=L"E_{\mathrm{dress}}")
    lines!(ax, ts, dressed_mass_energy;
           label=L"(m_0 - q\Phi_r(z))(\gamma - 1)",
           color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("dressed_mass_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the effective mass m_eff(t) = m_0 - qΦ_r(z(t)).
function plot_effective_mass(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    params = sim.sol.prob.p

    # Restrict Φ_r to Ω before interpolating it to the particle location.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    effective_mass = zeros(nt)

    # Sample Φ_r(z(t)) and form the effective mass
    #   m_eff = m_0 - q Φ_r(z).
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        Q_field = sim.sol.u[i].x[1]
        Φr = @view Q_field[mask, 3]
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        effective_mass[i] = m0 - q * Φrz
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Effective Mass}",
              xlabel=L"t",
              ylabel=L"m_{\mathrm{eff}}")
    lines!(ax, ts, effective_mass;
           label=L"m_0 - q\Phi_r(z)", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("effective_mass")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the particle energy together with its ingredients and the particle/worldline
# power term associated with the Hamiltonian-like particle energy
#   E_p = γ(m_0 + q Φ_r(z)) + V(z),
# whose source term is
#   P_part = q Π_r(z) / γ.
function plot_particle_energy(sim; paper=false, boundary=nothing)
    u0 = first(sim.sol.u)
    if !(hasproperty(u0, :x) && length(u0.x) >= 2)
        throw(ArgumentError("plot_particle_energy expects a coupled field+particle simulation with state entries `(field_state, particle_state)`, but received states of type $(typeof(u0))."))
    end

    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    h = hasproperty(params, :h) ? params.h : sim.params.dx
    V(x) = InitialData.Harmonic(x, k)

    # Restrict the regular field to Ω before sampling it on the worldline.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    particle_energy = zeros(nt)
    kinetic_term = zeros(nt)
    field_coupling_term = zeros(nt)
    singular_field_coupling_term = zeros(nt)
    potential_term = zeros(nt)
    non_radiative_energy = zeros(nt)
    regular_flux = zeros(nt)
    singular_phi_source_rate = zeros(nt)
    total_particle_rhs_rate = zeros(nt)
    particle_power = zeros(nt)
    non_radiative_sol = ParticlePotentialPlots.non_radiative_reference_solution(sim)
    z_arr = zeros(nt)
    v_arr = zeros(nt)

    # Build the particle energy directly from gathered regular fields and gathered
    # singular cell averages, so the worldline diagnostics use the same split-cell
    # FV sampling story as the evolution code.
    for i in eachindex(ts)
        m = sim.sol.u[i].x[2][1]
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        z_arr[i] = z
        v_arr[i] = v
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Φrz_fv = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        Φsz_fv = ODE.gather_muscl(ScalarField.cell_averages_Φs(x, h, z, v, q), z, params)
        kinetic_term[i] = m * γ
        field_coupling_term[i] = q * Φrz_fv * γ
        singular_field_coupling_term[i] = q * Φsz_fv * γ

        potential_term[i] = V(z)
        particle_energy[i] = kinetic_term[i] + 0.0field_coupling_term[i] +
                             0.0 * singular_field_coupling_term[i] + potential_term[i]

        z_nr = non_radiative_sol.u[i][2]
        v_nr = non_radiative_sol.u[i][3]
        γ_nr = inv(sqrt(1.0 - v_nr * v_nr))
        non_radiative_energy[i] = m0 * γ_nr + V(z_nr)

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        regular_flux[i] = -q_right_plus[1] * q_right_plus[2] +
                          q_left_minus[1] * q_left_minus[2]
        particle_power[i] = q * Πrz / γ
    end

    # Rebuild the exact FV average of the puncture source term added to dtΦ_r,
    # namely the split-cell average of -a ∂ᵥΦ_s, using finite-difference
    # accelerations and midpoint-averaged particle state to avoid interpolation spikes.
    dt_steps = diff(ts)
    t_mid = @. 0.5 * (ts[1:(end - 1)] + ts[2:end])
    z_mid = @. 0.5 * (z_arr[1:(end - 1)] + z_arr[2:end])
    v_mid = @. 0.5 * (v_arr[1:(end - 1)] + v_arr[2:end])
    γ_mid = @. inv(sqrt(1.0 - v_mid * v_mid))
    particle_power_mid = @. 0.5 * (particle_power[1:(end - 1)] + particle_power[2:end])
    singular_phi_source_rate_mid = zeros(length(t_mid))
    total_particle_rhs_rate_mid = zeros(length(t_mid))
    for i in eachindex(t_mid)
        z = z_mid[i]
        v = v_mid[i]
        γ = γ_mid[i]
        a = (v_arr[i + 1] - v_arr[i]) / dt_steps[i]
        dvΦsz_fv = ODE.gather_muscl(ScalarField.cell_averages_dvΦs(x, h, z, v, q), z,
                                    params)
        singular_phi_source_rate_mid[i] = q * γ * a * dvΦsz_fv
        total_particle_rhs_rate_mid[i] = particle_power_mid[i] +
                                         singular_phi_source_rate_mid[i]
    end

    # Compare against a locally smoothed midpoint derivative of the assembled
    # particle energy so the rate curve shares the midpoint grid without spurious
    # interpolation/differencing spikes.
    d_particle_energy_dt = time_derivative(particle_energy, ts, t_mid)

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(300, 380), figure_padding=(1, 5, 5, 1))
    colors = mytheme.palette.color[]

    ax_energy = Axis(fig[1, 1];
                     title=L"\textrm{Particle Energy}",
                     xlabel=L"t",
                     ylabel=L"E_p")
    lines!(ax_energy, ts, particle_energy .- 1;
           label=L"E_p(t)",
           color=colors[1], linewidth=1.5)
    lines!(ax_energy, ts, kinetic_term .- 1;
           label=L"m_0\gamma",
           color=colors[2], linestyle=:solid, linewidth=1.4)
    lines!(ax_energy, ts, field_coupling_term;
           label=L"q\Phi_r^{\mathrm{FV}}(z)\gamma",
           color=colors[3], linestyle=:solid, linewidth=1.4)
    lines!(ax_energy, ts, singular_field_coupling_term;
           label=L"q\Phi_s^{\mathrm{FV}}(z)\gamma",
           color=colors[6], linestyle=:dash, linewidth=1.4)
    lines!(ax_energy, ts, potential_term;
           label=L"V(z)",
           color=colors[4], linestyle=:solid, linewidth=1.4)
    lines!(ax_energy, ts, non_radiative_energy .- 1;
           label=L"E_{\mathrm{nonrad}}(t)=m_0\gamma_{\mathrm{nr}}+V(z_{\mathrm{nr}})",
           color=colors[5], linestyle=:solid, linewidth=1.4)
    # lines!(ax_energy, ts, radiated_energy;
    #        label=L"E_{\mathrm{rad}}(t)=\int_0^t F_r(\tau)\,d\tau",
    #        color=colors[5], linewidth=1.4)
    apply_scientific_notation_if_needed!(ax_energy)

    ax_power = Axis(fig[3, 1];
                    title=L"\textrm{Particle Power Term}",
                    xlabel=L"t",
                    ylabel=L"P_{\mathrm{part}}")
    lines!(ax_power, t_mid, particle_power_mid;
           label=L"q\Pi_r(z)/\gamma",
           color=colors[6], linewidth=1.5)
    lines!(ax_power, t_mid, singular_phi_source_rate_mid;
           label=L"q\gamma a\,\partial_v\Phi_s^{\mathrm{FV}}(z)",
           color=colors[4], linestyle=:dot, linewidth=1.4)
    lines!(ax_power, t_mid, total_particle_rhs_rate_mid;
           label=L"q\Pi_r(z)/\gamma + q\gamma a\,\partial_v\Phi_s^{\mathrm{FV}}(z)",
           color=colors[3], linestyle=:dash, linewidth=1.4)
    lines!(ax_power, t_mid, d_particle_energy_dt;
           label=L"\frac{dE_p}{dt}",
           color=colors[7], linestyle=:dash, linewidth=1.4)
    apply_scientific_notation_if_needed!(ax_power)

    Legend(fig[2, 1], ax_energy; orientation=:horizontal, tellwidth=false, tellheight=true)
    Legend(fig[4, 1], ax_power; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.1))
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Relative(0.1))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("particle_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Compare the secular loss of particle energy against the accumulated full-field
# storage-plus-radiation budget built from the outward total boundary flux
#   F_out(t) = j_tot(x_R) - j_tot(x_L),
#   j_tot = -(Π_r + Π_s)(Ψ_r + Ψ_s),
# and
#   E_rad(t) = ∫_0^t F_out(τ)dτ.
# The residual tested here is
#   R(t) = E_p(t) - E_p(0) + (E_field^full(t) - E_field^full(0)) + E_rad(t).
function plot_particle_energy_vs_radiated_energy(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    particle_energy = zeros(nt)
    total_field_energy = zeros(nt)
    total_flux = zeros(nt)

    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        Φr = @view Q_field[mask, 3]
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        particle_energy[i] = (m0 - q * Φrz) * (γ - 1.0) + V(z)

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end
        total_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs + 0.5 * (Πs^2 + Ψs^2)
        total_field_energy[i] = integrate_over_domain(total_density, h)

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)

        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))

        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        total_flux[i] = -Πtot_right_face * Ψtot_right_face +
                        Πtot_left_face * Ψtot_left_face
    end

    particle_energy_loss = particle_energy[1] .- particle_energy
    field_energy_change = total_field_energy .- total_field_energy[1]
    radiated_energy = ParticlePotentialPlots.cumulative_trapezoidal(ts, total_flux)
    recovered_energy = field_energy_change .+ radiated_energy
    diagnostic_residual = particle_energy .- particle_energy[1] .+ field_energy_change .+
                          radiated_energy

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(340, 430), figure_padding=(1, 5, 5, 1))
    colors = mytheme.palette.color[]

    ax_energy = Axis(fig[1, 1];
                     title=L"\textrm{Particle Loss vs Radiated Energy}",
                     xlabel=L"t",
                     ylabel=L"\mathrm{Energy}")
    lines!(ax_energy, ts, particle_energy_loss;
           label=L"E_p(0) - E_p(t)",
           color=colors[1], linewidth=1.5)
    lines!(ax_energy, ts, field_energy_change;
           label=L"E_{\mathrm{field}}^{\mathrm{full}}(t)-E_{\mathrm{field}}^{\mathrm{full}}(0)",
           color=colors[2], linestyle=:solid, linewidth=1.5)
    lines!(ax_energy, ts, radiated_energy;
           label=L"E_{\mathrm{rad}}(t)=\int_0^t F_{\mathrm{out}}(\tau)\,d\tau",
           color=colors[3], linestyle=:solid, linewidth=1.5)
    lines!(ax_energy, ts, recovered_energy;
           label=L"\Delta E_{\mathrm{field}}^{\mathrm{full}} + E_{\mathrm{rad}}",
           color=colors[4], linestyle=:solid, linewidth=1.5)
    apply_scientific_notation_if_needed!(ax_energy)

    ax_residual = Axis(fig[3, 1];
                       title=L"\textrm{Radiation Diagnostic Residual}",
                       xlabel=L"t",
                       ylabel=L"\mathrm{Energy}")
    lines!(ax_residual, ts, diagnostic_residual;
           label=L"E_p(t)-E_p(0)+\Delta E_{\mathrm{field}}^{\mathrm{full}}+E_{\mathrm{rad}}(t)",
           color=:black, linewidth=1.3)
    apply_scientific_notation_if_needed!(ax_residual)

    Legend(fig[2, 1], ax_energy; orientation=:horizontal, tellwidth=false, tellheight=true)
    Legend(fig[4, 1], ax_residual; orientation=:horizontal, tellwidth=false,
           tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.1))
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Relative(0.1))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("particle_energy_vs_radiated_energy")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the candidate particle Hamiltonian γ(m_0 - qΦ_r(z)) + V(z).
function plot_total_particle_energy_candidate(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Restrict the regular field to the same control volume Ω and interpolate it
    # to the particle position.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    total_particle_energy_candidate = zeros(nt)

    # Form the Hamiltonian-like particle contribution
    #   E_part,cand = γ(m_0 - q Φ_r(z)) + V(z).
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Φr = @view Q_field[mask, 3]
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        total_particle_energy_candidate[i] = γ * (m0 - q * Φrz) + V(z)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Particle Energy Candidate}",
              xlabel=L"t",
              ylabel=L"E_{\mathrm{part,cand}}")
    lines!(ax, ts, total_particle_energy_candidate;
           label=L"\gamma(m_0 - q\Phi_r(z)) + V(z)",
           color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("total_particle_energy_candidate")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the finite conserved-energy candidate E_{r+×} + γ(m_0 - qΦ_r(z)) + V(z).
function plot_conserved_total_energy_candidate(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Restrict all field terms to the finite control volume Ω.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    conserved_total_energy_candidate = zeros(nt)

    # Build
    #   E_cons = E_{r+×} + γ(m_0 - q Φ_r(z)) + V(z)
    # directly at each saved time.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        regular_cross_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs
        regular_cross_energy = integrate_over_domain(regular_cross_density, h)
        particle_energy_candidate = γ * (m0 - q * Φrz) + V(z)
        conserved_total_energy_candidate[i] = regular_cross_energy +
                                              particle_energy_candidate
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Conserved Quantity Candidate}",
              xlabel=L"t",
              ylabel=L"E_{\mathrm{cons}}")
    lines!(ax, ts, conserved_total_energy_candidate;
           label=L"E_{r+\times} + \gamma(m_0 - q\Phi_r(z)) + V(z)",
           color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("conserved_total_energy_candidate")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the boundary flux built only from the reconstructed regular field, -Π_r Ψ_r.
function plot_flux_regular(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    params = sim.sol.prob.p

    # Build the finite control volume Ω and identify the FV cells touching its
    # left and right boundaries.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)
    flux_regular = zeros(nt)

    # Reconstruct one-sided regular traces at the two control-volume faces and form
    #   F_r = j_r(x_R) - j_r(x_L),  j_r = -Π_r Ψ_r.
    for i in eachindex(ts)
        Q_field = sim.sol.u[i].x[1]
        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)

        Πr_left_face = q_left_minus[1]
        Ψr_left_face = q_left_minus[2]
        Πr_right_face = q_right_plus[1]
        Ψr_right_face = q_right_plus[2]
        flux_regular[i] = -Πr_right_face * Ψr_right_face + Πr_left_face * Ψr_left_face
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Regular Flux}",
              xlabel=L"t",
              ylabel=L"\mathcal{F}_r")
    lines!(ax, ts, flux_regular;
           label=L"\mathcal{F}_r", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("flux_regular")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the full total-field boundary flux -(Π_r + Π_s)(Ψ_r + Ψ_s).
function plot_flux_total(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    params = sim.sol.prob.p

    # Build the control volume Ω and the fixed positions of its left/right faces.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)
    flux_total = zeros(nt)

    # Combine reconstructed regular traces with puncture face traces to form
    #   F_tot = j_tot(x_R) - j_tot(x_L),
    #   j_tot = -(Π_r + Π_s)(Ψ_r + Ψ_s).
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]
        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)

        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))

        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        flux_total[i] = -Πtot_right_face * Ψtot_right_face + Πtot_left_face * Ψtot_left_face
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Total Flux}",
              xlabel=L"t",
              ylabel=L"\mathcal{F}_{\mathrm{tot}}")
    lines!(ax, ts, flux_total;
           label=L"\mathcal{F}_{\mathrm{tot}}", color=mytheme.palette.color[][1],
           linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("flux_total")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the boundary flux of the regular-plus-cross part, j_r + j_×.
function plot_flux_regular_cross(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    params = sim.sol.prob.p

    # Build the control volume Ω and its left/right FV faces.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)
    flux_regular_cross = zeros(nt)

    # Combine the regular flux with the cross pieces:
    #   j_r + j_× = -Π_rΨ_r - Π_rΨ_s - Π_sΨ_r.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]
        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)

        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))

        Πr_left_face = q_left_minus[1]
        Ψr_left_face = q_left_minus[2]
        Πr_right_face = q_right_plus[1]
        Ψr_right_face = q_right_plus[2]
        flux_regular_cross[i] = -(Πr_right_face * Ψr_right_face +
                                  Πr_right_face * Ψs_right_face +
                                  Πs_right_face * Ψr_right_face) +
                                (Πr_left_face * Ψr_left_face +
                                 Πr_left_face * Ψs_left_face +
                                 Πs_left_face * Ψr_left_face)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Regular + Cross Flux}",
              xlabel=L"t",
              ylabel=L"\mathcal{F}_{r+\times}")
    lines!(ax, ts, flux_regular_cross;
           label=L"\mathcal{F}_{r+\times}", color=mytheme.palette.color[][1],
           linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("flux_regular_cross")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the bulk source-work ∫_Ω (-Π_r R_2 - Ψ_r R_3) dx from the singular terms.
function plot_work_singular_terms(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q

    # Build the finite control volume Ω and the masked cell centers inside it.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    # Recover the particle acceleration a = dv/dt because the singular source
    # averages depend on it.
    v_arr = [sim.sol.u[i].x[2][3] for i in eachindex(ts)]
    a_arr = time_derivative(v_arr, ts)
    work_singular_terms = zeros(nt)

    # Compute the source-work density
    #   Π_r S_Π + Ψ_r S_Ψ
    # and integrate it over Ω.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        source_density = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end
        work_singular_terms[i] = h * sum(source_density)
    end

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Work of Singular Terms}",
              xlabel=L"t",
              ylabel=L"W_{\mathrm{src}}")
    lines!(ax, ts, work_singular_terms;
           label=L"\int_{\Omega}(-\Pi_r R_2 - \Psi_r R_3)\,dx",
           color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("work_singular_terms")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the explicit work slot for the exact total field, which is identically zero here.
function plot_work_total_field(sim; paper=false, boundary=nothing)
    _ = boundary
    ts = sim.sol.t

    # The exact total-field energy has no separate bulk-work term in this plot;
    # the bookkeeping slot is therefore identically zero.
    work_total_field = zeros(length(ts))

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Total-Field Work}",
              xlabel=L"t",
              ylabel=L"W_{\mathrm{tot}}")
    lines!(ax, ts, work_total_field;
           label=L"W_{\mathrm{tot}} = 0", color=mytheme.palette.color[][1], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("work_total_field")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the time derivatives of the main field and particle energy ingredients on one axis.
function plot_energy_derivatives(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the finite control volume Ω and the cell centers used to sample the
    # regular field and the puncture field.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    regular_field_energy = zeros(nt)
    cross_field_energy = zeros(nt)
    singular_field_energy = zeros(nt)
    regular_cross_field_energy = zeros(nt)
    total_field_energy = zeros(nt)
    potential_energy = zeros(nt)
    effective_mass = zeros(nt)
    dressed_mass_energy = zeros(nt)
    particle_energy = zeros(nt)
    regular_particle_energy = zeros(nt)
    conserved_total_energy_candidate = zeros(nt)

    # Compute every energy ingredient directly in this function so the derivative
    # plot reads as a self-contained balance-law worksheet.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        cross_density = @. Πr * Πs + Ψr * Ψs
        singular_density = @. 0.5 * (Πs^2 + Ψs^2)

        regular_field_energy[i] = integrate_over_domain(regular_density, h)
        cross_field_energy[i] = integrate_over_domain(cross_density, h)
        singular_field_energy[i] = integrate_over_domain(singular_density, h)
        regular_cross_field_energy[i] = integrate_over_domain(regular_density .+
                                                              cross_density,
                                                              h)
        total_field_energy[i] = integrate_over_domain(regular_density .+ cross_density .+
                                                      singular_density, h)

        potential_energy[i] = V(z)
        effective_mass[i] = m0 - q * Φrz
        dressed_mass_energy[i] = effective_mass[i] * (γ - 1.0)
        particle_energy[i] = dressed_mass_energy[i] + potential_energy[i]
        regular_particle_energy[i] = regular_field_energy[i] + particle_energy[i]
        conserved_total_energy_candidate[i] = regular_cross_field_energy[i] +
                                              γ * effective_mass[i] + potential_energy[i]
    end

    # Differentiate the time histories so the relative size of each ingredient is
    # visible on a common rate axis.
    d_regular_field_energy_dt = time_derivative(regular_field_energy, ts)
    d_cross_field_energy_dt = time_derivative(cross_field_energy, ts)
    d_singular_field_energy_dt = time_derivative(singular_field_energy, ts)
    d_regular_cross_field_energy_dt = time_derivative(regular_cross_field_energy, ts)
    d_total_field_energy_dt = time_derivative(total_field_energy, ts)
    d_potential_energy_dt = time_derivative(potential_energy, ts)
    d_effective_mass_dt = time_derivative(effective_mass, ts)
    d_dressed_mass_energy_dt = time_derivative(dressed_mass_energy, ts)
    d_particle_energy_dt = time_derivative(particle_energy, ts)
    d_regular_particle_energy_dt = time_derivative(regular_particle_energy, ts)
    d_conserved_total_energy_candidate_dt = time_derivative(conserved_total_energy_candidate,
                                                            ts)

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(300, 220), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Time Derivatives of Energies}",
              xlabel=L"t",
              ylabel=L"\partial_t E")

    # Plot all differentiated ingredients together so we can compare which terms
    # are active in the bookkeeping at a given time.
    colors = mytheme.palette.color[]
    lines!(ax, ts, d_regular_field_energy_dt;
           label=L"\partial_t E_r", color=colors[1], linewidth=1.5)
    lines!(ax, ts, d_cross_field_energy_dt;
           label=L"\partial_t E_\times", color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax, ts, d_singular_field_energy_dt;
           label=L"\partial_t E_s", color=colors[3], linestyle=:dot, linewidth=1.4)
    lines!(ax, ts, d_regular_cross_field_energy_dt;
           label=L"\partial_t(E_r + E_\times)", color=colors[4], linewidth=1.4)
    lines!(ax, ts, d_total_field_energy_dt;
           label=L"\partial_t E_{\mathrm{tot,field}}", color=colors[5], linestyle=:dashdot,
           linewidth=1.3)
    lines!(ax, ts, d_potential_energy_dt;
           label=L"\partial_t V(z)", color=colors[6], linewidth=1.3)
    lines!(ax, ts, d_effective_mass_dt;
           label=L"\partial_t(m_0 - q\Phi_r(z))", color=:gray55, linestyle=:dot,
           linewidth=1.2)
    lines!(ax, ts, d_dressed_mass_energy_dt;
           label=L"\partial_t E_{\mathrm{dress}}", color=:black, linestyle=:dashdot,
           linewidth=1.2)
    lines!(ax, ts, d_particle_energy_dt;
           label=L"\partial_t E_p", color=:gray30, linewidth=1.2)
    lines!(ax, ts, d_regular_particle_energy_dt;
           label=L"\partial_t(E_p + E_r)", color=:black, linestyle=:dash,
           linewidth=1.3)
    lines!(ax, ts, d_conserved_total_energy_candidate_dt;
           label=L"\partial_t E_{\mathrm{cons}}", color=:firebrick, linewidth=1.3)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.18))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("energy_derivatives")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the differential balance d(E_p + E_r)/dt = -F_r + W_src and its residual.
function plot_particle_regular_rate_balance(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the control volume Ω and the cell centers used in the source-work term.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    # Recover the particle acceleration needed by the singular source averages.
    v_arr = [sim.sol.u[i].x[2][3] for i in eachindex(ts)]
    a_arr = time_derivative(v_arr, ts)

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    regular_particle_energy = zeros(nt)
    flux_term_rate = zeros(nt)
    work_term_rate = zeros(nt)

    # Compute the three ingredients entering the rate balance:
    #   E_p + E_r,  -F_r,  and  W_src.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy = integrate_over_domain(regular_density, h)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        particle_energy = (m0 - q * Φrz) * (γ - 1.0) + V(z)
        regular_particle_energy[i] = regular_field_energy + particle_energy

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        flux_term_rate[i] = -(-q_right_plus[1] * q_right_plus[2] +
                              q_left_minus[1] * q_left_minus[2])

        source_density = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end
        work_term_rate[i] = h * sum(source_density)
    end

    # Differentiate the particle-plus-regular-field energy and compare it against
    # the separated RHS ingredients.
    d_regular_particle_energy_dt = time_derivative(regular_particle_energy, ts)
    rhs_rate = flux_term_rate .+ work_term_rate
    rate_residual = d_regular_particle_energy_dt .- rhs_rate

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(340, 235), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Rate Balance: } E_p + E_r",
              xlabel=L"t",
              ylabel=L"\mathrm{Rate}")

    colors = mytheme.palette.color[]
    lines!(ax, ts, d_regular_particle_energy_dt;
           label=L"\frac{d(E_p + E_r)}{dt}",
           color=colors[1], linewidth=1.5)
    lines!(ax, ts, flux_term_rate;
           label=L"-\mathcal{F}_r",
           color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax, ts, work_term_rate;
           label=L"W_{\mathrm{src}}",
           color=colors[3], linestyle=:dot, linewidth=1.4)
    lines!(ax, ts, rhs_rate;
           label=L"-\mathcal{F}_r + W_{\mathrm{src}}",
           color=colors[4], linestyle=:dashdot, linewidth=1.5)
    lines!(ax, ts, rate_residual;
           label=L"\mathrm{Residual}",
           color=:black, linewidth=1.2)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.2))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("particle_regular_rate_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the integrated balance Δ(E_p + E_r) against the accumulated flux and source work.
function plot_particle_regular_integral_balance(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the control volume Ω and the face/cell data needed in the integral law.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    v_arr = [sim.sol.u[i].x[2][3] for i in eachindex(ts)]
    a_arr = time_derivative(v_arr, ts)

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    regular_particle_energy = zeros(nt)
    flux_term_rate = zeros(nt)
    work_term_rate = zeros(nt)

    # Compute the same ingredients as in the rate plot, but keep their full time
    # histories so we can integrate them.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy = integrate_over_domain(regular_density, h)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        particle_energy = (m0 - q * Φrz) * (γ - 1.0) + V(z)
        regular_particle_energy[i] = regular_field_energy + particle_energy

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        flux_term_rate[i] = -(-q_right_plus[1] * q_right_plus[2] +
                              q_left_minus[1] * q_left_minus[2])

        source_density = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end
        work_term_rate[i] = h * sum(source_density)
    end

    # Accumulate the flux and source-work histories separately before adding them.
    energy_change = regular_particle_energy .- regular_particle_energy[1]
    cumulative_flux_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, flux_term_rate)
    cumulative_work_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, work_term_rate)
    rhs_integral = cumulative_flux_term .+ cumulative_work_term
    integral_residual = energy_change .- rhs_integral

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(340, 235), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Integral Balance: } E_p + E_r",
              xlabel=L"t",
              ylabel=L"\mathrm{Energy}")

    colors = mytheme.palette.color[]
    lines!(ax, ts, energy_change;
           label=L"\Delta(E_p + E_r)",
           color=colors[1], linewidth=1.5)
    lines!(ax, ts, cumulative_flux_term;
           label=L"\int_0^t -\mathcal{F}_r\,d\tau",
           color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax, ts, cumulative_work_term;
           label=L"\int_0^t W_{\mathrm{src}}\,d\tau",
           color=colors[3], linestyle=:dot, linewidth=1.4)
    lines!(ax, ts, rhs_integral;
           label=L"\int_0^t (-\mathcal{F}_r + W_{\mathrm{src}})\,d\tau",
           color=colors[4], linestyle=:dashdot, linewidth=1.5)
    lines!(ax, ts, integral_residual;
           label=L"\mathrm{Residual}",
           color=:black, linewidth=1.2)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.2))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("particle_regular_integral_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the conserved-energy candidate, its rate law, and its integral balance.
function plot_conserved_total_energy_balance(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the finite control volume Ω together with its boundary faces.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    regular_cross_field_energy = zeros(nt)
    total_particle_energy_candidate = zeros(nt)
    conserved_total_energy_candidate = zeros(nt)
    flux_term_rate = zeros(nt)
    work_term_rate = zeros(nt)

    # Compute the candidate conserved quantity together with the RHS terms in its
    # balance law:
    #   dE_cons/dt = -F_tot + q Π_r(z)/γ.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)

        regular_cross_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs
        regular_cross_field_energy[i] = integrate_over_domain(regular_cross_density, h)
        total_particle_energy_candidate[i] = γ * (m0 - q * Φrz) + V(z)
        conserved_total_energy_candidate[i] = regular_cross_field_energy[i] +
                                              total_particle_energy_candidate[i]

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))
        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        flux_total = -Πtot_right_face * Ψtot_right_face + Πtot_left_face * Ψtot_left_face

        flux_term_rate[i] = -flux_total
        work_term_rate[i] = q * Πrz / γ
    end

    # Build the differential and integral forms of the candidate balance law.
    d_conserved_total_energy_candidate_dt = time_derivative(conserved_total_energy_candidate,
                                                            ts)
    rhs_rate = flux_term_rate .+ work_term_rate
    rate_residual = d_conserved_total_energy_candidate_dt .- rhs_rate

    energy_change = conserved_total_energy_candidate .- conserved_total_energy_candidate[1]
    cumulative_flux_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, flux_term_rate)
    cumulative_work_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, work_term_rate)
    rhs_integral = cumulative_flux_term .+ cumulative_work_term
    integral_residual = energy_change .- rhs_integral

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(340, 540), figure_padding=(1, 5, 5, 1))
    colors = mytheme.palette.color[]

    ax_energy = Axis(fig[1, 1];
                     title=L"\textrm{Conserved Energy Candidate}",
                     ylabel=L"E")

    # Plot the two pieces entering E_cons together with their sum.
    lines!(ax_energy, ts, regular_cross_field_energy;
           label=L"E_{r+\times}",
           color=colors[1], linewidth=1.5)
    lines!(ax_energy, ts, total_particle_energy_candidate;
           label=L"\gamma(m_0 - q\Phi_r(z)) + V(z)",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax_energy, ts, conserved_total_energy_candidate;
           label=L"E_{\mathrm{cons}}",
           color=colors[4], linewidth=1.5)
    apply_scientific_notation_if_needed!(ax_energy)

    ax_rate = Axis(fig[3, 1];
                   title=L"\textrm{Rate Balance}",
                   ylabel=L"\mathrm{Rate}")

    # Compare dE_cons/dt with its flux and worldline-work ingredients.
    lines!(ax_rate, ts, d_conserved_total_energy_candidate_dt;
           label=L"\frac{dE_{\mathrm{cons}}}{dt}",
           color=colors[1], linewidth=1.5)
    lines!(ax_rate, ts, flux_term_rate;
           label=L"-\mathcal{F}_{\mathrm{tot}}",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax_rate, ts, work_term_rate;
           label=L"q\Pi_r(z)/\gamma",
           color=colors[3], linestyle=:dot, linewidth=1.4)
    lines!(ax_rate, ts, rhs_rate;
           label=L"-\mathcal{F}_{\mathrm{tot}} + q\Pi_r(z)/\gamma",
           color=colors[4], linestyle=:dashdot, linewidth=1.5)
    lines!(ax_rate, ts, rate_residual;
           label=L"\mathrm{Residual}",
           color=:black, linewidth=1.2, linestyle=:dot)
    apply_scientific_notation_if_needed!(ax_rate)

    ax_int = Axis(fig[5, 1];
                  title=L"\textrm{Integral Balance}",
                  xlabel=L"t",
                  ylabel=L"\mathrm{Energy}")

    # Compare the accumulated energy change against the accumulated RHS terms.
    lines!(ax_int, ts, energy_change;
           label=L"\Delta E_{\mathrm{cons}}",
           color=colors[1], linewidth=1.5)
    lines!(ax_int, ts, cumulative_flux_term;
           label=L"\int_0^t -\mathcal{F}_{\mathrm{tot}}\,d\tau",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax_int, ts, cumulative_work_term;
           label=L"\int_0^t q\Pi_r(z)/\gamma\,d\tau",
           color=colors[3], linestyle=:dot, linewidth=1.4)
    lines!(ax_int, ts, rhs_integral;
           label=L"\int_0^t \left(-\mathcal{F}_{\mathrm{tot}} + q\Pi_r(z)/\gamma\right)\,d\tau",
           color=colors[4], linestyle=:dashdot, linewidth=1.5)
    lines!(ax_int, ts, integral_residual;
           label=L"\mathrm{Residual}",
           color=:black, linewidth=1.2, linestyle=:dot)
    apply_scientific_notation_if_needed!(ax_int)

    Legend(fig[2, 1], ax_energy; orientation=:horizontal, tellwidth=false, tellheight=true)
    Legend(fig[4, 1], ax_rate; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    Legend(fig[6, 1], ax_int; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.08))
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Relative(0.08))
    rowsize!(fig.layout, 5, Auto())
    rowsize!(fig.layout, 6, Relative(0.08))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("conserved_total_energy_balance_candidate")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the time derivative of the conserved total-energy candidate together with q Π_r(z).
function plot_total_energy_derivative_and_work(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the finite control volume Ω and the cell centers used to evaluate
    # the regular field and the puncture field.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    conserved_total_energy_candidate = zeros(nt)
    work_term = zeros(nt)

    # Compute
    #   E_cons = E_{r+×} + γ(m_0 - q Φ_r(z)) + V(z)
    # and the worldline term requested for comparison,
    #   q Π_r(z).
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = 1#inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end
        Πsz, _ = ODE.gather_muscl(hcat(ScalarField.cell_averages_Πs(params.x, h, z, v, q),
                                       ScalarField.cell_averages_Ψs(params.x, h, z, v, q)),
                                  z, params)
        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)

        total_energy_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs +
                                  0.5 * (Πs^2 + Ψs^2)
        total_energy_density = integrate_over_domain(total_energy_density, h)
        particle_energy_candidate = γ * (m0 - q * Φrz) + V(z)

        conserved_total_energy_candidate[i] = total_energy_density +
                                              particle_energy_candidate
        work_term[i] = -(q / γ) * (Πrz + Πsz)
    end

    # Differentiate the total-energy candidate in time and compare it directly
    # against the requested worldline term q Π_r(z).
    d_conserved_total_energy_candidate_dt = time_derivative(conserved_total_energy_candidate,
                                                            ts)

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(300, 220), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Total-Energy Derivative and Work}",
              xlabel=L"t",
              ylabel=L"\mathrm{Rate}")

    colors = mytheme.palette.color[]
    lines!(ax, ts, d_conserved_total_energy_candidate_dt;
           label=L"\frac{dE_{\mathrm{cons}}}{dt}",
           color=colors[1], linewidth=1.5)
    lines!(ax, ts, work_term;
           label=L"q\Pi_r(z)",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.16))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("total_energy_derivative_and_work")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the integrated balance test
#   R(t) = E_{\mathrm{tot}}(t) - E_{\mathrm{tot}}(0) + \int_0^t W(\tau)\,d\tau
# using the same total-energy candidate and work term as the derivative plot above.
function plot_total_energy_integrated_work_balance(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the finite control volume Ω and the cell centers used to evaluate
    # the regular field and the puncture field.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]

    conserved_total_energy_candidate = zeros(nt)
    work_term = zeros(nt)

    # Compute the same total-energy candidate and worldline work term used in the
    # differential diagnostic, but keep them in integral form here.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)

        regular_cross_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs
        regular_cross_energy = integrate_over_domain(regular_cross_density, h)
        particle_energy_candidate = γ * (m0 - q * Φrz) + V(z)

        conserved_total_energy_candidate[i] = regular_cross_energy +
                                              particle_energy_candidate
        work_term[i] = q * Πrz
    end

    # Form the integrated residual from the image:
    #   R(t) = E_tot(t) - E_tot(0) + ∫_0^t W(τ) dτ.
    energy_change = conserved_total_energy_candidate .- conserved_total_energy_candidate[1]
    cumulative_work_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, work_term)
    integrated_residual = energy_change .+ cumulative_work_term

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(320, 235), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Integrated Total-Energy Test}",
              xlabel=L"t",
              ylabel=L"\mathrm{Energy}")

    colors = mytheme.palette.color[]
    # Plot the two pieces entering the integrated test together with the residual
    # so it is easy to see whether the cancellation is good.
    lines!(ax, ts, energy_change;
           label=L"E_{\mathrm{tot}}(t) - E_{\mathrm{tot}}(0)",
           color=colors[1], linewidth=1.5)
    lines!(ax, ts, cumulative_work_term;
           label=L"\int_0^t W(\tau)\,d\tau",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax, ts, integrated_residual;
           label=L"R(t)",
           color=:black, linewidth=1.3)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.18))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("total_energy_integrated_work_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the balance law for E_{r+\times}, including the boundary flux, the
# worldline source term -q Π_r(z), and the acceleration-dependent puncture term.
function plot_regular_cross_balance_law(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    params = sim.sol.prob.p

    # Build the finite control volume Ω together with its two boundary faces.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    # Recover the particle acceleration a_p = dv/dt because the puncture source
    # averages depend on it.
    v_arr = [sim.sol.u[i].x[2][3] for i in eachindex(ts)]
    a_arr = time_derivative(v_arr, ts)

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    regular_cross_field_energy = zeros(nt)
    flux_term_rate = zeros(nt)
    point_source_rate = zeros(nt)
    acceleration_term_rate = zeros(nt)

    # Compute
    #   E_{r+\times} = ∫ [1/2(Π_r^2 + Ψ_r^2) + Π_rΠ_p + Ψ_rΨ_p] dx
    # and the three RHS ingredients appearing in the picture:
    #   [Π_rΨ_r + Π_rΨ_p + Π_pΨ_r]_{x_L}^{x_R},
    #   -q Π_r(z),
    #   -a_p ∫ (Π_r ∂_vΠ_p + Ψ_r ∂_vΨ_p) dx.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        regular_cross_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs
        regular_cross_field_energy[i] = integrate_over_domain(regular_cross_density, h)

        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        point_source_rate[i] = -q * Πrz

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))

        # In this module the stored flux uses the outward convention
        #   F_{r+\times} = -(f_R - f_L),
        # so the bracketed term from the picture is -F_{r+\times}.
        flux_regular_cross = -(q_right_plus[1] * q_right_plus[2] +
                               q_right_plus[1] * Ψs_right_face +
                               Πs_right_face * q_right_plus[2]) +
                             (q_left_minus[1] * q_left_minus[2] +
                              q_left_minus[1] * Ψs_left_face +
                              Πs_left_face * q_left_minus[2])
        flux_term_rate[i] = -flux_regular_cross

        # The existing singular source averages are exactly the acceleration-
        # dependent puncture term in cell-averaged form.
        acceleration_density = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            acceleration_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end
        acceleration_term_rate[i] = h * sum(acceleration_density)
    end

    # Compare the derivative of E_{r+\times} against the separated RHS pieces and
    # their sum.
    d_regular_cross_field_energy_dt = time_derivative(regular_cross_field_energy, ts)
    rhs_rate = flux_term_rate .+ point_source_rate .+ acceleration_term_rate
    rate_residual = d_regular_cross_field_energy_dt .- rhs_rate

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(360, 250), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Balance Law: } E_{r+\times}",
              xlabel=L"t",
              ylabel=L"\mathrm{Rate}")

    colors = mytheme.palette.color[]
    lines!(ax, ts, d_regular_cross_field_energy_dt;
           label=L"\frac{dE_{r+\times}}{dt}",
           color=colors[1], linewidth=1.5)
    lines!(ax, ts, flux_term_rate;
           label=L"\left[\Pi_r\Psi_r + \Pi_r\Psi_p + \Pi_p\Psi_r\right]_{x_L}^{x_R}",
           color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax, ts, point_source_rate;
           label=L"-q\Pi_r(z)",
           color=colors[3], linestyle=:dot, linewidth=1.4)
    lines!(ax, ts, acceleration_term_rate;
           label=L"-a_p \int_\Omega (\Pi_r\partial_v\Pi_p + \Psi_r\partial_v\Psi_p)\,dx",
           color=colors[4], linestyle=:dashdot, linewidth=1.4)
    lines!(ax, ts, rhs_rate;
           label=L"\mathrm{RHS}",
           color=colors[5], linewidth=1.5)
    lines!(ax, ts, rate_residual;
           label=L"\mathrm{Residual}",
           color=:black, linewidth=1.2)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.22))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("regular_cross_balance_law")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Plot the integrated residual test for E_p + E_r when the boundary term is
# assumed negligible:
#   R_{p+r}(t) = E_{p+r}(t) - E_{p+r}(0)
#              - ∫_0^t F_r^{disc}(τ) dτ + ∫_0^t W_{\mathrm{src}}^{disc}(τ) dτ.
function plot_regular_integrated_residual(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    v_arr = [sim.sol.u[i].x[2][3] for i in eachindex(ts)]
    a_arr = time_derivative(v_arr, ts)

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    total_energy_candidate = zeros(nt)
    work_term = zeros(nt)
    flux_term = zeros(nt)
    particle_energy = zeros(nt)
    gamma = zeros(nt)
    phir = zeros(nt)
    pot = zeros(nt)

    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        γ = inv(sqrt(1.0 - v * v))
        gamma[i] = γ

        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy = integrate_over_domain(regular_density, h)

        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        phir[i] = 0m0 - (γ - 1) * q * Φrz
        particle_energy[i] = (γ - 1) * (m0 - q * Φrz) + V(z)
        pot[i] = V(z)
        total_energy_candidate[i] = regular_field_energy + 0particle_energy[i]

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        flux_regular = -q_right_plus[1] * q_right_plus[2] +
                       q_left_minus[1] * q_left_minus[2]
        flux_term[i] = -flux_regular

        source_density = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end

        work_term[i] = h * sum(source_density)
    end

    energy_change = total_energy_candidate .- total_energy_candidate[1]
    cumulative_flux_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, flux_term)
    cumulative_work_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, work_term)
    integrated_residual = energy_change .- cumulative_flux_term .- cumulative_work_term

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(340, 235), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Integrated Residual: } E_p + E_r",
              xlabel=L"t",
              ylabel=L"\mathrm{Energy}")

    lines!(ax, ts, integrated_residual;
           label=L"R_{p+r}(t)",
           color=:black, linewidth=1.3)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.2))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("regular_integrated_residual")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_regular_integrated_residual_resolutions(sim1, sim2, sim4; paper=false,
                                                      boundary=nothing)
    ts1 = sim1.sol.t
    nt1 = length(ts1)
    x1 = sim1.params.x
    h1 = sim1.params.dx
    q1 = sim1.params.q
    m01 = hasproperty(sim1.params, :m0) ? sim1.params.m0 : sim1.params.m
    k1 = sim1.params.k
    params1 = sim1.sol.prob.p
    V1(x) = InitialData.Harmonic(x, k1)

    cutoff1 = isnothing(boundary) ? ts1[end] : boundary
    mask1 = (x1 .>= -cutoff1) .& (x1 .<= cutoff1)
    if iszero(ts1[end]) || !any(mask1)
        mask1 = trues(length(x1))
    end
    x_filtered1 = x1[mask1]
    left_idx1 = findfirst(mask1)
    right_idx1 = findlast(mask1)

    v_arr1 = [sim1.sol.u[i].x[2][3] for i in eachindex(ts1)]
    a_arr1 = time_derivative(v_arr1, ts1)

    L1 = Equations.left_eigenmatrix(params1.equation)
    R1 = Equations.right_eigenmatrix(params1.equation)

    total_energy_candidate1 = zeros(nt1)
    work_term1 = zeros(nt1)
    flux_term1 = zeros(nt1)
    particle_energy1 = zeros(nt1)
    gamma1 = zeros(nt1)
    phir1 = zeros(nt1)
    pot1 = zeros(nt1)

    for i in eachindex(ts1)
        z = sim1.sol.u[i].x[2][2]
        v = v_arr1[i]
        a = a_arr1[i]
        γ = inv(sqrt(1.0 - v * v))
        gamma1[i] = γ

        Q_field = sim1.sol.u[i].x[1]
        Πr = @view Q_field[mask1, 1]
        Ψr = @view Q_field[mask1, 2]
        Φr = @view Q_field[mask1, 3]

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy = integrate_over_domain(regular_density, h1)

        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params1)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params1)
        phir1[i] = 0m01 - (γ - 1) * q1 * Φrz
        particle_energy1[i] = (γ - 1) * (m01 - q1 * Φrz) + V1(z)
        pot1[i] = V1(z)
        total_energy_candidate1[i] = regular_field_energy + 0particle_energy1[i]

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx1, params1.N, L1, R1,
                                               params1.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx1, params1.N, L1, R1,
                                               params1.slope_limiter)
        flux_regular = -q_right_plus[1] * q_right_plus[2] +
                       q_left_minus[1] * q_left_minus[2]
        flux_term1[i] = -flux_regular

        source_density = zeros(length(x_filtered1))
        @inbounds for (idx, x_i) in pairs(x_filtered1)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h1, z, v, q1),
                          -a * ScalarField.cell_average_dvΨs(x_i, h1, z, v, q1)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end

        work_term1[i] = h1 * sum(source_density)
    end

    energy_change1 = total_energy_candidate1 .- total_energy_candidate1[1]
    cumulative_flux_term1 = ParticlePotentialPlots.cumulative_trapezoidal(ts1, flux_term1)
    cumulative_work_term1 = ParticlePotentialPlots.cumulative_trapezoidal(ts1, work_term1)
    residual1 = energy_change1 .- cumulative_flux_term1 .- cumulative_work_term1

    ts2 = sim2.sol.t
    nt2 = length(ts2)
    x2 = sim2.params.x
    h2 = sim2.params.dx
    q2 = sim2.params.q
    m02 = hasproperty(sim2.params, :m0) ? sim2.params.m0 : sim2.params.m
    k2 = sim2.params.k
    params2 = sim2.sol.prob.p
    V2(x) = InitialData.Harmonic(x, k2)

    cutoff2 = isnothing(boundary) ? ts2[end] : boundary
    mask2 = (x2 .>= -cutoff2) .& (x2 .<= cutoff2)
    if iszero(ts2[end]) || !any(mask2)
        mask2 = trues(length(x2))
    end
    x_filtered2 = x2[mask2]
    left_idx2 = findfirst(mask2)
    right_idx2 = findlast(mask2)

    v_arr2 = [sim2.sol.u[i].x[2][3] for i in eachindex(ts2)]
    a_arr2 = time_derivative(v_arr2, ts2)

    L2 = Equations.left_eigenmatrix(params2.equation)
    R2 = Equations.right_eigenmatrix(params2.equation)

    total_energy_candidate2 = zeros(nt2)
    work_term2 = zeros(nt2)
    flux_term2 = zeros(nt2)
    particle_energy2 = zeros(nt2)
    gamma2 = zeros(nt2)
    phir2 = zeros(nt2)
    pot2 = zeros(nt2)

    for i in eachindex(ts2)
        z = sim2.sol.u[i].x[2][2]
        v = v_arr2[i]
        a = a_arr2[i]
        γ = inv(sqrt(1.0 - v * v))
        gamma2[i] = γ

        Q_field = sim2.sol.u[i].x[1]
        Πr = @view Q_field[mask2, 1]
        Ψr = @view Q_field[mask2, 2]
        Φr = @view Q_field[mask2, 3]

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy = integrate_over_domain(regular_density, h2)

        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params2)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params2)
        phir2[i] = 0m02 - (γ - 1) * q2 * Φrz
        particle_energy2[i] = (γ - 1) * (m02 - q2 * Φrz) + V2(z)
        pot2[i] = V2(z)
        total_energy_candidate2[i] = regular_field_energy + 0particle_energy2[i]

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx2, params2.N, L2, R2,
                                               params2.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx2, params2.N, L2, R2,
                                               params2.slope_limiter)
        flux_regular = -q_right_plus[1] * q_right_plus[2] +
                       q_left_minus[1] * q_left_minus[2]
        flux_term2[i] = -flux_regular

        source_density = zeros(length(x_filtered2))
        @inbounds for (idx, x_i) in pairs(x_filtered2)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h2, z, v, q2),
                          -a * ScalarField.cell_average_dvΨs(x_i, h2, z, v, q2)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end

        work_term2[i] = h2 * sum(source_density)
    end

    energy_change2 = total_energy_candidate2 .- total_energy_candidate2[1]
    cumulative_flux_term2 = ParticlePotentialPlots.cumulative_trapezoidal(ts2, flux_term2)
    cumulative_work_term2 = ParticlePotentialPlots.cumulative_trapezoidal(ts2, work_term2)
    residual2_full = energy_change2 .- cumulative_flux_term2 .- cumulative_work_term2

    ts4 = sim4.sol.t
    nt4 = length(ts4)
    x4 = sim4.params.x
    h4 = sim4.params.dx
    q4 = sim4.params.q
    m04 = hasproperty(sim4.params, :m0) ? sim4.params.m0 : sim4.params.m
    k4 = sim4.params.k
    params4 = sim4.sol.prob.p
    V4(x) = InitialData.Harmonic(x, k4)

    cutoff4 = isnothing(boundary) ? ts4[end] : boundary
    mask4 = (x4 .>= -cutoff4) .& (x4 .<= cutoff4)
    if iszero(ts4[end]) || !any(mask4)
        mask4 = trues(length(x4))
    end
    x_filtered4 = x4[mask4]
    left_idx4 = findfirst(mask4)
    right_idx4 = findlast(mask4)

    v_arr4 = [sim4.sol.u[i].x[2][3] for i in eachindex(ts4)]
    a_arr4 = time_derivative(v_arr4, ts4)

    L4 = Equations.left_eigenmatrix(params4.equation)
    R4 = Equations.right_eigenmatrix(params4.equation)

    total_energy_candidate4 = zeros(nt4)
    work_term4 = zeros(nt4)
    flux_term4 = zeros(nt4)
    particle_energy4 = zeros(nt4)
    gamma4 = zeros(nt4)
    phir4 = zeros(nt4)
    pot4 = zeros(nt4)

    for i in eachindex(ts4)
        z = sim4.sol.u[i].x[2][2]
        v = v_arr4[i]
        a = a_arr4[i]
        γ = inv(sqrt(1.0 - v * v))
        gamma4[i] = γ

        Q_field = sim4.sol.u[i].x[1]
        Πr = @view Q_field[mask4, 1]
        Ψr = @view Q_field[mask4, 2]
        Φr = @view Q_field[mask4, 3]

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy = integrate_over_domain(regular_density, h4)

        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params4)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params4)
        phir4[i] = 0m04 - (γ - 1) * q4 * Φrz
        particle_energy4[i] = (γ - 1) * (m04 - q4 * Φrz) + V4(z)
        pot4[i] = V4(z)
        total_energy_candidate4[i] = regular_field_energy + 0particle_energy4[i]

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx4, params4.N, L4, R4,
                                               params4.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx4, params4.N, L4, R4,
                                               params4.slope_limiter)
        flux_regular = -q_right_plus[1] * q_right_plus[2] +
                       q_left_minus[1] * q_left_minus[2]
        flux_term4[i] = -flux_regular

        source_density = zeros(length(x_filtered4))
        @inbounds for (idx, x_i) in pairs(x_filtered4)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h4, z, v, q4),
                          -a * ScalarField.cell_average_dvΨs(x_i, h4, z, v, q4)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end

        work_term4[i] = h4 * sum(source_density)
    end

    energy_change4 = total_energy_candidate4 .- total_energy_candidate4[1]
    cumulative_flux_term4 = ParticlePotentialPlots.cumulative_trapezoidal(ts4, flux_term4)
    cumulative_work_term4 = ParticlePotentialPlots.cumulative_trapezoidal(ts4, work_term4)
    residual4_full = energy_change4 .- cumulative_flux_term4 .- cumulative_work_term4

    residual2 = [residual2_full[2i - 1] for i in 1:nt1]
    residual4 = [residual4_full[4i - 3] for i in 1:nt1]

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(340, 235), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Integrated Residual: } E_p + E_r",
              xlabel=L"t",
              ylabel=L"\mathrm{Energy}")

    colors = mytheme.palette.color[]
    lines!(ax, ts1, residual1;
           label=L"\textrm{low res}",
           color=colors[2], linestyle=:solid, linewidth=1.6)
    lines!(ax, ts1, residual2;
           label=L"\textrm{mid res}",
           color=colors[3], linestyle=:dash, linewidth=1.4)
    lines!(ax, ts1, residual4;
           label=L"\textrm{high res}",
           color=colors[4], linestyle=:dot, linewidth=1.4)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.2))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("regular_integrated_residual_resolutions")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Copy of the regular integrated residual plot with the additional particle/worldline
# contribution needed by the E_p + E_r balance:
#   P_p = q[(1 - 1/γ) Π_r(z) + v Ψ_r(z)].
function plot_regular_integrated_residual_test(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)

    v_arr = [sim.sol.u[i].x[2][3] for i in eachindex(ts)]
    a_arr = time_derivative(v_arr, ts)

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    total_energy_candidate = zeros(nt)
    flux_term = zeros(nt)
    work_term = zeros(nt)
    particle_term = zeros(nt)

    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = v_arr[i]
        a = a_arr[i]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        regular_density = @. 0.5 * (Πr^2 + Ψr^2)
        regular_field_energy = integrate_over_domain(regular_density, h)

        Πrz, Ψrz = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        particle_energy = (γ - 1) * (m0 - q * Φrz) + V(z)
        total_energy_candidate[i] = regular_field_energy + particle_energy

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        flux_regular = -q_right_plus[1] * q_right_plus[2] +
                       q_left_minus[1] * q_left_minus[2]
        flux_term[i] = -flux_regular

        source_density = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            S_Pi, S_Psi = -a * ScalarField.cell_average_dvΠs(x_i, h, z, v, q),
                          -a * ScalarField.cell_average_dvΨs(x_i, h, z, v, q)
            source_density[idx] = Πr[idx] * S_Pi + Ψr[idx] * S_Psi
        end
        work_term[i] = h * sum(source_density)
        particle_term[i] = q * ((1.0 - inv(γ)) * Πrz + v * Ψrz)
    end

    energy_change = total_energy_candidate .- total_energy_candidate[1]
    cumulative_flux_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, flux_term)
    cumulative_work_term = ParticlePotentialPlots.cumulative_trapezoidal(ts, work_term)
    cumulative_particle_term = ParticlePotentialPlots.cumulative_trapezoidal(ts,
                                                                             particle_term)
    rhs_integral = cumulative_flux_term .+ cumulative_work_term .+ cumulative_particle_term
    integrated_residual = energy_change .- rhs_integral

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(360, 245), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Integrated Residual Test: } E_p + E_r",
              xlabel=L"t",
              ylabel=L"\mathrm{Energy}")

    colors = mytheme.palette.color[]
    lines!(ax, ts, energy_change;
           label=L"\Delta(E_p + E_r)",
           color=colors[1], linewidth=1.5)
    lines!(ax, ts, cumulative_flux_term;
           label=L"\int_0^t -\mathcal{F}_{r}(\tau)\,d\tau",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax, ts, cumulative_work_term;
           label=L"\int_0^t W_{\mathrm{src}}(\tau)\,d\tau",
           color=colors[3], linestyle=:dot, linewidth=1.5)
    lines!(ax, ts, cumulative_particle_term;
           label=L"\int_0^t P_{\mathrm{part}}(\tau)\,d\tau",
           color=colors[4], linestyle=:dashdot, linewidth=1.5)
    lines!(ax, ts, integrated_residual;
           label=L"R_{p+r}^{\mathrm{test}}(t)",
           color=:black, linewidth=1.3)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.22))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("regular_integrated_residual_test")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Backward-compatible alias used by prepare_energy_balance_ingredient_plots.
function plot_regular_cross_integrated_residual(sim; paper=false, boundary=nothing)
    return plot_regular_integrated_residual(sim; paper=paper, boundary=boundary)
end

# Plot the exact field, particle, and total rate balances using the total-field work term.
function plot_total_field_rate_balance(sim; paper=false, boundary=nothing)
    _require_coupled_field_particle_sim(sim, "plot_total_field_rate_balance")

    source_ts = sim.sol.t
    length(source_ts) >= 2 ||
        throw(ArgumentError("Need at least 2 sampled times to build the total-field rate-balance plot; got only $(length(source_ts))."))
    nt = length(source_ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    params = sim.sol.prob.p

    # Build the finite control volume Ω and all face indices used in the balance law.
    cutoff = isnothing(boundary) ? source_ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(source_ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h
    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    field_energy_total = zeros(nt)
    flux_out = zeros(nt)
    work_trace_plain = zeros(nt)

    for i in eachindex(source_ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]

        # Rebuild the singular field on Ω and integrate the full field energy.
        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end
        total_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs + 0.5 * (Πs^2 + Ψs^2)
        field_energy_total[i] = h * sum(total_density)

        # Correct the cell containing the particle by splitting the puncture field at z(t).
        if left_face_x <= z <= right_face_x
            local_cell = clamp(floor(Int, (z - left_face_x) / h) + 1, 1, length(x_filtered))
            cell_center = x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                frac_L = (z - cell_left_face) / h
                frac_R = one(frac_L) - frac_L
                Πs_left, Πs_right, Ψs_left, Ψs_right = ScalarField.singular_particle_traces(z,
                                                                                            v,
                                                                                            q)
                Πr_cell = Πr[local_cell]
                Ψr_cell = Ψr[local_cell]
                split_cell_density = 0.5 *
                                     (frac_L *
                                      ((Πr_cell + Πs_left)^2 + (Ψr_cell + Ψs_left)^2) +
                                      frac_R *
                                      ((Πr_cell + Πs_right)^2 + (Ψr_cell + Ψs_right)^2))
                field_energy_total[i] += h *
                                         (split_cell_density - total_density[local_cell])
            end
        end

        # Compute the outward total-field flux through the two boundary faces of Ω.
        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))
        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        flux_out[i] = -Πtot_right_face * Ψtot_right_face + Πtot_left_face * Ψtot_left_face

        # Compute the worldline work term -q Π̄_tot^tr(z) from the split trace reconstruction.
        Πtrace_avg = zero(eltype(Πr))
        if left_face_x <= z <= right_face_x
            local_cell = clamp(floor(Int, (z - left_face_x) / h) + 1, 1, length(x_filtered))
            cell_center = x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                frac_L = (z - cell_left_face) / h
                frac_R = one(frac_L) - frac_L
                Πs_left, Πs_right, _, _ = ScalarField.singular_particle_traces(z, v, q)
                global_cell = left_idx + local_cell - 1
                q_particle_minus, q_particle_plus = ODE.reconstruct_cell(Q_field,
                                                                         global_cell,
                                                                         params.N, L, R,
                                                                         params.slope_limiter)
                Πr_trace = frac_R * q_particle_minus[1] + frac_L * q_particle_plus[1]
                Πtot_left_trace = Πr_trace + Πs_left
                Πtot_right_trace = Πr_trace + Πs_right
                Πtrace_avg = 0.5 * (Πtot_left_trace + Πtot_right_trace)
            end
        end
        work_trace_plain[i] = -q * Πtrace_avg
    end

    # Convert the raw histories into the finite-difference and midpoint rates shown in the plot.
    dt_steps = diff(source_ts)
    t_mid = @. 0.5 * (source_ts[1:(end - 1)] + source_ts[2:end])
    d_field_energy_dt = diff(field_energy_total) ./ dt_steps
    flux_mid = @. 0.5 * (flux_out[1:(end - 1)] + flux_out[2:end])
    work_mid_plain = @. 0.5 * (work_trace_plain[1:(end - 1)] + work_trace_plain[2:end])
    rhs_rate_plain = .-flux_mid .+ work_mid_plain
    residual_plain = d_field_energy_dt .- rhs_rate_plain

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(360, 240), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              title=L"\textrm{Rate Balance: Total Field}",
              xlabel=L"t",
              ylabel=L"\mathrm{Rate}")

    colors = mytheme.palette.color[]
    lines!(ax, t_mid, d_field_energy_dt;
           label=L"\Delta E_{\mathrm{field}}/\Delta t",
           color=colors[1], linewidth=1.5)
    lines!(ax, t_mid, .-flux_mid;
           label=L"-\mathcal{F}_{\mathrm{out}}",
           color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax, t_mid, work_mid_plain;
           label=L"-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=colors[4], linestyle=:dashdot, linewidth=1.4)
    lines!(ax, t_mid, rhs_rate_plain;
           label=L"-\mathcal{F}_{\mathrm{out}}-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=:black, linewidth=1.2)
    lines!(ax, t_mid, residual_plain;
           label=L"\mathrm{Residual}\ (q)",
           color=:gray35, linestyle=:dash, linewidth=1.1)
    apply_scientific_notation_if_needed!(ax)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.28))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)

    savename = simple_name("total_field_rate_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

function plot_interactive_field_energy_decomposition_and_rate_balance(sim; boundary=nothing)
    _require_coupled_field_particle_sim(sim,
                                        "plot_interactive_field_energy_decomposition_and_rate_balance")
    source_ts = sim.sol.t
    length(source_ts) >= 2 ||
        throw(ArgumentError("Need at least 2 sampled times to build the interactive total-field rate-balance plot; got only $(length(source_ts))."))

    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    params = sim.sol.prob.p
    tf = source_ts[end]

    # Build the finite control volume Ω used for the rate balance.
    cutoff = isnothing(boundary) ? source_ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(source_ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    # Build the full snapshot domain used by the left interactive panel.
    plot_mask = (x .>= -tf) .& (x .<= tf)
    if iszero(tf) || !any(plot_mask)
        plot_mask = trues(length(x))
    end
    plot_x_filtered = x[plot_mask]
    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    nt = length(source_ts)
    field_energy_total = zeros(nt)
    flux_out = zeros(nt)
    work_trace_plain = zeros(nt)
    plot_x_data = Vector{Vector{Float64}}(undef, nt)
    plot_regular_data = Vector{Vector{Float64}}(undef, nt)
    plot_cross_data = Vector{Vector{Float64}}(undef, nt)
    plot_singular_data = Vector{Vector{Float64}}(undef, nt)
    plot_total_data = Vector{Vector{Float64}}(undef, nt)
    left_density_min = Inf
    left_density_max = -Inf

    for i in eachindex(source_ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        Q_field = sim.sol.u[i].x[1]

        # Reconstruct the full field-energy density on Ω for the rate balance.
        Πr_balance = @view Q_field[mask, 1]
        Ψr_balance = @view Q_field[mask, 2]
        Πs_balance = zeros(length(x_filtered))
        Ψs_balance = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs_balance[idx], Ψs_balance[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end
        total_density = @. 0.5 * (Πr_balance^2 + Ψr_balance^2) +
                           Πr_balance * Πs_balance + Ψr_balance * Ψs_balance +
                           0.5 * (Πs_balance^2 + Ψs_balance^2)
        field_energy_total[i] = h * sum(total_density)

        # Split the particle cell so the total-field energy matches the exact puncture jump.
        if left_face_x <= z <= right_face_x
            local_cell = clamp(floor(Int, (z - left_face_x) / h) + 1, 1, length(x_filtered))
            cell_center = x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                frac_L = (z - cell_left_face) / h
                frac_R = one(frac_L) - frac_L
                Πs_left, Πs_right, Ψs_left, Ψs_right = ScalarField.singular_particle_traces(z,
                                                                                            v,
                                                                                            q)
                Πr_cell = Πr_balance[local_cell]
                Ψr_cell = Ψr_balance[local_cell]
                split_cell_density = 0.5 *
                                     (frac_L *
                                      ((Πr_cell + Πs_left)^2 + (Ψr_cell + Ψs_left)^2) +
                                      frac_R *
                                      ((Πr_cell + Πs_right)^2 + (Ψr_cell + Ψs_right)^2))
                field_energy_total[i] += h *
                                         (split_cell_density - total_density[local_cell])
            end
        end

        # Compute the total flux leaving Ω through the two boundary faces.
        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))
        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        flux_out[i] = -Πtot_right_face * Ψtot_right_face + Πtot_left_face * Ψtot_left_face

        # Compute the plain worldline work term -q Π̄_tot^tr(z).
        Πtrace_avg = zero(eltype(Πr_balance))
        if left_face_x <= z <= right_face_x
            local_cell = clamp(floor(Int, (z - left_face_x) / h) + 1, 1, length(x_filtered))
            cell_center = x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                frac_L = (z - cell_left_face) / h
                frac_R = one(frac_L) - frac_L
                Πs_left, Πs_right, _, _ = ScalarField.singular_particle_traces(z, v, q)
                global_cell = left_idx + local_cell - 1
                q_particle_minus, q_particle_plus = ODE.reconstruct_cell(Q_field,
                                                                         global_cell,
                                                                         params.N, L, R,
                                                                         params.slope_limiter)
                Πr_trace = frac_R * q_particle_minus[1] + frac_L * q_particle_plus[1]
                Πtot_left_trace = Πr_trace + Πs_left
                Πtot_right_trace = Πr_trace + Πs_right
                Πtrace_avg = 0.5 * (Πtot_left_trace + Πtot_right_trace)
            end
        end
        work_trace_plain[i] = -q * Πtrace_avg

        # Rebuild the plotted density pieces on the full snapshot domain used by the left panel.
        Πr_plot = @view Q_field[plot_mask, 1]
        Ψr_plot = @view Q_field[plot_mask, 2]
        Πs_plot = zeros(length(plot_x_filtered))
        Ψs_plot = zeros(length(plot_x_filtered))
        @inbounds for (idx, x_i) in pairs(plot_x_filtered)
            Πs_plot[idx], Ψs_plot[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                                         ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end
        regular_density = @. 0.5 * (Πr_plot^2 + Ψr_plot^2)
        cross_density = @. Πr_plot * Πs_plot + Ψr_plot * Ψs_plot
        singular_density = @. 0.5 * (Πs_plot^2 + Ψs_plot^2)
        total_density_plot = regular_density .+ cross_density .+ singular_density

        # Turn the cell-centered densities into the piecewise-constant curves shown in the panel.
        n_plot = length(plot_x_filtered)
        plot_left_face_x = plot_x_filtered[1] - 0.5 * h
        plot_right_face_x = plot_x_filtered[end] + 0.5 * h
        particle_cell = nothing
        if plot_left_face_x <= z <= plot_right_face_x
            local_cell = clamp(floor(Int, (z - plot_left_face_x) / h) + 1, 1, n_plot)
            cell_center = plot_x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                particle_cell = local_cell
            end
        end

        x_plot = Float64[]
        regular_plot = Float64[]
        cross_plot = Float64[]
        singular_plot = Float64[]
        total_plot = Float64[]
        for j in eachindex(plot_x_filtered)
            cell_center = plot_x_filtered[j]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            reg_val = regular_density[j]
            if !isnothing(particle_cell) && j == particle_cell
                Πs_left, Πs_right, Ψs_left, Ψs_right = ScalarField.singular_particle_traces(z,
                                                                                            v,
                                                                                            q)
                cross_left = Πr_plot[j] * Πs_left + Ψr_plot[j] * Ψs_left
                cross_right = Πr_plot[j] * Πs_right + Ψr_plot[j] * Ψs_right
                singular_left = 0.5 * (Πs_left^2 + Ψs_left^2)
                singular_right = 0.5 * (Πs_right^2 + Ψs_right^2)
                append!(x_plot, (cell_left_face, z, z, cell_right_face))
                append!(regular_plot, (reg_val, reg_val, reg_val, reg_val))
                append!(cross_plot, (cross_left, cross_left, cross_right, cross_right))
                append!(singular_plot,
                        (singular_left, singular_left, singular_right, singular_right))
                append!(total_plot,
                        (reg_val + cross_left + singular_left,
                         reg_val + cross_left + singular_left,
                         reg_val + cross_right + singular_right,
                         reg_val + cross_right + singular_right))
            else
                append!(x_plot, (cell_left_face, cell_right_face))
                append!(regular_plot, (reg_val, reg_val))
                append!(cross_plot, (cross_density[j], cross_density[j]))
                append!(singular_plot, (singular_density[j], singular_density[j]))
                append!(total_plot, (total_density_plot[j], total_density_plot[j]))
            end
        end

        plot_x_data[i] = x_plot
        plot_regular_data[i] = regular_plot
        plot_cross_data[i] = cross_plot
        plot_singular_data[i] = singular_plot
        plot_total_data[i] = total_plot

        snapshot_min = min(minimum(regular_plot), minimum(cross_plot),
                           minimum(singular_plot),
                           minimum(total_plot))
        snapshot_max = max(maximum(regular_plot), maximum(cross_plot),
                           maximum(singular_plot),
                           maximum(total_plot))
        left_density_min = min(left_density_min, snapshot_min)
        left_density_max = max(left_density_max, snapshot_max)
    end

    # Convert the histories to the midpoint rates displayed on the right panel.
    dt_steps = diff(source_ts)
    t_mid = @. 0.5 * (source_ts[1:(end - 1)] + source_ts[2:end])
    d_field_energy_dt = diff(field_energy_total) ./ dt_steps
    flux_mid = @. 0.5 * (flux_out[1:(end - 1)] + flux_out[2:end])
    work_mid_plain = @. 0.5 * (work_trace_plain[1:(end - 1)] + work_trace_plain[2:end])
    rhs_rate_plain = .-flux_mid .+ work_mid_plain
    residual_plain = d_field_energy_dt .- rhs_rate_plain

    y_pad = left_density_max > left_density_min ?
            0.05 * (left_density_max - left_density_min) :
            max(abs(left_density_max), 1.0) * 0.05
    left_density_limits = (left_density_min - y_pad, left_density_max + y_pad)
    activate_glmakie!()
    colors = [:dodgerblue3, :darkorange2, :seagreen4, :firebrick3, :mediumpurple3]

    fig = Figure(; size=(920, 620), figure_padding=(10, 10, 10, 10))

    ax_left = Axis(fig[1, 1];
                   title="Field Energy Densities",
                   xlabel=L"x",
                   ylabel=L"\mathcal{E}(x,t)",
                   titlesize=26,
                   xlabelsize=24,
                   ylabelsize=24,
                   xticklabelsize=20,
                   yticklabelsize=20)
    ax_right = Axis(fig[1, 2];
                    title=L"\textrm{Rate Balance: Total Field}",
                    xlabel=L"t",
                    ylabel=L"\mathrm{Rate}",
                    titlesize=26,
                    xlabelsize=24,
                    ylabelsize=24,
                    xticklabelsize=20,
                    yticklabelsize=20)

    slider_layout = fig[3, 1:2] = GridLayout()
    Label(slider_layout[1, 1], "Time index"; fontsize=24)
    time_slider = Slider(slider_layout[1, 2];
                         range=1:length(source_ts),
                         startvalue=1,
                         horizontal=true)
    time_label = lift(time_slider.value) do i
        "t = $(round(source_ts[i], digits=6))"
    end
    Label(slider_layout[1, 3], time_label; fontsize=24)

    current_index = time_slider.value
    current_time = lift(current_index) do i
        source_ts[i]
    end
    current_energy_change_ratio = lift(current_index) do i
        (field_energy_total[i] - field_energy_total[1]) / field_energy_total[1]
    end

    total_points = lift(i -> Point2f.(plot_x_data[i], plot_total_data[i]), current_index)
    regular_points = lift(i -> Point2f.(plot_x_data[i], plot_regular_data[i]),
                          current_index)
    cross_points = lift(i -> Point2f.(plot_x_data[i], plot_cross_data[i]), current_index)
    singular_points = lift(i -> Point2f.(plot_x_data[i], plot_singular_data[i]),
                           current_index)

    lines!(ax_left, total_points;
           label=L"\mathcal{E}_{\mathrm{tot}}", color=:black, linewidth=2.0)
    lines!(ax_left, regular_points;
           label=L"\mathcal{E}_r", color=colors[1], linewidth=1.6)
    lines!(ax_left, cross_points;
           label=L"\mathcal{E}_\times", color=colors[2], linewidth=1.6,
           linestyle=:dash)
    lines!(ax_left, singular_points;
           label=L"\mathcal{E}_s", color=colors[4], linewidth=1.6,
           linestyle=:dot)
    if !isnothing(boundary) && boundary < sim.sol.t[end]
        vlines!(ax_left, [-boundary, boundary];
                color=:firebrick, linestyle=:dash, linewidth=1.2)
    end
    on(lift((ti, er) -> (ti, er), current_time, current_energy_change_ratio)) do (ti, er)
        ax_left.title[] = "Field Energy Densities at t = $(round(ti, digits=6)) | ΔE_field / E_field(0) = $(round(er, digits=6))"
    end
    xlims!(ax_left, -tf, tf)
    ylims!(ax_left, left_density_limits...)

    lines!(ax_right, t_mid, d_field_energy_dt;
           label=L"\Delta E_{\mathrm{field}}/\Delta t",
           color=colors[1], linewidth=1.5)
    lines!(ax_right, t_mid, .-flux_mid;
           label=L"-\mathcal{F}_{\mathrm{out}}",
           color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax_right, t_mid, work_mid_plain;
           label=L"-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=colors[4], linestyle=:dashdot, linewidth=1.4)
    lines!(ax_right, t_mid, rhs_rate_plain;
           label=L"-\mathcal{F}_{\mathrm{out}}-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=:black, linewidth=1.2)
    lines!(ax_right, t_mid, residual_plain;
           label=L"\mathrm{Residual}\ (q)",
           color=:gray35, linestyle=:dash, linewidth=1.1)
    vlines!(ax_right, lift(ti -> [ti], current_time);
            color=:firebrick, linestyle=:dash, linewidth=1.2)
    xlims!(ax_right, first(source_ts), last(source_ts))

    Legend(fig[2, 1], ax_left; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, labelsize=22, colgap=20, rowgap=10, padding=(10, 10, 10, 10))
    Legend(fig[2, 2], ax_right; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=2, labelsize=22, colgap=20, rowgap=10,
           padding=(10, 10, 10, 10))
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Auto())
    rowsize!(fig.layout, 3, Fixed(42))
    colsize!(fig.layout, 1, Relative(0.5))
    colsize!(fig.layout, 2, Relative(0.5))
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 10)

    return fig
end

function plot_interactive_field_particle_energy_dashboard_legacy(sim; boundary=nothing)
    _require_coupled_field_particle_sim(sim,
                                        "plot_interactive_field_particle_energy_dashboard")
    source_ts = sim.sol.t
    length(source_ts) >= 2 ||
        throw(ArgumentError("Need at least 2 sampled times to build the interactive field-particle dashboard; got only $(length(source_ts))."))

    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)
    tf = source_ts[end]

    # Build the finite control volume Ω once so every field and particle quantity uses the same domain.
    cutoff = isnothing(boundary) ? source_ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(source_ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    # Build the full snapshot domain used by the left dashboard panel.
    plot_mask = (x .>= -tf) .& (x .<= tf)
    if iszero(tf) || !any(plot_mask)
        plot_mask = trues(length(x))
    end
    plot_x_filtered = x[plot_mask]
    plot_balance_mask = (plot_x_filtered .>= -cutoff) .& (plot_x_filtered .<= cutoff)
    if iszero(source_ts[end]) || !any(plot_balance_mask)
        plot_balance_mask = trues(length(plot_x_filtered))
    end
    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)
    a = left_face_x
    b = right_face_x

    nt = length(source_ts)
    field_energy_total = zeros(nt)
    flux_out = zeros(nt)
    work_trace_plain = zeros(nt)
    regular_field_energy = zeros(nt)
    cross_field_energy = zeros(nt)
    puncture_field_energy = zeros(nt)
    z_arr = zeros(nt)
    v_arr = zeros(nt)
    γ_arr = zeros(nt)
    Πrz_arr = zeros(nt)
    Ψrz_arr = zeros(nt)
    Φrz_fv_arr = zeros(nt)
    particle_energy = zeros(nt)
    kinetic_term = zeros(nt)
    field_coupling_term = zeros(nt)
    singular_field_coupling_term = zeros(nt)
    singular_phi_source_rate = zeros(nt)
    total_particle_rhs_rate = zeros(nt)
    Φsz_fv_arr = zeros(nt)
    potential_term = zeros(nt)
    non_radiative_energy = zeros(nt)
    non_radiative_sol = ParticlePotentialPlots.non_radiative_reference_solution(sim)
    plot_x_data = Vector{Vector{Float64}}(undef, nt)
    plot_regular_data = Vector{Vector{Float64}}(undef, nt)
    plot_cross_data = Vector{Vector{Float64}}(undef, nt)
    plot_singular_data = Vector{Vector{Float64}}(undef, nt)
    plot_total_data = Vector{Vector{Float64}}(undef, nt)
    left_density_min = Inf
    left_density_max = -Inf

    for i in eachindex(source_ts)
        m = sim.sol.u[i].x[2][1]
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        z_arr[i] = z
        v_arr[i] = v
        γ_arr[i] = γ
        Q_field = sim.sol.u[i].x[1]

        # Reconstruct the singular field on Ω and compute the full field-energy pieces.
        Πr_balance = @view Q_field[mask, 1]
        Ψr_balance = @view Q_field[mask, 2]
        Φr_balance = @view Q_field[mask, 3]
        Πs_balance = zeros(length(x_filtered))
        Ψs_balance = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs_balance[idx], Ψs_balance[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end
        regular_density_balance = @. 0.5 * (Πr_balance^2 + Ψr_balance^2)
        cross_density_balance = @. Πr_balance * Πs_balance + Ψr_balance * Ψs_balance
        total_density = @. regular_density_balance + cross_density_balance +
                           0.5 * (Πs_balance^2 + Ψs_balance^2)
        regular_field_energy[i] = integrate_over_domain(regular_density_balance, h)
        cross_field_energy[i] = integrate_over_domain(cross_density_balance, h)
        puncture_field_energy[i] = 0.125 * q^2 * (1.0 + v * v) *
                                   ((z - a) / (1.0 + v)^2 + (b - z) / (1.0 - v)^2)
        total_field_energy[i] = regular_field_energy[i] + cross_field_energy[i] +
                                puncture_field_energy[i]
        field_energy_total[i] = h * sum(total_density)

        # Split the cell containing the particle for the exact total-field energy used in the rate panel.
        if left_face_x <= z <= right_face_x
            local_cell = clamp(floor(Int, (z - left_face_x) / h) + 1, 1, length(x_filtered))
            cell_center = x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                frac_L = (z - cell_left_face) / h
                frac_R = one(frac_L) - frac_L
                Πs_left, Πs_right, Ψs_left, Ψs_right = ScalarField.singular_particle_traces(z,
                                                                                            v,
                                                                                            q)
                Πr_cell = Πr_balance[local_cell]
                Ψr_cell = Ψr_balance[local_cell]
                split_cell_density = 0.5 *
                                     (frac_L *
                                      ((Πr_cell + Πs_left)^2 + (Ψr_cell + Ψs_left)^2) +
                                      frac_R *
                                      ((Πr_cell + Πs_right)^2 + (Ψr_cell + Ψs_right)^2))
                field_energy_total[i] += h *
                                         (split_cell_density - total_density[local_cell])
            end
        end

        # Compute the boundary flux and the worldline work used in the total-field rate balance.
        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))
        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        flux_out[i] = -Πtot_right_face * Ψtot_right_face + Πtot_left_face * Ψtot_left_face

        Πtrace_avg = zero(eltype(Πr_balance))
        if left_face_x <= z <= right_face_x
            local_cell = clamp(floor(Int, (z - left_face_x) / h) + 1, 1, length(x_filtered))
            cell_center = x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                frac_L = (z - cell_left_face) / h
                frac_R = one(frac_L) - frac_L
                Πs_left, Πs_right, _, _ = ScalarField.singular_particle_traces(z, v, q)
                global_cell = left_idx + local_cell - 1
                q_particle_minus, q_particle_plus = ODE.reconstruct_cell(Q_field,
                                                                         global_cell,
                                                                         params.N, L, R,
                                                                         params.slope_limiter)
                Πr_trace = frac_R * q_particle_minus[1] + frac_L * q_particle_plus[1]
                Πtot_left_trace = Πr_trace + Πs_left
                Πtot_right_trace = Πr_trace + Πs_right
                Πtrace_avg = 0.5 * (Πtot_left_trace + Πtot_right_trace)
            end
        end
        work_trace_plain[i] = -q * Πtrace_avg

        # Sample the regular Π and Ψ fields on the worldline for the power diagnostics.
        Πrz_arr[i], Ψrz_arr[i] = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        kinetic_term[i] = m * γ

        Φrz_fv_arr[i] = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        field_coupling_term[i] = 0#-q * Φrz_fv_arr[i] * γ

        Φsz_fv_arr[i] = ODE.gather_muscl(ScalarField.cell_averages_Φs(x, h, z, v, q), z,
                                         params)
        singular_field_coupling_term[i] = -q * Φsz_fv_arr[i] * γ

        potential_term[i] = V(z)
        particle_energy[i] = kinetic_term[i] + field_coupling_term[i] +
                             singular_field_coupling_term[i] + potential_term[i]

        z_nr = non_radiative_sol.u[i][2]
        v_nr = non_radiative_sol.u[i][3]
        γ_nr = inv(sqrt(1.0 - v_nr * v_nr))
        non_radiative_energy[i] = m0 * γ_nr + V(z_nr)

        # Rebuild the left-panel snapshot data on the full plotting domain.
        Πr_plot = @view Q_field[plot_mask, 1]
        Ψr_plot = @view Q_field[plot_mask, 2]
        Πs_plot = zeros(length(plot_x_filtered))
        Ψs_plot = zeros(length(plot_x_filtered))
        @inbounds for (idx, x_i) in pairs(plot_x_filtered)
            Πs_plot[idx], Ψs_plot[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                                         ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end
        regular_density_plot = @. 0.5 * (Πr_plot^2 + Ψr_plot^2)
        cross_density_plot = @. Πr_plot * Πs_plot + Ψr_plot * Ψs_plot
        singular_density_plot = @. 0.5 * (Πs_plot^2 + Ψs_plot^2)
        total_density_plot = regular_density_plot .+ cross_density_plot .+
                             singular_density_plot

        n_plot = length(plot_x_filtered)
        plot_left_face_x = plot_x_filtered[1] - 0.5 * h
        plot_right_face_x = plot_x_filtered[end] + 0.5 * h
        particle_cell = nothing
        if plot_left_face_x <= z <= plot_right_face_x
            local_cell = clamp(floor(Int, (z - plot_left_face_x) / h) + 1, 1, n_plot)
            cell_center = plot_x_filtered[local_cell]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            if cell_left_face <= z <= cell_right_face
                particle_cell = local_cell
            end
        end

        x_plot = Float64[]
        regular_plot = Float64[]
        cross_plot = Float64[]
        singular_plot = Float64[]
        total_plot = Float64[]
        for j in eachindex(plot_x_filtered)
            cell_center = plot_x_filtered[j]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h
            reg_val = regular_density_plot[j]
            if !isnothing(particle_cell) && j == particle_cell
                Πs_left, Πs_right, Ψs_left, Ψs_right = ScalarField.singular_particle_traces(z,
                                                                                            v,
                                                                                            q)
                cross_left = Πr_plot[j] * Πs_left + Ψr_plot[j] * Ψs_left
                cross_right = Πr_plot[j] * Πs_right + Ψr_plot[j] * Ψs_right
                singular_left = 0.5 * (Πs_left^2 + Ψs_left^2)
                singular_right = 0.5 * (Πs_right^2 + Ψs_right^2)
                append!(x_plot, (cell_left_face, z, z, cell_right_face))
                append!(regular_plot, (reg_val, reg_val, reg_val, reg_val))
                append!(cross_plot, (cross_left, cross_left, cross_right, cross_right))
                append!(singular_plot,
                        (singular_left, singular_left, singular_right, singular_right))
                append!(total_plot,
                        (reg_val + cross_left + singular_left,
                         reg_val + cross_left + singular_left,
                         reg_val + cross_right + singular_right,
                         reg_val + cross_right + singular_right))
            else
                append!(x_plot, (cell_left_face, cell_right_face))
                append!(regular_plot, (reg_val, reg_val))
                append!(cross_plot, (cross_density_plot[j], cross_density_plot[j]))
                append!(singular_plot, (singular_density_plot[j], singular_density_plot[j]))
                append!(total_plot, (total_density_plot[j], total_density_plot[j]))
            end
        end
        plot_x_data[i] = x_plot
        plot_regular_data[i] = regular_plot
        plot_cross_data[i] = cross_plot
        plot_singular_data[i] = singular_plot
        plot_total_data[i] = total_plot

        snapshot_min = min(minimum(regular_plot), minimum(cross_plot),
                           minimum(singular_plot),
                           minimum(total_plot))
        snapshot_max = max(maximum(regular_plot), maximum(cross_plot),
                           maximum(singular_plot),
                           maximum(total_plot))
        left_density_min = min(left_density_min, snapshot_min)
        left_density_max = max(left_density_max, snapshot_max)
    end

    # Rebuild the same discrete puncture source average that enters dtΦ_r so the
    # dashboard can show the extra term from
    #   d/dt[γ(m0 - qΦ_r^FV)+V] = -qΠ_r/γ + qγa∂ᵥΦ_s.
    # Use nodal finite differences for a(t_i) so the diagnostic follows the same
    # spike-resistant choice as the particle-rate plot.
    a_arr = zeros(nt)
    if nt >= 2
        a_arr[1] = (v_arr[2] - v_arr[1]) / (source_ts[2] - source_ts[1])
        a_arr[end] = (v_arr[end] - v_arr[end - 1]) / (source_ts[end] - source_ts[end - 1])
    end
    for i in 2:(nt - 1)
        a_arr[i] = (v_arr[i + 1] - v_arr[i - 1]) / (source_ts[i + 1] - source_ts[i - 1])
    end
    for i in eachindex(source_ts)
        z = z_arr[i]
        v = v_arr[i]
        γ = γ_arr[i]
        a = a_arr[i]
        dvΦsz_fv = ODE.gather_muscl(ScalarField.cell_averages_dvΦs(x, h, z, v, q), z,
                                    params)
        singular_phi_source_rate[i] = q * γ * a * dvΦsz_fv
        total_particle_rhs_rate[i] = -q * Πrz_arr[i] / γ + singular_phi_source_rate[i]
    end

    # Convert the histories to the midpoint rates used in the right-hand total-field balance panel.
    dt_steps = diff(source_ts)
    t_mid = @. 0.5 * (source_ts[1:(end - 1)] + source_ts[2:end])
    d_field_energy_dt = diff(field_energy_total) ./ dt_steps
    flux_mid = @. 0.5 * (flux_out[1:(end - 1)] + flux_out[2:end])
    work_mid_plain = @. 0.5 * (work_trace_plain[1:(end - 1)] + work_trace_plain[2:end])
    rhs_rate_plain = .-flux_mid .+ work_mid_plain
    residual_plain = d_field_energy_dt .- rhs_rate_plain

    y_pad = left_density_max > left_density_min ?
            0.05 * (left_density_max - left_density_min) :
            max(abs(left_density_max), 1.0) * 0.05
    left_density_limits = (left_density_min - y_pad, left_density_max + y_pad)

    activate_glmakie!()
    colors = [:dodgerblue3, :darkorange2, :seagreen4, :firebrick3, :mediumpurple3]

    fig = Figure(; size=(1260, 860), figure_padding=(10, 6, 10, 10))

    ax_left = Axis(fig[1, 1];
                   title="Field Energy Densities",
                   xlabel=L"x",
                   ylabel=L"\mathcal{E}(x,t)",
                   titlesize=30,
                   xlabelsize=26,
                   ylabelsize=26,
                   xticklabelsize=22,
                   yticklabelsize=22)
    ax_right = Axis(fig[1, 2];
                    title=L"\textrm{Rate Balance: Total Field}",
                    xlabel=L"t",
                    ylabel=L"\mathrm{Rate}",
                    titlesize=30,
                    xlabelsize=26,
                    ylabelsize=26,
                    xticklabelsize=22,
                    yticklabelsize=22)
    ax_particle = Axis(fig[3, 1:2];
                       title=L"\textrm{Particle Energy}",
                       xlabel=L"t",
                       ylabel=L"E_p",
                       titlesize=30,
                       xlabelsize=26,
                       ylabelsize=26,
                       xticklabelsize=22,
                       yticklabelsize=22)

    diagnostics_layout = fig[1:4, 3] = GridLayout()
    Label(diagnostics_layout[1, 1], "Diagnostics";
          fontsize=30, halign=:left, tellwidth=false)

    slider_layout = fig[5, 1:2] = GridLayout()
    Label(slider_layout[1, 1], "Time index"; fontsize=24)
    time_slider = Slider(slider_layout[1, 2];
                         range=1:length(source_ts),
                         startvalue=1,
                         horizontal=true)
    time_label = lift(time_slider.value) do i
        "t = $(round(source_ts[i], digits=6))"
    end
    Label(slider_layout[1, 3], time_label; fontsize=24)

    current_index = time_slider.value
    current_time = lift(current_index) do i
        source_ts[i]
    end
    current_energy_change_ratio = lift(current_index) do i
        (field_energy_total[i] - field_energy_total[1]) / field_energy_total[1]
    end
    diagnostics_text = lift(current_index, current_time) do i, ti
        period_str = iszero(sim.params.ω) ? "N/A" :
                     PlottingUtils.format_val(ti / (2π / sim.params.ω))
        zp_str = iszero(sim.params.A) ? "N/A" :
                 PlottingUtils.format_val(z_arr[i] / sim.params.A)
        particle_ratio_str = iszero(particle_energy[1]) ? "N/A" :
                             PlottingUtils.format_val((particle_energy[i] -
                                                       particle_energy[1]) /
                                                      particle_energy[1])
        field_ratio_str = iszero(total_field_energy[1]) ? "N/A" :
                          PlottingUtils.format_val((total_field_energy[i] -
                                                    total_field_energy[1]) /
                                                   total_field_energy[1])
        particle_lines = ["Particle",
                          "t = $(PlottingUtils.format_val(ti))",
                          "t/T = $(period_str)",
                          "z = $(PlottingUtils.format_val(z_arr[i]))",
                          "z_p = z_i/A = $(zp_str)",
                          "v = $(PlottingUtils.format_val(v_arr[i]))",
                          "v_p = v_i = $(PlottingUtils.format_val(v_arr[i]))",
                          "γ = $(PlottingUtils.format_val(γ_arr[i]))",
                          "E_p = $(PlottingUtils.format_val(particle_energy[i]))",
                          "ΔE_p = $(PlottingUtils.format_val(particle_energy[i] - particle_energy[1]))",
                          "ΔE_p / E_p(0) = $(particle_ratio_str)",
                          "-q Π_r(z)/γ = $(PlottingUtils.format_val(-q * Πrz_arr[i] / γ_arr[i]))",
                          "qγa∂ᵥΦ_s^FV(z) = $(PlottingUtils.format_val(singular_phi_source_rate[i]))",
                          "Particle RHS = $(PlottingUtils.format_val(total_particle_rhs_rate[i]))",
                          "m₀γ = $(PlottingUtils.format_val(m0 * γ_arr[i]))",
                          "-q Φ_r^FV(z)γ = $(PlottingUtils.format_val(-q * Φrz_fv_arr[i] * γ_arr[i]))",
                          "-q Φ_s^FV(z)γ = $(PlottingUtils.format_val(-q * Φsz_fv_arr[i] * γ_arr[i]))",
                          "(m₀ - q Φ_r^FV(z))γ = $(PlottingUtils.format_val((m0 - q * Φrz_fv_arr[i]) * γ_arr[i]))",
                          "(m₀ - q Φ_r^FV(z) - q Φ_s^FV(z))γ = $(PlottingUtils.format_val((m0 - q * Φrz_fv_arr[i] - q * Φsz_fv_arr[i]) * γ_arr[i]))",
                          "V(z) = $(PlottingUtils.format_val(potential_term[i]))",
                          "E_nonrad = $(PlottingUtils.format_val(non_radiative_energy[i]))"]
        field_lines = ["Field",
                       "E_r = $(PlottingUtils.format_val(regular_field_energy[i]))",
                       "E_× = $(PlottingUtils.format_val(cross_field_energy[i]))",
                       "E_punct (particle-field) = $(PlottingUtils.format_val(puncture_field_energy[i]))",
                       "E_field_full = $(PlottingUtils.format_val(total_field_energy[i]))",
                       "ΔE_field_full = $(PlottingUtils.format_val(total_field_energy[i] - total_field_energy[1]))",
                       "ΔE_field_full / E_field_full(0) = $(field_ratio_str)",
                       "E_field_full(0) = $(PlottingUtils.format_val(total_field_energy[1]))"]
        join(vcat(particle_lines, "", field_lines), "\n")
    end

    total_points = lift(i -> Point2f.(plot_x_data[i],
                                      plot_total_data[i] .- q^2 / 8),
                        current_index)
    regular_points = lift(i -> Point2f.(plot_x_data[i], plot_regular_data[i]),
                          current_index)
    cross_points = lift(i -> Point2f.(plot_x_data[i], plot_cross_data[i]), current_index)
    singular_points = lift(i -> Point2f.(plot_x_data[i], plot_singular_data[i]),
                           current_index)

    lines!(ax_left, total_points;
           label=L"\mathcal{E}_{\mathrm{tot}}", color=:black, linewidth=2.0)
    lines!(ax_left, regular_points;
           label=L"\mathcal{E}_r", color=colors[1], linewidth=1.6)
    lines!(ax_left, cross_points;
           label=L"\mathcal{E}_\times", color=colors[2], linewidth=1.6,
           linestyle=:dash)
    lines!(ax_left, singular_points;
           label=L"\mathcal{E}_s", color=colors[4], linewidth=1.6,
           linestyle=:dot)
    if !isnothing(boundary) && boundary < sim.sol.t[end]
        vlines!(ax_left, [-boundary, boundary];
                color=:firebrick, linestyle=:dash, linewidth=1.2)
    end
    on(lift((ti, er) -> (ti, er), current_time, current_energy_change_ratio)) do (ti, er)
        ax_left.title[] = "Field Energy Densities at t = $(round(ti, digits=6)) | ΔE_field / E_field(0) = $(round(er, digits=6))"
    end
    xlims!(ax_left, -tf, tf)
    ylims!(ax_left, left_density_limits...)

    lines!(ax_right, t_mid, d_field_energy_dt;
           label=L"\Delta E_{\mathrm{field}}/\Delta t",
           color=colors[1], linewidth=1.5)
    lines!(ax_right, t_mid, .-flux_mid;
           label=L"-\mathcal{F}_{\mathrm{out}}",
           color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax_right, t_mid, work_mid_plain;
           label=L"-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=colors[4], linestyle=:dashdot, linewidth=1.4)
    lines!(ax_right, t_mid, rhs_rate_plain;
           label=L"-\mathcal{F}_{\mathrm{out}}-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=:black, linewidth=1.2)
    lines!(ax_right, t_mid, residual_plain;
           label=L"\mathrm{Residual}\ (q)",
           color=:gray35, linestyle=:dash, linewidth=1.1)
    vlines!(ax_right, lift(ti -> [ti], current_time);
            color=:firebrick, linestyle=:dash, linewidth=1.2)
    xlims!(ax_right, first(source_ts), last(source_ts))

    lines!(ax_particle, source_ts, particle_energy .- 1;
           label=L"E_p(t)",
           color=colors[1], linewidth=1.5)
    lines!(ax_particle, source_ts, kinetic_term .- 1;
           label=L"m_0\gamma",
           color=colors[2], linewidth=1.4)
    lines!(ax_particle, source_ts, field_coupling_term;
           label=L"-q\Phi_r^{\mathrm{FV}}(z)\gamma",
           color=colors[3], linewidth=1.4)
    lines!(ax_particle, source_ts, singular_field_coupling_term;
           label=L"-q\Phi_s^{\mathrm{FV}}(z)\gamma",
           color=:gray35, linestyle=:dash, linewidth=1.4)
    lines!(ax_particle, source_ts, potential_term;
           label=L"V(z)",
           color=colors[4], linewidth=1.4)
    lines!(ax_particle, source_ts, non_radiative_energy .- 1;
           label=L"E_{\mathrm{nonrad}}(t)",
           color=colors[5], linewidth=1.4)
    vlines!(ax_particle, lift(ti -> [ti], current_time);
            color=:firebrick, linestyle=:dash, linewidth=1.2)
    xlims!(ax_particle, first(source_ts), last(source_ts))

    Label(diagnostics_layout[2, 1], diagnostics_text;
          halign=:left, valign=:top, justification=:left, tellwidth=false,
          tellheight=false, fontsize=26)

    Legend(fig[2, 1], ax_left; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, labelsize=22, colgap=20, rowgap=10, padding=(10, 10, 10, 10))
    Legend(fig[2, 2], ax_right; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=2, labelsize=22, colgap=20, rowgap=10,
           padding=(10, 10, 10, 10))
    Legend(fig[4, 1:2], ax_particle; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=3, labelsize=22, colgap=20, rowgap=10,
           padding=(10, 10, 10, 10))

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Auto())
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Auto())
    rowsize!(fig.layout, 5, Fixed(42))
    colsize!(fig.layout, 1, Relative(0.38))
    colsize!(fig.layout, 2, Relative(0.38))
    colsize!(fig.layout, 3, Relative(0.24))
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 10)

    rowsize!(diagnostics_layout, 1, Auto())
    rowsize!(diagnostics_layout, 2, Auto())
    colsize!(diagnostics_layout, 1, Relative(1.0))

    return fig
end

function plot_interactive_field_particle_energy_dashboard(sim; boundary=nothing)
    _require_coupled_field_particle_sim(sim,
                                        "plot_interactive_field_particle_energy_dashboard")
    source_ts = sim.sol.t
    length(source_ts) >= 2 ||
        throw(ArgumentError("Need at least 2 sampled times to build the interactive field-particle dashboard; got only $(length(source_ts))."))

    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)
    tf = source_ts[end]

    cutoff = isnothing(boundary) ? source_ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(source_ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    plot_mask = (x .>= -tf) .& (x .<= tf)
    if iszero(tf) || !any(plot_mask)
        plot_mask = trues(length(x))
    end
    plot_x_filtered = x[plot_mask]
    plot_balance_mask = (plot_x_filtered .>= -cutoff) .& (plot_x_filtered .<= cutoff)
    if iszero(source_ts[end]) || !any(plot_balance_mask)
        plot_balance_mask = trues(length(plot_x_filtered))
    end
    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)
    a = left_face_x
    b = right_face_x

    nt = length(source_ts)
    field_energy_total = zeros(nt)
    flux_out = zeros(nt)
    work_trace_plain = zeros(nt)
    regular_field_energy = zeros(nt)
    cross_field_energy = zeros(nt)
    z_arr = zeros(nt)
    v_arr = zeros(nt)
    γ_arr = zeros(nt)
    Πrz_arr = zeros(nt)
    Ψrz_arr = zeros(nt)
    Φrz_fv_arr = zeros(nt)
    particle_energy = zeros(nt)
    kinetic_term = zeros(nt)
    field_coupling_term = zeros(nt)
    singular_field_coupling_term = zeros(nt)
    singular_phi_source_rate = zeros(nt)
    total_particle_rhs_rate = zeros(nt)
    Φsz_fv_arr = zeros(nt)
    potential_term = zeros(nt)
    non_radiative_energy = zeros(nt)
    non_radiative_sol = ParticlePotentialPlots.non_radiative_reference_solution(sim)
    plot_x_data = Vector{Vector{Float64}}(undef, nt)
    plot_regular_data = Vector{Vector{Float64}}(undef, nt)
    plot_cross_data = Vector{Vector{Float64}}(undef, nt)
    plot_singular_data = Vector{Vector{Float64}}(undef, nt)
    plot_total_data = Vector{Vector{Float64}}(undef, nt)
    left_density_min = Inf
    left_density_max = -Inf

    for i in eachindex(source_ts)
        m = sim.sol.u[i].x[2][1]
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        z_arr[i] = z
        v_arr[i] = v
        Q_field = sim.sol.u[i].x[1]

        Πr_plot = @view Q_field[plot_mask, 1]
        Ψr_plot = @view Q_field[plot_mask, 2]
        Πs_plot = ScalarField.cell_averages_Πs(plot_x_filtered, h, z, v, q)
        Ψs_plot = ScalarField.cell_averages_Ψs(plot_x_filtered, h, z, v, q)
        regular_density_plot = @. 0.5 * (Πr_plot^2 + Ψr_plot^2)
        cross_density_plot = @. Πr_plot * Πs_plot + Ψr_plot * Ψs_plot
        singular_density_plot = @. 0.5 * (Πs_plot^2 + Ψs_plot^2) - (q^2 / 8.0)
        total_density_plot = regular_density_plot .+ cross_density_plot .+
                             singular_density_plot

        regular_field_energy[i] = 0.0
        cross_field_energy[i] = 0.0
        field_energy_total[i] = 0.0

        left_face_balance = x_filtered[1] - 0.5 * h
        right_face_balance = x_filtered[end] + 0.5 * h

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))
        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        flux_out[i] = -Πtot_right_face * Ψtot_right_face + Πtot_left_face * Ψtot_left_face

        Πrz_arr[i], Ψrz_arr[i] = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Φrz_fv_arr[i] = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        Φsz_fv_arr[i] = ODE.gather_muscl(ScalarField.cell_average_Φs, z, params, v, q)
        γ_arr[i] = γ
        kinetic_term[i] = m0 * γ
        field_coupling_term[i] = q * Φrz_fv_arr[i] * γ
        singular_field_coupling_term[i] = q * Φsz_fv_arr[i] * γ
        potential_term[i] = V(z)
        particle_energy[i] = kinetic_term[i] + field_coupling_term[i] +
                             singular_field_coupling_term[i] + potential_term[i]

        Πtrace_avg = zero(eltype(Πr_plot))
        if left_face_balance <= z <= right_face_balance
            balance_particle_cell = clamp(floor(Int, (z - left_face_balance) / h) + 1,
                                          1, length(x_filtered))
            balance_cell_center = x_filtered[balance_particle_cell]
            balance_cell_left_face = balance_cell_center - 0.5 * h
            balance_cell_right_face = balance_cell_center + 0.5 * h
            if balance_cell_left_face <= z <= balance_cell_right_face
                frac_left = (z - balance_cell_left_face) / h
                frac_right = one(frac_left) - frac_left
                Πs_left, Πs_right, _, _ = ScalarField.singular_particle_traces(z, v, q)
                global_cell = left_idx + balance_particle_cell - 1
                q_particle_minus, q_particle_plus = ODE.reconstruct_cell(Q_field,
                                                                         global_cell,
                                                                         params.N, L, R,
                                                                         params.slope_limiter)
                Πr_trace = frac_right * q_particle_minus[1] + frac_left * q_particle_plus[1]
                Πtrace_avg = 0.5 * ((Πr_trace + Πs_left) + (Πr_trace + Πs_right))
            end
        end
        work_trace_plain[i] = -q * Πtrace_avg

        z_nr = non_radiative_sol.u[i][2]
        v_nr = non_radiative_sol.u[i][3]
        γ_nr = inv(sqrt(1.0 - v_nr * v_nr))
        non_radiative_energy[i] = m0 * γ_nr + V(z_nr)

        plot_particle_cell = nothing
        plot_left_face_x = plot_x_filtered[1] - 0.5 * h
        plot_right_face_x = plot_x_filtered[end] + 0.5 * h
        if plot_left_face_x <= z <= plot_right_face_x
            plot_particle_cell = clamp(floor(Int, (z - plot_left_face_x) / h) + 1,
                                       1, length(plot_x_filtered))
            plot_cell_center = plot_x_filtered[plot_particle_cell]
            plot_cell_left_face = plot_cell_center - 0.5 * h
            plot_cell_right_face = plot_cell_center + 0.5 * h
            if !(plot_cell_left_face <= z <= plot_cell_right_face)
                plot_particle_cell = nothing
            end
        end
        Πs_left = 0.0
        Πs_right = 0.0
        Ψs_left = 0.0
        Ψs_right = 0.0
        if !isnothing(plot_particle_cell)
            Πs_left, Πs_right, Ψs_left, Ψs_right = ScalarField.singular_particle_traces(z,
                                                                                        v,
                                                                                        q)
        end

        x_plot = Float64[]
        regular_plot = Float64[]
        cross_plot = Float64[]
        singular_plot = Float64[]
        total_plot = Float64[]
        for j in eachindex(plot_x_filtered)
            cell_center = plot_x_filtered[j]
            cell_left_face = cell_center - 0.5 * h
            cell_right_face = cell_center + 0.5 * h

            if !isnothing(plot_particle_cell) && j == plot_particle_cell
                regular_value = regular_density_plot[j]
                cross_left = Πr_plot[j] * Πs_left + Ψr_plot[j] * Ψs_left
                cross_right = Πr_plot[j] * Πs_right + Ψr_plot[j] * Ψs_right
                singular_left = 0.5 * (Πs_left^2 + Ψs_left^2) - (q^2 / 8.0)
                singular_right = 0.5 * (Πs_right^2 + Ψs_right^2) - (q^2 / 8.0)
                total_left = regular_value + cross_left + singular_left
                total_right = regular_value + cross_right + singular_right
                left_width = z - cell_left_face
                right_width = cell_right_face - z

                if plot_balance_mask[j]
                    regular_field_energy[i] += h * regular_value
                    cross_field_energy[i] += left_width * cross_left +
                                             right_width * cross_right
                    field_energy_total[i] += left_width * total_left +
                                             right_width * total_right
                end

                append!(x_plot, (cell_left_face, z, z, cell_right_face))
                append!(regular_plot,
                        (regular_value, regular_value, regular_value, regular_value))
                append!(cross_plot, (cross_left, cross_left, cross_right, cross_right))
                append!(singular_plot,
                        (singular_left, singular_left, singular_right, singular_right))
                append!(total_plot,
                        (total_left, total_left, total_right, total_right))
            else
                if plot_balance_mask[j]
                    regular_field_energy[i] += h * regular_density_plot[j]
                    cross_field_energy[i] += h * cross_density_plot[j]
                    field_energy_total[i] += h * total_density_plot[j]
                end
                append!(x_plot, (cell_left_face, cell_right_face))
                append!(regular_plot, (regular_density_plot[j], regular_density_plot[j]))
                append!(cross_plot, (cross_density_plot[j], cross_density_plot[j]))
                append!(singular_plot,
                        (singular_density_plot[j], singular_density_plot[j]))
                append!(total_plot, (total_density_plot[j], total_density_plot[j]))
            end
        end

        plot_x_data[i] = x_plot
        plot_regular_data[i] = regular_plot
        plot_cross_data[i] = cross_plot
        plot_singular_data[i] = singular_plot
        plot_total_data[i] = total_plot

        snapshot_min = min(minimum(regular_plot),
                           minimum(cross_plot),
                           minimum(singular_plot),
                           minimum(total_plot))
        snapshot_max = max(maximum(regular_plot),
                           maximum(cross_plot),
                           maximum(singular_plot),
                           maximum(total_plot))
        left_density_min = min(left_density_min, snapshot_min)
        left_density_max = max(left_density_max, snapshot_max)
    end

    a_arr = zeros(nt)
    if nt >= 2
        a_arr[1] = (v_arr[2] - v_arr[1]) / (source_ts[2] - source_ts[1])
        a_arr[end] = (v_arr[end] - v_arr[end - 1]) / (source_ts[end] - source_ts[end - 1])
    end
    for i in 2:(nt - 1)
        a_arr[i] = (v_arr[i + 1] - v_arr[i - 1]) / (source_ts[i + 1] - source_ts[i - 1])
    end
    for i in eachindex(source_ts)
        z = z_arr[i]
        v = v_arr[i]
        γ = γ_arr[i]
        a = a_arr[i]
        dvΦsz_fv = ODE.gather_muscl(ScalarField.cell_average_dvΦs, z, params, v, q)
        singular_phi_source_rate[i] = q * γ * a * dvΦsz_fv
        total_particle_rhs_rate[i] = q * Πrz_arr[i] / γ + singular_phi_source_rate[i]
    end

    dt_steps = diff(source_ts)
    t_mid = @. 0.5 * (source_ts[1:(end - 1)] + source_ts[2:end])
    d_field_energy_dt = diff(field_energy_total) ./ dt_steps
    flux_mid = @. 0.5 * (flux_out[1:(end - 1)] + flux_out[2:end])
    work_mid_plain = @. 0.5 * (work_trace_plain[1:(end - 1)] + work_trace_plain[2:end])
    rhs_rate_plain = .-flux_mid .+ work_mid_plain
    residual_plain = d_field_energy_dt .- rhs_rate_plain

    y_pad = left_density_max > left_density_min ?
            0.05 * (left_density_max - left_density_min) :
            max(abs(left_density_max), 1.0) * 0.05
    left_density_limits = (left_density_min - y_pad, left_density_max + y_pad)

    activate_glmakie!()
    colors = [:dodgerblue3, :darkorange2, :seagreen4, :firebrick3, :mediumpurple3]

    fig = Figure(; size=(1260, 860), figure_padding=(10, 6, 10, 10))

    ax_left = Axis(fig[1, 1];
                   title="Field Energy Densities",
                   xlabel=L"x",
                   ylabel=L"\mathcal{E}(x,t)",
                   titlesize=30,
                   xlabelsize=26,
                   ylabelsize=26,
                   xticklabelsize=22,
                   yticklabelsize=22)
    ax_right = Axis(fig[1, 2];
                    title=L"\textrm{Rate Balance: Total Field}",
                    xlabel=L"t",
                    ylabel=L"\mathrm{Rate}",
                    titlesize=30,
                    xlabelsize=26,
                    ylabelsize=26,
                    xticklabelsize=22,
                    yticklabelsize=22)
    ax_particle = Axis(fig[3, 1:2];
                       title=L"\textrm{Particle Energy}",
                       xlabel=L"t",
                       ylabel=L"E_p",
                       titlesize=30,
                       xlabelsize=26,
                       ylabelsize=26,
                       xticklabelsize=22,
                       yticklabelsize=22)

    diagnostics_layout = fig[1:4, 3] = GridLayout()
    Label(diagnostics_layout[1, 1], "Diagnostics";
          fontsize=30, halign=:left, tellwidth=false)

    slider_layout = fig[5, 1:2] = GridLayout()
    Label(slider_layout[1, 1], "Time index"; fontsize=24)
    time_slider = Slider(slider_layout[1, 2];
                         range=1:length(source_ts),
                         startvalue=1,
                         horizontal=true)
    time_label = lift(time_slider.value) do i
        "t = $(round(source_ts[i], digits=6))"
    end
    Label(slider_layout[1, 3], time_label; fontsize=24)

    current_index = time_slider.value
    current_time = lift(current_index) do i
        source_ts[i]
    end
    current_energy_change_ratio = lift(current_index) do i
        (field_energy_total[i] - field_energy_total[1]) / field_energy_total[1]
    end
    diagnostics_text = lift(current_index, current_time) do i, ti
        period_str = iszero(sim.params.ω) ? "N/A" :
                     PlottingUtils.format_val(ti / (2π / sim.params.ω))
        zp_str = iszero(sim.params.A) ? "N/A" :
                 PlottingUtils.format_val(z_arr[i] / sim.params.A)
        particle_ratio_str = iszero(particle_energy[1]) ? "N/A" :
                             PlottingUtils.format_val((particle_energy[i] -
                                                       particle_energy[1]) /
                                                      particle_energy[1])
        field_ratio_str = iszero(field_energy_total[1]) ? "N/A" :
                          PlottingUtils.format_val((field_energy_total[i] -
                                                    field_energy_total[1]) /
                                                   field_energy_total[1])
        particle_lines = ["Particle",
                          "t = $(PlottingUtils.format_val(ti))",
                          "t/T = $(period_str)",
                          "z = $(PlottingUtils.format_val(z_arr[i]))",
                          "z_p = z_i/A = $(zp_str)",
                          "v = $(PlottingUtils.format_val(v_arr[i]))",
                          "v_p = v_i = $(PlottingUtils.format_val(v_arr[i]))",
                          "γ = $(PlottingUtils.format_val(γ_arr[i]))",
                          "E_p = $(PlottingUtils.format_val(particle_energy[i]))",
                          "ΔE_p = $(PlottingUtils.format_val(particle_energy[i] - particle_energy[1]))",
                          "ΔE_p / E_p(0) = $(particle_ratio_str)",
                          "q Π_r(z)/γ = $(PlottingUtils.format_val(q * Πrz_arr[i] / γ_arr[i]))",
                          "qγa∂ᵥΦ_s^FV(z) = $(PlottingUtils.format_val(singular_phi_source_rate[i]))",
                          "Particle RHS = $(PlottingUtils.format_val(total_particle_rhs_rate[i]))",
                          "mγ = $(PlottingUtils.format_val(kinetic_term[i]))",
                          "q Φ_r^FV(z)γ = $(PlottingUtils.format_val(field_coupling_term[i]))",
                          "q Φ_s^FV(z)γ = $(PlottingUtils.format_val(singular_field_coupling_term[i]))",
                          "V(z) = $(PlottingUtils.format_val(potential_term[i]))",
                          "E_nonrad = $(PlottingUtils.format_val(non_radiative_energy[i]))"]
        field_lines = ["Field",
                       "E_r = $(PlottingUtils.format_val(regular_field_energy[i]))",
                       "E_× = $(PlottingUtils.format_val(cross_field_energy[i]))",
                       "E_field_full = $(PlottingUtils.format_val(field_energy_total[i]))",
                       "ΔE_field_full = $(PlottingUtils.format_val(field_energy_total[i] - field_energy_total[1]))",
                       "ΔE_field_full / E_field_full(0) = $(field_ratio_str)",
                       "E_field_full(0) = $(PlottingUtils.format_val(field_energy_total[1]))"]
        join(vcat(particle_lines, "", field_lines), "\n")
    end

    total_points = lift(i -> Point2f.(plot_x_data[i], plot_total_data[i]), current_index)
    regular_points = lift(i -> Point2f.(plot_x_data[i], plot_regular_data[i]),
                          current_index)
    cross_points = lift(i -> Point2f.(plot_x_data[i], plot_cross_data[i]), current_index)
    singular_points = lift(i -> Point2f.(plot_x_data[i], plot_singular_data[i]),
                           current_index)

    lines!(ax_left, total_points;
           label=L"\mathcal{E}_{\mathrm{tot}}", color=:black, linewidth=2.0)
    lines!(ax_left, regular_points;
           label=L"\mathcal{E}_r", color=colors[1], linewidth=1.6)
    lines!(ax_left, cross_points;
           label=L"\mathcal{E}_\times", color=colors[2], linewidth=1.6,
           linestyle=:dash)
    lines!(ax_left, singular_points;
           label=L"\mathcal{E}_s", color=colors[4], linewidth=1.6,
           linestyle=:dot)
    if !isnothing(boundary) && boundary < sim.sol.t[end]
        vlines!(ax_left, [-boundary, boundary];
                color=:firebrick, linestyle=:dash, linewidth=1.2)
    end
    on(lift((ti, er) -> (ti, er), current_time, current_energy_change_ratio)) do (ti, er)
        ax_left.title[] = "Field Energy Densities at t = $(round(ti, digits=6)) | ΔE_field / E_field(0) = $(round(er, digits=6))"
    end
    xlims!(ax_left, -tf, tf)
    ylims!(ax_left, left_density_limits...)

    lines!(ax_right, t_mid, d_field_energy_dt;
           label=L"\Delta E_{\mathrm{field}}/\Delta t",
           color=colors[1], linewidth=1.5)
    lines!(ax_right, t_mid, .-flux_mid;
           label=L"-\mathcal{F}_{\mathrm{out}}",
           color=colors[2], linestyle=:dash, linewidth=1.4)
    lines!(ax_right, t_mid, work_mid_plain;
           label=L"-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=colors[4], linestyle=:dashdot, linewidth=1.4)
    lines!(ax_right, t_mid, rhs_rate_plain;
           label=L"-\mathcal{F}_{\mathrm{out}}-q\,\overline{\Pi}_{\mathrm{tot}}^{\mathrm{tr}}(z)",
           color=:black, linewidth=1.2)
    lines!(ax_right, t_mid, residual_plain;
           label=L"\mathrm{Residual}\ (q)",
           color=:gray35, linestyle=:dash, linewidth=1.1)
    vlines!(ax_right, lift(ti -> [ti], current_time);
            color=:firebrick, linestyle=:dash, linewidth=1.2)
    xlims!(ax_right, first(source_ts), last(source_ts))

    lines!(ax_particle, source_ts, particle_energy .- 1;
           label=L"E_p(t)",
           color=colors[1], linewidth=1.5)
    lines!(ax_particle, source_ts, kinetic_term .- 1;
           label=L"m\gamma",
           color=colors[2], linewidth=1.4)
    lines!(ax_particle, source_ts, field_coupling_term;
           label=L"q\Phi_r^{\mathrm{FV}}(z)\gamma",
           color=colors[3], linewidth=1.4)
    lines!(ax_particle, source_ts, singular_field_coupling_term;
           label=L"q\Phi_s^{\mathrm{FV}}(z)\gamma",
           color=:gray35, linestyle=:dash, linewidth=1.4)
    lines!(ax_particle, source_ts, potential_term;
           label=L"V(z)",
           color=colors[4], linewidth=1.4)
    lines!(ax_particle, source_ts, non_radiative_energy .- 1;
           label=L"E_{\mathrm{nonrad}}(t)",
           color=colors[5], linewidth=1.4)
    vlines!(ax_particle, lift(ti -> [ti], current_time);
            color=:firebrick, linestyle=:dash, linewidth=1.2)
    xlims!(ax_particle, first(source_ts), last(source_ts))

    Label(diagnostics_layout[2, 1], diagnostics_text;
          halign=:left, valign=:top, justification=:left, tellwidth=false,
          tellheight=false, fontsize=26)

    Legend(fig[2, 1], ax_left; orientation=:horizontal, tellwidth=false, tellheight=true,
           nbanks=2, labelsize=22, colgap=20, rowgap=10, padding=(10, 10, 10, 10))
    Legend(fig[2, 2], ax_right; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=2, labelsize=22, colgap=20, rowgap=10,
           padding=(10, 10, 10, 10))
    Legend(fig[4, 1:2], ax_particle; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=3, labelsize=22, colgap=20, rowgap=10,
           padding=(10, 10, 10, 10))

    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Auto())
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, Auto())
    rowsize!(fig.layout, 5, Fixed(42))
    colsize!(fig.layout, 1, Relative(0.38))
    colsize!(fig.layout, 2, Relative(0.38))
    colsize!(fig.layout, 3, Relative(0.24))
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 10)

    rowsize!(diagnostics_layout, 1, Auto())
    rowsize!(diagnostics_layout, 2, Auto())
    colsize!(diagnostics_layout, 1, Relative(1.0))

    return fig
end

# Plot the exact field, particle, and total rate balances using the total-field work term.
function plot_exact_total_energy_balance(sim; paper=false, boundary=nothing)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    h = sim.params.dx
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    params = sim.sol.prob.p
    V(x) = InitialData.Harmonic(x, k)

    # Build the finite control volume Ω and the geometry of its two boundary faces.
    cutoff = isnothing(boundary) ? ts[end] : boundary
    mask = (x .>= -cutoff) .& (x .<= cutoff)
    if iszero(ts[end]) || !any(mask)
        mask = trues(length(x))
    end
    x_filtered = x[mask]
    left_idx = findfirst(mask)
    right_idx = findlast(mask)
    left_face_x = x[left_idx] - 0.5 * h
    right_face_x = x[right_idx] + 0.5 * h

    L = Equations.left_eigenmatrix(params.equation)
    R = Equations.right_eigenmatrix(params.equation)

    field_energy_total = zeros(nt)
    particle_energy_total = zeros(nt)
    total_energy = zeros(nt)
    flux_total_values = zeros(nt)
    field_work_total = zeros(nt)
    particle_work_total = zeros(nt)

    # Compute the exact total-field energy, the particle energy paired with the
    # total field, and the three rate laws shown in the figure.
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Q_field = sim.sol.u[i].x[1]
        Πr = @view Q_field[mask, 1]
        Ψr = @view Q_field[mask, 2]
        Φr = @view Q_field[mask, 3]

        Πs = zeros(length(x_filtered))
        Ψs = zeros(length(x_filtered))
        @inbounds for (idx, x_i) in pairs(x_filtered)
            Πs[idx], Ψs[idx] = ScalarField.cell_average_Πs(x_i, h, z, v, q),
                               ScalarField.cell_average_Ψs(x_i, h, z, v, q)
        end

        total_density = @. 0.5 * (Πr^2 + Ψr^2) + Πr * Πs + Ψr * Ψs + 0.5 * (Πs^2 + Ψs^2)
        field_energy_total[i] = integrate_over_domain(total_density, h)

        q_left_minus, _ = ODE.reconstruct_cell(Q_field, left_idx, params.N, L, R,
                                               params.slope_limiter)
        _, q_right_plus = ODE.reconstruct_cell(Q_field, right_idx, params.N, L, R,
                                               params.slope_limiter)
        Πs_left_face, Ψs_left_face = ScalarField.singular_face_trace(left_face_x,
                                                                     z,
                                                                     v,
                                                                     q,
                                                                     one(v))
        Πs_right_face, Ψs_right_face = ScalarField.singular_face_trace(right_face_x,
                                                                       z,
                                                                       v,
                                                                       q,
                                                                       -one(v))
        Πtot_left_face = q_left_minus[1] + Πs_left_face
        Ψtot_left_face = q_left_minus[2] + Ψs_left_face
        Πtot_right_face = q_right_plus[1] + Πs_right_face
        Ψtot_right_face = q_right_plus[2] + Ψs_right_face
        flux_total_values[i] = -Πtot_right_face * Ψtot_right_face +
                               Πtot_left_face * Ψtot_left_face

        Φrz = ODE.gather_muscl(@view(Q_field[:, 3]), z, params)
        Πrz, _ = ODE.gather_muscl(@view(Q_field[:, 1:2]), z, params)
        Πsz, _ = ODE.gather_muscl(hcat(ScalarField.cell_averages_Πs(x, h, z, v, q),
                                       ScalarField.cell_averages_Ψs(x, h, z, v, q)),
                                  z, params)
        Φsz = ODE.gather_muscl(ScalarField.cell_averages_Φs(x, h, z, v, q), z, params)
        Φtotz = Φrz + Φsz
        Πtotz = Πrz + Πsz

        # Only attach the particle energy while the particle lies inside Ω,
        # matching the earlier exact-total diagnostic.
        if first(x_filtered) <= z <= last(x_filtered)
            particle_energy_total[i] = γ * (m0 - q * Φtotz) + V(z)
        end
        total_energy[i] = field_energy_total[i] + particle_energy_total[i]

        field_work_total[i] = q * Πtotz / γ
        particle_work_total[i] = -field_work_total[i]
    end

    # Differentiate the field, particle, and total energies and compare them to
    # the corresponding exact RHS rates.
    d_field_energy_dt = time_derivative(field_energy_total, ts)
    d_particle_energy_dt = time_derivative(particle_energy_total, ts)
    d_total_energy_dt = time_derivative(total_energy, ts)

    rhs_field_rate = .-flux_total_values .+ field_work_total
    rhs_particle_rate = particle_work_total
    rhs_total_rate = .-flux_total_values

    field_rate_residual = d_field_energy_dt .- rhs_field_rate
    particle_rate_residual = d_particle_energy_dt .- rhs_particle_rate
    total_rate_residual = d_total_energy_dt .- rhs_total_rate

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(560, 620), figure_padding=(6, 8, 8, 6))
    colors = mytheme.palette.color[]

    panel11 = GridLayout(fig[1, 1])
    panel12 = GridLayout(fig[1, 2])
    panel21 = GridLayout(fig[3, 1])
    panel22 = GridLayout(fig[3, 2])

    ax_energy = Axis(panel11[1, 1];
                     title=L"\textrm{Exact Total Energies}",
                     ylabel=L"E")

    # Show the exact field energy, particle energy, and their sum.
    lines!(ax_energy, ts, field_energy_total;
           label=L"E_{\mathrm{field}}",
           color=colors[1], linewidth=1.5)
    lines!(ax_energy, ts, particle_energy_total;
           label=L"E_{\mathrm{part}}",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax_energy, ts, total_energy;
           label=L"E_{\mathrm{field}} + E_{\mathrm{part}}",
           color=colors[4], linewidth=1.5)

    ax_field = Axis(panel12[1, 1];
                    title=L"\textrm{Field Balance}",
                    ylabel=L"\mathrm{Rate}")

    # Compare the field-energy derivative against its exact flux-plus-work law.
    lines!(ax_field, ts, d_field_energy_dt;
           label=L"\frac{dE_{\mathrm{field}}}{dt}",
           color=colors[1], linewidth=1.5)
    lines!(ax_field, ts, rhs_field_rate;
           label=L"-\mathcal{F}_{\mathrm{tot}} + q\Pi_{\mathrm{tot}}(z)/\gamma",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax_field, ts, field_rate_residual;
           label=L"\mathrm{Residual}",
           color=:black, linestyle=:dot, linewidth=1.2)

    ax_particle = Axis(panel21[1, 1];
                       title=L"\textrm{Particle Balance}",
                       xlabel=L"t",
                       ylabel=L"\mathrm{Rate}")

    # Compare the particle-energy derivative against the opposite worldline work.
    lines!(ax_particle, ts, d_particle_energy_dt;
           label=L"\frac{dE_{\mathrm{part}}}{dt}",
           color=colors[1], linewidth=1.5)
    lines!(ax_particle, ts, rhs_particle_rate;
           label=L"-q\Pi_{\mathrm{tot}}(z)/\gamma",
           color=colors[3], linestyle=:dash, linewidth=1.5)
    lines!(ax_particle, ts, particle_rate_residual;
           label=L"\mathrm{Residual}",
           color=:black, linestyle=:dot, linewidth=1.2)

    ax_total = Axis(panel22[1, 1];
                    title=L"\textrm{Total Balance}",
                    xlabel=L"t",
                    ylabel=L"\mathrm{Rate}")

    # Compare the exact total-energy derivative against the pure boundary flux.
    lines!(ax_total, ts, d_total_energy_dt;
           label=L"\frac{d(E_{\mathrm{field}}+E_{\mathrm{part}})}{dt}",
           color=colors[1], linewidth=1.5)
    lines!(ax_total, ts, rhs_total_rate;
           label=L"-\mathcal{F}_{\mathrm{tot}}",
           color=colors[2], linestyle=:dash, linewidth=1.5)
    lines!(ax_total, ts, total_rate_residual;
           label=L"\mathrm{Residual}",
           color=:black, linestyle=:dot, linewidth=1.2)

    apply_scientific_notation_if_needed!(ax_energy)
    apply_scientific_notation_if_needed!(ax_field)
    apply_scientific_notation_if_needed!(ax_particle)
    apply_scientific_notation_if_needed!(ax_total)

    Legend(panel11[2, 1], ax_energy; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=2)
    Legend(panel12[2, 1], ax_field; orientation=:horizontal, tellwidth=false,
           tellheight=true, nbanks=1)
    Legend(panel21[2, 1], ax_particle; orientation=:horizontal, tellwidth=false,
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

    savename = simple_name("exact_total_energy_balance")
    save(joinpath(get_fig_path(), savename), fig)
    if paper
        save(joinpath(get_paper_fig_path(), savename), fig)
    end
    return fig
end

# Generate the full set of one-series ingredient plots plus the balance-summary figures.
function prepare_energy_balance_ingredient_plots(sim; paper=false, boundary=nothing)
    plot_regular_field_energy(sim; paper=paper, boundary=boundary)
    plot_cross_field_energy(sim; paper=paper, boundary=boundary)
    plot_singular_field_energy(sim; paper=paper, boundary=boundary)
    plot_total_field_energy(sim; paper=paper, boundary=boundary)
    plot_regular_cross_field_energy(sim; paper=paper, boundary=boundary)
    plot_particle_full_field_energy_change_decomposition(sim; paper=paper,
                                                         boundary=boundary)
    plot_potential_energy(sim; paper=paper, boundary=boundary)
    plot_effective_mass(sim; paper=paper, boundary=boundary)
    plot_dressed_mass_energy(sim; paper=paper, boundary=boundary)
    plot_particle_energy(sim; paper=paper, boundary=boundary)
    plot_total_particle_energy_candidate(sim; paper=paper, boundary=boundary)
    plot_conserved_total_energy_candidate(sim; paper=paper, boundary=boundary)
    plot_flux_regular(sim; paper=paper, boundary=boundary)
    plot_flux_total(sim; paper=paper, boundary=boundary)
    plot_flux_regular_cross(sim; paper=paper, boundary=boundary)
    plot_work_singular_terms(sim; paper=paper, boundary=boundary)
    plot_work_total_field(sim; paper=paper, boundary=boundary)
    plot_energy_derivatives(sim; paper=paper, boundary=boundary)
    plot_particle_regular_rate_balance(sim; paper=paper, boundary=boundary)
    plot_total_field_rate_balance(sim; paper=paper, boundary=boundary)
    plot_particle_regular_integral_balance(sim; paper=paper, boundary=boundary)
    plot_conserved_total_energy_balance(sim; paper=paper, boundary=boundary)
    plot_total_energy_derivative_and_work(sim; paper=paper, boundary=boundary)
    plot_total_energy_integrated_work_balance(sim; paper=paper, boundary=boundary)
    plot_regular_cross_balance_law(sim; paper=paper, boundary=boundary)
    plot_regular_cross_integrated_residual(sim; paper=paper, boundary=boundary)
    plot_exact_total_energy_balance(sim; paper=paper, boundary=boundary)
    return nothing
end

end # module
