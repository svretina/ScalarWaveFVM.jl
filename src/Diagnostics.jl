module Diagnostics

using ..Equations
using ..ODE
using ..ScalarField
using StaticArrays: SVector

export semidiscrete_energy_terms, particle_observables, runtime_health,
       remainder_energy_balance

@inline _energy_dot(q, p, c) = q[1] * p[1] / c^2 + q[2] * p[2]
@inline _energy_flux(q, equation, c) =
    0.5 * _energy_dot(q, Equations.principal_part(equation) * q, c)

@inline function _face_flux(Q, i, params, L, R)
    _, q_left = ODE.reconstruct_cell(Q, i, params.N, L, R,
                                     params.slope_limiter)
    q_right, _ = ODE.reconstruct_cell(Q, i + 1, params.N, L, R,
                                      params.slope_limiter)
    return params.numerical_flux(q_left, q_right, params.equation)
end

"""
    semidiscrete_energy_terms(Q, source, params)

Return the exact periodic semidiscrete energy decomposition for the MUSCL
operator used by `ODE.muscl_optimized!`:

`rate = source_power - dissipation + reconstruction_residual`.

The energy is `h/2 Σ(Π²/c² + Ψ²)`. Godunov and Rusanov both have wave-speed
dissipation `c/2 ΔQᵀHΔQ`; the residual measures the exact change caused by
feeding reconstructed rather than cell-average states to the numerical flux.
"""
function semidiscrete_energy_terms(Q, source, params)
    N, h = params.N, params.h
    equation = params.equation
    c = equation.velocity
    flux = params.numerical_flux
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)
    limiter = params.slope_limiter

    source_power = zero(eltype(Q))
    dissipation = zero(eltype(Q))
    reconstruction_residual = zero(eltype(Q))

    for i in 1:N
        ip1 = mod(i, N) + 1
        Qi = SVector(Q[i, 1], Q[i, 2])
        Qip1 = SVector(Q[ip1, 1], Q[ip1, 2])
        _, qL = ODE.reconstruct_cell(Q, i, N, L, R, limiter)
        qR, _ = ODE.reconstruct_cell(Q, i + 1, N, L, R, limiter)
        actual_flux = flux(qL, qR, equation)
        ΔQ = Qip1 - Qi
        cell_flux = Equations.principal_part(equation) * (Qi + Qip1) / 2 - c * ΔQ / 2

        dissipation += 0.5 * c * _energy_dot(ΔQ, ΔQ, c)
        reconstruction_residual += _energy_dot(ΔQ, actual_flux - cell_flux, c)
        source_power += h * _energy_dot(Qi,
                                        SVector(source[i, 1], source[i, 2]), c)
    end

    rate = source_power - dissipation + reconstruction_residual
    return (; rate, source_power, dissipation, reconstruction_residual,
            boundary_flux=zero(rate))
end

"""
    semidiscrete_energy_terms(Q, source, params, cells)

Return the exact semidiscrete energy decomposition on a contiguous, interior
range of finite-volume cells. `boundary_flux_left` and `boundary_flux_right`
are positive outwards from the selected control volume. With
`g(Q)=QᵀHAQ/2` and numerical face flux `F̂`, they are
`g(Q_a)-Q_aᵀHF̂_{a-1/2}` and `Q_bᵀHF̂_{b+1/2}-g(Q_b)`. This gives
`rate = source_power - boundary_flux - dissipation + reconstruction_residual`.
"""
function semidiscrete_energy_terms(Q, source, params,
                                   cells::AbstractUnitRange{<:Integer})
    first(cells) > 1 && last(cells) < params.N ||
        throw(ArgumentError("the control volume must be interior to the grid"))

    h = params.h
    equation = params.equation
    c = equation.velocity
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)
    source_power = zero(eltype(Q))
    dissipation = zero(eltype(Q))
    reconstruction_residual = zero(eltype(Q))

    for i in cells
        Qi = SVector(Q[i, 1], Q[i, 2])
        source_power += h * _energy_dot(Qi,
                                        SVector(source[i, 1], source[i, 2]), c)
    end

    for i in first(cells):(last(cells) - 1)
        Qi = SVector(Q[i, 1], Q[i, 2])
        Qip1 = SVector(Q[i + 1, 1], Q[i + 1, 2])
        actual_flux = _face_flux(Q, i, params, L, R)
        ΔQ = Qip1 - Qi
        cell_flux = Equations.principal_part(equation) * (Qi + Qip1) / 2 -
                    c * ΔQ / 2
        dissipation += 0.5 * c * _energy_dot(ΔQ, ΔQ, c)
        reconstruction_residual += _energy_dot(ΔQ,
                                                actual_flux - cell_flux, c)
    end

    first_cell, last_cell = first(cells), last(cells)
    Qleft = SVector(Q[first_cell, 1], Q[first_cell, 2])
    Qright = SVector(Q[last_cell, 1], Q[last_cell, 2])
    flux_left = _face_flux(Q, first_cell - 1, params, L, R)
    flux_right = _face_flux(Q, last_cell, params, L, R)
    boundary_flux_left = _energy_flux(Qleft, equation, c) -
                         _energy_dot(Qleft, flux_left, c)
    boundary_flux_right = _energy_dot(Qright, flux_right, c) -
                          _energy_flux(Qright, equation, c)
    boundary_flux = boundary_flux_left + boundary_flux_right
    rate = source_power - boundary_flux - dissipation +
           reconstruction_residual
    return (; rate, source_power, dissipation, reconstruction_residual,
            boundary_flux, boundary_flux_left, boundary_flux_right)
end

function particle_observables(Q_field, particle, params; m0=params.m0, τ0=0.0)
    m, z, v, τ = particle
    Φr, Πr, Ψr = ODE.regular_fields_at_particle(Q_field, z, params)
    reconstructed_mass = ODE.dynamical_mass(m0, params.q, Φr, τ; τ0=τ0)
    γ = inv(sqrt(1 - v^2))
    return (; m, z, v, τ, γ, Φr, Πr, Ψr, reconstructed_mass,
            mass_constraint=m - reconstructed_mass,
            scalar_force_combination=Ψr + v * Πr)
end

function runtime_health(states)
    masses = (u.x[2][1] for u in states)
    velocities = (u.x[2][3] for u in states)
    minimum_mass = minimum(masses)
    maximum_speed = maximum(abs, velocities)
    finite = all(u -> all(isfinite, u.x[1]) && all(isfinite, u.x[2]), states)
    return (; minimum_mass, maximum_speed, finite,
            valid=finite && minimum_mass > 0 && maximum_speed < 1)
end

@inline function _cumulative_trapezoid(t, y)
    out = zeros(promote_type(eltype(t), eltype(y)), length(t))
    for i in 2:length(t)
        out[i] = out[i - 1] + 0.5 * (t[i] - t[i - 1]) * (y[i - 1] + y[i])
    end
    return out
end

"""Evaluate the semidiscrete remainder-energy balance over a solution."""
function remainder_energy_balance(sim; domain=nothing)
    t = sim.sol.t
    params = sim.sol.prob.p
    h = params.h
    c = params.equation.velocity
    nt = length(t)
    energy = zeros(nt)
    source_power = zeros(nt)
    boundary_flux = zeros(nt)
    boundary_flux_left = zeros(nt)
    boundary_flux_right = zeros(nt)
    dissipation = zeros(nt)
    reconstruction = zeros(nt)
    energy_rate = zeros(nt)
    source = zeros(eltype(first(sim.sol.u).x[1]), size(first(sim.sol.u).x[1]))
    dQ = similar(source)

    if domain === nothing
        cells = 1:params.N
        actual_domain = nothing
    else
        xmin, xmax = domain
        xmin < xmax || throw(ArgumentError("domain must satisfy xmin < xmax"))
        indices = findall(x -> xmin <= x <= xmax, params.x)
        isempty(indices) && throw(ArgumentError("domain contains no grid cells"))
        cells = first(indices):last(indices)
        actual_domain = (params.x[first(cells)] - 0.5h,
                         params.x[last(cells)] + 0.5h)
    end

    for j in eachindex(t)
        Q = sim.sol.u[j].x[1]
        particle = sim.sol.u[j].x[2]
        m, z, v = particle[1], particle[2], particle[3]
        _, Πr, Ψr = ODE.regular_fields_at_particle(Q, z, params)
        feedback = hasproperty(params, :particle_feedback) ? params.particle_feedback : true
        rates = ODE.particle_rates(m, z, v, params.q, Πr, Ψr,
                                   params.dxPotential;
                                   particle_feedback=feedback)
        fill!(source, zero(eltype(source)))
        if !hasproperty(params, :field_source) || params.field_source
            for i in 1:params.N
                ODE.add_singular_sources!(source, i, params.x[i], h, z, v,
                                          params.q, rates.dv)
            end
        end
        terms = domain === nothing ? semidiscrete_energy_terms(Q, source, params) :
                semidiscrete_energy_terms(Q, source, params, cells)
        ODE.muscl_optimized!(dQ, Q, params, t[j])
        dQ .+= source
        energy[j] = 0.5h * sum(@view(Q[cells, 1]) .^ 2 ./ c^2 .+
                               @view(Q[cells, 2]) .^ 2)
        energy_rate[j] = h * sum(@view(Q[cells, 1]) .* @view(dQ[cells, 1]) ./ c^2 .+
                                 @view(Q[cells, 2]) .* @view(dQ[cells, 2]))
        source_power[j] = terms.source_power
        boundary_flux[j] = terms.boundary_flux
        if domain !== nothing
            boundary_flux_left[j] = terms.boundary_flux_left
            boundary_flux_right[j] = terms.boundary_flux_right
        end
        dissipation[j] = terms.dissipation
        reconstruction[j] = terms.reconstruction_residual
    end

    integrated_source = _cumulative_trapezoid(t, source_power)
    integrated_boundary_flux = _cumulative_trapezoid(t, boundary_flux)
    integrated_boundary_flux_left = _cumulative_trapezoid(t, boundary_flux_left)
    integrated_boundary_flux_right = _cumulative_trapezoid(t, boundary_flux_right)
    integrated_dissipation = _cumulative_trapezoid(t, dissipation)
    integrated_reconstruction = _cumulative_trapezoid(t, reconstruction)
    energy_change = energy .- first(energy)
    predicted_change = integrated_source .- integrated_boundary_flux .-
                       integrated_dissipation .+ integrated_reconstruction
    residual = energy_change .- predicted_change
    predicted_rate = source_power .- boundary_flux .- dissipation .+
                     reconstruction
    power_residual = energy_rate .- predicted_rate
    return (; t, domain=actual_domain, cells, energy, energy_change, source_power,
            boundary_flux, boundary_flux_left, boundary_flux_right, dissipation,
            reconstruction, energy_rate, predicted_rate, power_residual,
            integrated_source,
            integrated_boundary_flux,
            integrated_boundary_flux_left, integrated_boundary_flux_right,
            integrated_dissipation, integrated_reconstruction, predicted_change,
            residual)
end

end
