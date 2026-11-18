module Energies

export newton_cotes, field_energy_density, field_momentum_density,
       field_energy, field_momentum, particle_energy, kinetic_energy

function newton_cotes(y::AbstractVector{T}, h::T; order::Int=1) where {T}
    n = length(y) - 1
    n >= order || throw(ArgumentError("insufficient points for order $order"))
    n % order == 0 || throw(ArgumentError("interval count must be divisible by $order"))

    if order == 1
        total = 0.5 * (y[1] + y[end])
        @inbounds for i in 2:n
            total += y[i]
        end
        return h * total
    elseif order == 2
        total = y[1] + y[end]
        @inbounds for i in 2:n
            total += (iseven(i) ? 4 : 2) * y[i]
        end
        return h * total / 3
    elseif order == 3
        total = y[1] + y[end]
        @inbounds for i in 2:n
            total += (mod(i - 1, 3) == 0 ? 2 : 3) * y[i]
        end
        return 3h * total / 8
    end
    throw(ArgumentError("supported Newton-Cotes orders are 1, 2, and 3"))
end

field_energy_density(Q; c=1) = @. 0.5 * (Q[:, 1]^2 / c^2 + Q[:, 2]^2)
field_momentum_density(Q; c=1) = @. -Q[:, 1] * Q[:, 2] / c^2

field_energy(Q, h; c=1) = h * sum(field_energy_density(Q; c))
field_momentum(Q, h; c=1) = h * sum(field_momentum_density(Q; c))

@inline function particle_energy(m, v, Vz=zero(m))
    return m / sqrt(one(v) - v^2) + Vz
end

@inline function kinetic_energy(m, v)
    γ = inv(sqrt(one(v) - v^2))
    return m * (γ - one(γ))
end

end
