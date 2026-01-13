module NumericalFluxes

using ..Equations
using LinearAlgebra

## Scalar Advection Equation 1D
function flux_godunov(ul, ur, equation::LinearScalarAdvectionEquation1D)
    v_normal = equation.advection_velocity
    if v_normal >= 0.0
        return v_normal * ul
    else
        return v_normal * ur
    end
end

#### Scalar Wave Equation 1D
####
####
####

# Godunov flux function
# LeVeque page 83 Eq. (4.56)
# Monotone: YES
function flux_godunov(ql, qr, equation::LinearScalarWaveEquation1D)
    return Equations.Aplus(equation) * ql + Equations.Aminus(equation) * qr
end

# Roe flux function
# LeVeque page 84 Eq. (4.61)
function flux_roe(ql, qr, equation::LinearScalarWaveEquation1D)
    A = Equations.principal_part(equation)
    absA = Equations.Aabs(equation)
    return 0.5 * (A * ql + A * qr) - 0.5 * absA * (qr - ql)
end

# Lax Wendroff flux function
# LeVeque page 118 Eq. (6.48)
function flux_LxW1(ql, qr, equation::LinearScalarWaveEquation1D; kwargs...)
    cfl = get(kwargs, :cfl, 1.0)
    A = Equations.principal_part(equation)
    return 0.5 * A * (ql + qr) - 0.5 * cfl * A * A * (qr - ql)
end

# Lax Wendroff flux function
# LeVeque page 118 Eq. (6.49)
function flux_LxW2(ql, qr, equation::LinearScalarWaveEquation1D; kwargs...)
    cfl = get(kwargs, :cfl, 1.0)
    Ap = Equations.Aplus(equation)
    Am = Equations.Aminus(equation)
    absA = Equations.Aabs(equation)
    return Ap * ql + Am * qr + 0.5 * absA * (I - cfl * absA) * (qr - ql)
end

# Lax Friedrichs flux function
# Mishra page 99 Eq. (7.35)
function flux_LxF(ql, qr, equation::LinearScalarWaveEquation1D; kwargs...)
    cfl = get(kwargs, :cfl, 1.0)
    A = Equations.principal_part(equation)
    return 0.5 * A * (qr + ql) - 0.5 * cfl * (qr - ql)
end

# Rusanov flux function
# https://www.uio.no/studier/emner/matnat/math/MAT-IN9240/h17/pensumliste/numcl_notes.pdf
# Mishra page 99 Eq. (7.36)
# Monotone: YES
# is the same as godunov flux for wave equation: Exercise 7.3 Mishra
function flux_rusanov(ql, qr, equation::LinearScalarWaveEquation1D; kwargs...)
    A = Equations.principal_part(equation)
    return 0.5 * (A * (ql + qr)) - 0.5equation.velocity * (qr - ql)
end

@inline function physical_flux(q, equation::LinearScalarWaveEquation1D)
    # F(q) = A * q
    return Equations.principal_part(equation) * q
end

function flux_hll(ql, qr, equation::LinearScalarWaveEquation1D)
    # 1. Retrieve Wave Speed c
    c = equation.velocity

    # 2. Define Wave Speeds estimates (Exact for linear case)
    # S_L: Speed of the fastest wave moving to the left
    # S_R: Speed of the fastest wave moving to the right
    S_L = -c
    S_R = c

    # 3. Calculate Physical Fluxes
    fl = physical_flux(ql, equation)
    fr = physical_flux(qr, equation)

    # 4. HLL Logic
    if S_L >= 0.0
        # Flow is entirely supersonic to the right (impossible for wave eq where c > 0)
        return fl
    elseif S_R <= 0.0
        # Flow is entirely supersonic to the left (impossible for wave eq where c > 0)
        return fr
    else
        # Subsonic / Transonic case (The standard case for wave eq)
        # F_hll = (Sr*Fl - Sl*Fr + Sl*Sr*(Qr - Ql)) / (Sr - Sl)

        inv_diff = one(c) / (S_R - S_L)

        # Note: We must be careful with matrix/vector multiplication order
        # Sl, Sr are scalars. fl, fr, ql, qr are Vectors (SVector).

        return inv_diff * (S_R * fl - S_L * fr + S_L * S_R * (qr - ql))
    end
end

end #end of module
