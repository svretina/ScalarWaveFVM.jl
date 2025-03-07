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

# same as rusanov
# function flux_HLLE(ql, qr, equation::LinearScalarWaveEquation1D; kwargs...)

#     A = Equations.principal_part(equation)
#     return (c * A * ql + c * A * qr - c * c(qr - ql)) / (2c)
# end

end #end of module
