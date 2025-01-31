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

## Scalar Wave Equation 1D
# function flux_godunov(ΔQ, equation::LinearScalarWaveEquation1D, direction::Int)
#     if direction > 0
#         return -Equations.Aplus(equation) * ΔQ
#     elseif direction < 0
#         return Equations.Aminus(equation) * ΔQ
#     end
# end

function flux_godunov(ql, qr, equation::LinearScalarWaveEquation1D)
    return Equations.Aplus(equation) * ql + Equations.Aminus(equation) * qr
end

# function flux_godunov(ql, qr, equation::LinearScalarWaveEquation1D, direction::Int)
#     ΔQ = qr - ql
#     return flux_godunov(ΔQ, equation, direction)
# end

# flux function page 84 Eq. (4.61)
function flux_roe(ql, qr, equation::LinearScalarWaveEquation1D)
    A = Equations.principal_part(equation)
    absA = Equations.Aabs(equation)
    return 0.5 * (A * ql + A * qr) - 0.5 * absA * (qr - ql)
end

# flux function page 118 Eq. (6.48)
function flux_lax1(ql, qr, equation::LinearScalarWaveEquation1D; kwargs...)
    cfl = get(kwargs, :cfl, 1.0)
    A = Equations.principal_part(equation)
    return 0.5 * A * (ql + qr) - 0.5 * cfl * A * A * (qr - ql)
end

# flux function page 118 Eq. (6.49)
function flux_lax2(ql, qr, equation::LinearScalarWaveEquation1D; kwargs...)
    cfl = get(kwargs, :cfl, 1.0)
    Ap = Equations.Aplus(equation)
    Am = Equations.Aminus(equation)
    absA = Equations.Aabs(equation)
    return Ap * ql + Am * qr + 0.5 * absA * (I - cfl * absA) * (qr - ql)
end

end #end of module
