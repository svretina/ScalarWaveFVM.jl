module NumericalFluxes

using ..Equations

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
function flux_godunov(ΔQ, equation::LinearScalarWaveEquation1D, direction::Int)
    if direction > 0
        return -Equations.Aplus(equation) * ΔQ
    elseif direction < 0
        return Equations.Aminus(equation) * ΔQ
    end
end

function flux_godunov(ql, qr, equation::LinearScalarWaveEquation1D, direction::Int)
    ΔQ = qr - ql
    return flux_godunov(ΔQ, equation, direction)
end

end #end of module
