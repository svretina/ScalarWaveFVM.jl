module NumericalFluxes

# numerical fluxes
function flux_central(ul, ur, flux)
    return (flux(ul) + flux(ur)) / 2
end

function flux_upwind(ul, ur)
    return equation.advection_velocity * (ur - ul)
end

function dissipation_LxF()
end

# is this with cfl = dx/dt? or the advection speed?
function flux_lax_friedrichs(flux, u, h, dt)
    return flux_central(ul, ur, flux) - (cfl / 2) * (u[i + 1] - u[i])
end

#monotone fluxes
function flux_godunov(ul, ur, equation)
    v_normal = equation.advection_velocity
    if v_normal >= 0.0
        return v_normal * ul
    else
        return v_normal * ur
    end
end

function flux_rusanov end

function flux_engquist_osher(ul, ur, equation)
    return flux_central(ul, ur) - 0.5 * abs(equation.advection_velocity) * (ul - ur)
end

end
