module Reconstructions

function piecewise_constant(u)
    return u
end

function piecewise_linear(xj, uj, ul, ur, params)
    p(x, u) = u + params.slope_func(u, ul, ur, params.h) * (x - xj)
    return p
end

#Slopes
function constant_slope(u, ul, ur, h)
    return 0.0
end

function central_slope(u, ul, ur, h)
    return (ur - ul) / (2h)
end

function backward_slope(u, ul, ur, h)
    return (u - ul) / h
end

function forward_slope(u, ul, ur, h)
    return (ur - u) / h
end

# function central_slope(u, i, h)
#     return (u[i + 1] - u[i - 1]) / (2h)
# end

# function backward_slope(u, i, h)
#     return (u[i] - u[i - 1]) / h
# end

# function forward_slope(u, i, h)
#     return (u[i + 1] - u[i]) / h
# end

#Slope Limiters
function minmod(u, i, h)
    upwind = forward_slope(u, i, h)
    downwind = backward_slope(u, i, h)
    return minmod(upwind, downwind)
end

function minmod(σ1, σ2)
    if sign(σ1) == sign(σ2)
        return min(σ1, σ2)
    else
        return zero(σ1)
    end
end

function maxmod(σ1, σ2)
    if sign(σ1) == sign(σ2)
        return max(σ1, σ2)
    else
        return zero(σ1)
    end
end

function superbee(u, i, h)
    σ1 = backward_slope(u, i, h)
    σ2 = forward_slope(u, i, h)
    σL = minmod(2σ1, σ2)
    σR = minmod(σ1, 2σ2)
    return maxmod(σL, σR)
end

end
