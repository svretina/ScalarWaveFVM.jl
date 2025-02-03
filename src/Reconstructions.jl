module Reconstructions

export downwind_slope, upwind_slope, forward_slope, centered_slope

@inline function piecewise_constant(u)
    return u
end

@inline function piecewise_linear(xi, ui, ul, ur, params)
    p(x) = ui + params.slope_func(ui, ul, ur, params.h) * (x - xi)
    return p
end

@inline function reconstruct_at_faces(qi, ql, qr, h)
    σ = compute_slope(qi, ql, qr, h)
    slope_term = 0.5σ * h
    qL = qi - slope_term
    qR = qi + slope_term
    return qL, qR
end

@inline function compute_slope(qi, ql, qr, h)
    σF = forward_slope(qi, ql, qr, h)
    σB = backward_slope(qi, ql, qr, h)
    σC = central_slope(qi, ql, qr, h)
    σ = slope_limiter(σF, σB)
    return σ
end

# Slopes
@inline function constant_slope(u, ul, ur, h)
    return zero(u)
end

@inline function centered_slope(u, ul, ur, h)
    return (ur - ul) / (2h)
end

@inline function upwind_slope(u, ul, ur, h)
    return (u - ul) / h
end

@inline function forward_slope(u, ul, ur, h)
    return (ur - u) / h
end

@inline function downwind_slope(u, ul, ur, h)
    return forward_slope(u, ul, ur, h)
end

end
