module Reconstructions

export downwind_slope, upwind_slope, forward_slope, centered_slope, compute_limited_slope

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

@inline function compute_limited_slope(qi, ql, qr, h, slope_limiter)
    σD = downwind_slope(qi, ql, qr, h)
    σU = upwind_slope(qi, ql, qr, h)
    # σC = central_slope(qi, ql, qr, h)
    θ = σU / σD
    σ = slope_limiter(θ) * σD
    return σ
end

# Slopes
@inline function constant_slope(qi, ql, qr, h)
    return zero(qi)
end

@inline function centered_slope(qi, ql, qr, h)
    return (qr - ql) / (2h)
end

@inline function upwind_slope(qi, ql, qr, h)
    return (qi - ql) / h
end

@inline function downwind_slope(qi, ql, qr, h)
    return (qr - qi) / h
end

end
