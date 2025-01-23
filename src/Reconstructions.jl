module Reconstructions

@inline function piecewise_constant(u)
    return u
end

@inline function piecewise_linear(xi, ui, ul, ur, params)
    p(x) = ui + params.slope_func(ui, ul, ur, params.h) * (x - xi)
    return p
end

#Slopes
@inline function constant_slope(u, ul, ur, h)
    return zero(u)
end

@inline function central_slope(u, ul, ur, h)
    return (ur - ul) / (2h)
end

@inline function backward_slope(u, ul, ur, h)
    return (u - ul) / h
end

@inline function forward_slope(u, ul, ur, h)
    return (ur - u) / h
end

end
