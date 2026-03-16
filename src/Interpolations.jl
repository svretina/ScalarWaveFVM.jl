module Interpolations

export LinearInterpolation1D
export CubicInterpolation1D
export CubicHermiteInterpolation1D

@inline function _uniform_step(x::AbstractVector{T}) where {T<:Real}
    return abs(x[2] - x[1])
end

@inline function _segment_index_uniform(xnew::T, x::AbstractVector{T}, h::T) where {T<:Real}
    if xnew == x[begin]
        return 1
    end
    i0 = floor(Int64, (xnew - x[begin]) / h) + 1
    return max(1, min(i0, length(x) - 1))
end

# assumes uniform grid h=const
function LinearInterpolation1D(xnew::T, x::AbstractVector{T},
                               y::AbstractVector{T}) where {T<:Real}
    h = _uniform_step(x)
    i0 = _segment_index_uniform(xnew, x, h)
    return y[i0] * ((x[i0 + 1] - xnew) / h) + y[i0 + 1] * ((xnew - x[i0]) / h)
end

"""
    CubicInterpolation1D(xnew, x, y)

Fourth-order local stencil (4 points) cubic Lagrange interpolation on a uniform grid.
Falls back to linear interpolation if fewer than 4 points are available.
"""
function CubicInterpolation1D(xnew::T, x::AbstractVector{T},
                              y::AbstractVector{T}) where {T<:Real}
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have the same length"))
    n >= 2 || throw(ArgumentError("at least 2 points are required"))
    n < 4 && return LinearInterpolation1D(xnew, x, y)

    h = _uniform_step(x)
    i0 = _segment_index_uniform(xnew, x, h)

    # Center the 4-point stencil around the active segment.
    i = clamp(i0, 2, n - 2)
    τ = (xnew - x[i]) / h

    @inbounds y_m1 = y[i - 1]
    @inbounds y_0 = y[i]
    @inbounds y_p1 = y[i + 1]
    @inbounds y_p2 = y[i + 2]

    # Uniform-grid cubic Lagrange basis with nodes τ = -1, 0, 1, 2
    l_m1 = -τ * (τ - 1) * (τ - 2) / 6
    l_0 = (τ + 1) * (τ - 1) * (τ - 2) / 2
    l_p1 = -(τ + 1) * τ * (τ - 2) / 2
    l_p2 = (τ + 1) * τ * (τ - 1) / 6

    return l_m1 * y_m1 + l_0 * y_0 + l_p1 * y_p1 + l_p2 * y_p2
end

@inline function _endpoint_or_center_slope(y::AbstractVector{T}, i::Int, h::T) where {T<:Real}
    n = length(y)
    if i <= 1
        @inbounds return (y[2] - y[1]) / h
    elseif i >= n
        @inbounds return (y[n] - y[n - 1]) / h
    else
        @inbounds return (y[i + 1] - y[i - 1]) / (2h)
    end
end

@inline function _cubic_hermite_segment(xnew::T, xi::T, h::T,
                                        y0::T, y1::T, m0::T, m1::T) where {T<:Real}
    t = (xnew - xi) / h
    t2 = t * t
    t3 = t2 * t
    h00 = 2t3 - 3t2 + 1
    h10 = t3 - 2t2 + t
    h01 = -2t3 + 3t2
    h11 = t3 - t2
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
end

"""
    CubicHermiteInterpolation1D(xnew, x, y)

Piecewise cubic Hermite interpolation on a uniform grid with finite-difference slopes:
centered in the interior and one-sided at boundaries.
Falls back to linear interpolation if fewer than 3 points are available.
"""
function CubicHermiteInterpolation1D(xnew::T, x::AbstractVector{T},
                                     y::AbstractVector{T}) where {T<:Real}
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have the same length"))
    n >= 2 || throw(ArgumentError("at least 2 points are required"))
    n < 3 && return LinearInterpolation1D(xnew, x, y)

    h = _uniform_step(x)
    i0 = _segment_index_uniform(xnew, x, h)

    @inbounds y0 = y[i0]
    @inbounds y1 = y[i0 + 1]
    m0 = _endpoint_or_center_slope(y, i0, h)
    m1 = _endpoint_or_center_slope(y, i0 + 1, h)

    return _cubic_hermite_segment(xnew, x[i0], h, y0, y1, m0, m1)
end

"""
    CubicHermiteInterpolation1D(xnew, x, y, dy)

Piecewise cubic Hermite interpolation on a uniform grid using user-provided derivatives `dy`.
"""
function CubicHermiteInterpolation1D(xnew::T, x::AbstractVector{T},
                                     y::AbstractVector{T},
                                     dy::AbstractVector{T}) where {T<:Real}
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have the same length"))
    n == length(dy) || throw(ArgumentError("dy must have the same length as x"))
    n >= 2 || throw(ArgumentError("at least 2 points are required"))

    h = _uniform_step(x)
    i0 = _segment_index_uniform(xnew, x, h)

    @inbounds y0 = y[i0]
    @inbounds y1 = y[i0 + 1]
    @inbounds m0 = dy[i0]
    @inbounds m1 = dy[i0 + 1]

    return _cubic_hermite_segment(xnew, x[i0], h, y0, y1, m0, m1)
end

end # end of module
