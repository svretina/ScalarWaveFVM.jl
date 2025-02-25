module Interpolations

# assumes uniform grid h=const
function LinearInterpolation1D(xnew::T, x::AbstractVector{T},
                               y::AbstractVector{T}) where {T<:Real}
    h = abs(x[2]-x[1])
    if xnew != x[begin]
        i0 = ceil(Int64, (xnew - x[begin]) / h)
    else
        i0 = 1
    end
    return y[i0] * ((x[i0 + 1] - xnew) / h) + y[i0 + 1] * ((xnew - x[i0]) / h)
end




end # end of module
