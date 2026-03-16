module Equations

export LinearScalarAdvectionEquation1D
export LinearScalarWaveEquation1D

###
### Linear Scalar Advection equation
###

struct LinearScalarAdvectionEquation1D{RealT<:Real}
    advection_velocity::RealT
end

@inline function flux(u, equation::LinearScalarAdvectionEquation1D)
    a = equation.advection_velocity
    return a * u
end

###
### LINEAR SCALAR WAVE EQUATION IN FIRST ORDER REDUCTION
###

using StaticArrays

struct LinearScalarWaveEquation1D{RealT<:Real}
    velocity::RealT
end

# Calculate 1D flux in for a single point

@inline function principal_part(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    return @SMatrix [zero(c) -c * c
                     -one(c) zero(c)]
end

# @inline function right_eigenvectors(equation::LinearScalarWaveEquation1D)
#     c = equation.velocity
#     r1 = @SVector [c, one(c)] # eigenvalue = -c # left going
#     r2 = @SVector [-c, one(c)] # eigenvalue = c # right going
#     return r1, r2
# end

# @inline function left_eigenvectors(equation::LinearScalarWaveEquation1D)
#     c = equation.velocity
#     l1 = @SVector [one(c) / 2c, one(c) / 2]
#     l2 = @SVector [-one(c) / 2c, one(c) / 2]
#     return l1, l2
# end

@inline function eigenvalue_matrix(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    # return @SMatrix [c zero(c)
    #                  zero(c) -c]
    return @SMatrix [-c zero(c)
                     zero(c) c]
end

@inline function right_eigenmatrix(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    # return @SMatrix [one(c)/c one(c)
    #                  -one(c)/c one(c)]
    return @SMatrix [c -c
                     one(c) one(c)]
end

@inline function left_eigenmatrix(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    return @SMatrix [one(c)/2c one(c)/2
                     -one(c)/2c one(c)/2]
end

@inline function characteristic_variables(q, equation::LinearScalarWaveEquation1D)
    L = left_eigenmatrix(equation) # R^-1 # inverse of R
    w = L * q
    # q = [Π, Ψ]
    # w[1] right going Π - Ψ
    # w[2] left going Π + Ψ
    return w
end

@inline function Aplus(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    # A = @SMatrix [c/2 -c * c/2
    #               -one(c)/2 c/2]

    A = @SMatrix [c/2 -c * c/2
                  -one(c)/2 c/2]

    return A
end

@inline function Aminus(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    # A = @SMatrix [-c/2 -c * c/2
    #               -one(c)/2 -c/2]
    A = @SMatrix [-c/2 -c * c/2
                  -one(c)/2 -c/2]

    return A
end

@inline function Aabs(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    return @SMatrix [c zero(c)
                     zero(c) c]
end

end # end of module
