module Equations

export LinearScalarAdvectionEquation1D
export LinearScalarWaveEquation1D

###
### Linear Scalar Advection equation
###

struct LinearScalarAdvectionEquation1D{RealT<:Real}
    advection_velocity::RealT
end

###
### LINEAR SCALAR WAVE EQUATION IN FIRST ORDER REDUCTION
###

using StaticArrays

struct LinearScalarWaveEquation1D{RealT<:Real}
    velocity::RealT
end

# Calculate 1D flux in for a single point
@inline function flux(u, equation::LinearScalarAdvectionEquation1D)
    a = equation.advection_velocity
    return a * u
end

@inline function principal_part(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    A = zeros(typeof(c), 2, 2)
    A[1, 2] = -c * c
    A[2, 1] = -one(c)
    return A
end

@inline function right_eigenvectors(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    r1 = @SVector [c, one(c)] # eigenvalue = -c # left going
    r2 = @SVector [-c, one(c)] # eigenvalue = c # right going
    return r1, r2
end

@inline function left_eigenvectors(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    l1 = @SVector [one(c) / 2c, one(c) / 2]
    l2 = @SVector [-one(c) / 2c, one(c) / 2]
    return l1, l2
end

@inline function characteristic_variables(q, equation::LinearScalarWaveEquation1D)
    R = @SMatrix [c -c
                  one(c) one(c)]
    L = @SMatrix [one(c)/(2c) 0.5
                  -one(c)/(2c) 0.5] # R^-1 # inverse of R
    w = L * q
    # q = [Π, Ψ]
    # w[1] left going
    # w[2] right going
    return w
end

@inline function Aplus(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
    A = @SMatrix [c/2 -c * c/2
                  -one(c)/2 c/2]
    return A
end

@inline function Aminus(equation::LinearScalarWaveEquation1D)
    c = equation.velocity
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
