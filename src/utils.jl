"""
Numerically Integrate array `y` with constant spacing `h` using Simpson's rule.

For odd number of points, last interval is computed using the trapezoid rule.
"""
function simps(y, h::Number)
    n  = length(y)-1
    tail = 0
    if  n % 2 == 1
        n -= 1
        tail = (y[end]+y[end-1])/2*h
    end
    s = sum(y[1:2:n] + 4y[2:2:n]+y[3:2:n+1])
    h/3*s + tail
end
#simps(y, h::Base.TwicePrecision{T}) where {T<:AbstractFloat} = simps(y, h.hi) 
simps(y, x::Range{T}) where {T} = simps(y, step(x))

# convenience
"""Squeeze all dimensions of length 1"""
Base.squeeze(A::AbstractArray) = squeeze(A, tuple(find(size(A).==1)...))