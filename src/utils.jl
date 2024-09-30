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

"""
Sum a series until convergence.
"""
function converge_series(f; nmin=10, nmax=100_000, rtol=1e-9)
    has_converged(a, s, k) = k <= nmax ? 
                                (nmin < k < nmax) && abs(a/s) < rtol :
                                 throw(ErrorException("Failed to converge in $nmax, $k")) # TODO: in 1.? we can return nothing
    k=1
    a_k = f(k)
    s_k = a_k
    while !has_converged(a_k, s_k, k)
        k += 1
        a_k = f(k)
        s_k += a_k
    end
    s_k
end

# convenience
"""Squeeze all dimensions of length 1"""
Base.squeeze(A::AbstractArray) = squeeze(A, tuple(find(size(A).==1)...))