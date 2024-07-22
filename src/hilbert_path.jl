import Base: ==, <, isless, show, print, isvalid

export HilbertPath, order, isvalid

# Hilbert paths
"""
    consecutive(array)

True if any element is identical to the next. Naive implementation
"""
consecutive(a) = any(a[2:end] .== a[1:end-1])

"""Hilbert Path of order N-2"""
struct HilbertPath{N}
    p::NTuple{N,String}
end
HilbertPath(args::Vararg{String,N}) where {N} = HilbertPath{N}(args)
order(p::HilbertPath{N}) where {N} = N-2
show(io::IO, p::HilbertPath{N}) where {N} = print(io, "HilbertPath{$N}$(p.p)")
print(io::IO, p::HilbertPath) = print(io, ("<", join(p.p, ","), ">")...)

function ==(p1::HilbertPath, p2::HilbertPath)
    order(p1) == order(p2) && p1.p == p2.p
end
# make sortable
isless(p1::HilbertPath{N}, p2::HilbertPath{N}) where {N} = p1.p < p2.p

