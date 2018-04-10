
# unpack arguments
# unpack 'system' arguments
R1(g::TimeGrid{<:Real,3}, s::System, p::HilbertPath{5}) = R1(
    g, energy(s, p), lineshape(s, p))


R1(g::TimeGrid{T,3}, s::System) where {T} = sum(p->R1(g, s, p), hilbert_paths(s, 3))
# allow unpacking both ways.
R1(t1, t2, t3, ws::NTuple{3,Real}, ls::NTuple{6,Function}) = R1(t1, t2, t3, ws..., ls...)
R1(g::TimeGrid{<:Real,3}, ws, ls) = R1(grid(g)..., ws..., ls...)
R1(t1::Array{<:Real,3}, t2::Array{<:Real,3}, t3::Array{<:Real,3}, a...) = R1.(t1, t2, t3, a...)
