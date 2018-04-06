
# unpack arguments
# unpack 'system' arguments
function R1(g::TimeGrid{<:Real,3}, s::System, p::HilbertPath{5})
    R1(g,
       energy(s, p),
       lineshape(s, p)
    )
end

# try another approach ... this doesn't work: it tries to broadcast the tuple too.
# rely on '@eval'?
R1(t1, t2, t3, ws::NTuple{3,Real}, ls::NTuple{6,Function}) = R1(t1, t2, t, ws..., ls...)
R1(g::TimeGrid{<:Real,3}, ws, ls) = R1.(grid(g), ws, ls)
#=
R1(g::TimeGrid{T,3}, ws::NTuple{3,Real}, ls::NTuple{6,Function}) where {T} = R1(grid(g), ws, ls)
R1(g::TimeGrid{T,3}, s::System) where {T} = sum(p->R1(g, s, p), hilbert_paths(s, 3))

# Unpack timegrid arguments
R1(t::NTuple{3,AbstractArray}, w::NTuple{3,Real}, g::NTuple{6,Function}) = 
    R1(t, w..., g...)
R1(t::NTuple{3,AbstractArray}, v...) = R1.(t[1], t[2], t[3], v...)
=#
