# core API

# unpack arguments
# unpack 'system' arguments
function R1(g::TimeGrid{<:Real,3}, s::System, p::HilbertPath{5})
    R1(g,
       energy(s, p),
       lineshape(s, p)
    )
end

# try another approach
R1(t1, t2, t3, ws::NTuple{3,Real}, ls::NTuple{6,Function}) = R1(t1, t2, t3, ws..., ls...)
R1(g::TimeGrid{<:Real,3}, ws, ls) = R1.(grid(g)..., ws, ls)

