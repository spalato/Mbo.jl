import Base: isvalid
using IterTools

export System, states, energy, energy!, isground, dipole, dipole!, isallowed,
    allallowed, lineshape, lineshape!, amplitude, hilbert_paths

# TODO: add rotating frames
# TODO: check performance if immutable
# TODO: use NamedArray? # compare NamedArray to triangular array and custom index lookup
struct System 
    #indices::Array{String, 1}  #is there an 'ordered set'? (with unique elements?)
    grounds::Set{String} # ground states (or set)
    energies::Dict{String, Float64}  
    dipoles::Dict{Tuple{String, String}, Float64} # 
    lineshapes::Dict{Tuple{String, String}, Function}
    # TODO: some constructors
end
System() = System(Set{String}(),
                  Dict{String, Float64}(),
                  Dict{Tuple{String, String}, Float64}(),
                  Dict{Tuple{String,String}, Function}()
                  )
                  
# states
states(s::System) = keys(s.energies)
energy(s::System, tag) = s.energies[tag]
energy!(s::System, tag, v) = s.energies[tag] = v
function energy(s::System, p::HilbertPath)
    ref, a, b, c, _ = p.p
    (energy(s, a), energy(s, b), energy(s, c))
end

isground(s::System, tag) = tag in s.grounds

# dipoles
dipole(s::System, tag1, tag2) = s.dipoles[tag1, tag2]
dipole!(s::System, tag1, tag2, v) = (s.dipoles[tag1,tag2]=v; s.dipoles[tag2,tag1]=v) 
isallowed(s::System, tag1, tag2) = dipole(s, tag1, tag2) > 0
allallowed(s::System, tags) = all(isallowed.(s, tags[2:end], tags[1:end-1]))
dipole(s::System, p::HilbertPath) = dipole.(s, p.p[1:end-1], p.p[2:end])
amplitude(s::System, p::HilbertPath) = prod(dipole(s, p)...) 

# lineshapes
lineshape(s::System, tag1, tag2) = s.lineshapes[tag1, tag2]
lineshape!(s, tag1, tag2, f) = s.lineshapes[tag1, tag2] = f
function lineshape(s::System, p::HilbertPath)
    ref, a, b, c, _ = p.p
    return (lineshape(s, a, a), lineshape(s, b, b), lineshape(s, c, c),
            lineshape(s, b, a), lineshape(s, c, a), lineshape(s, c, b))
end

# hilbert paths
isvalid(s::System, p::HilbertPath) = !consecutive(p.p) && allallowed(s, p)

# TODO: test
# We'll need to profile this...
# this is clunky. Maybe Lazy.jl can help
# or make an iterator with internal graph representation and adjacency lists...
# or unroll manually with Base.Cartesian-fu
function hilbert_paths(s::System, order::Int)
    prev = Set{HilbertPath}() # keep previously seen
    return Channel() do c
        for g in s.grounds
            for i in product(fill(states(s), order)...)
                p = HilbertPath((g, i..., g))
                if isvalid(s, p) && !(p in prev)
                    push!(prev, p)
                    push!(c, p)
                end
            end
        end 
    end
end