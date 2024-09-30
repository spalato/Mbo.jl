import Base: isvalid
using IterTools

export System, states, energy, energy!, isground, dipole, dipole!, isallowed,
    allallowed, lineshape, lineshape!, amplitude, hilbert_paths

# TODO: add rotating frames
struct System 
    ground::String # ground state
    energies::Dict{String, Float64}  
    dipoles::Dict{Tuple{String, String}, Float64} 
    lineshapes::Dict{Tuple{String, String}, Any}
end
#= function System(ground_tag::String, g_energy=0.0)
    s = System()
    s.ground=ground_tag
    energy!(s, ground_tag, g_energy)
    s
end =#
function System(ground_tag::String="ground", g_energy=0.0)
    s = System(
        ground_tag,
        Dict{String, Float64}(),
        Dict{Tuple{String, String}, Float64}(),
        Dict{Tuple{String,String}, Function}()
    )
    energy!(s, ground_tag, g_energy)
    s
end
# states
states(s::System) = keys(s.energies)
energy(s::System, tag) = s.energies[tag]
energy!(s::System, tag, v) = s.energies[tag] = v
# I can most likely streamline this, (as well as lineshapes)
function energy(s::System, p::HilbertPath{3})
    # WARNING: NOT USING REF ENERGY
    ref, a, _ = p.p
    (energy(s, a),)
end
function energy(s::System, p::HilbertPath{5})
    # WARNING: NOT USING REF ENERGY
    ref, a, b, c, _ = p.p
    (energy(s, a), energy(s, b), energy(s, c))
end

isground(s::System, tag) = tag == s.ground

# dipoles
dipole(s::System, tag1, tag2) = s.dipoles[tag1, tag2]
dipole!(s::System, tag1, tag2, v) = (s.dipoles[tag1,tag2]=v; s.dipoles[tag2,tag1]=v) 
isallowed(s::System, tag1, tag2) = dipole(s, tag1, tag2) > 0
# Note (s,) "protects" it from the broadcasting done by the .
allallowed(s::System, tags) = all(isallowed.((s,), tags[2:end], tags[1:end-1]))
allallowed(s::System, p::HilbertPath) = allallowed(s, p.p)
# Note (s,) "protects" it from the broadcasting done by the .
dipole(s::System, p::HilbertPath) = dipole.((s,), p.p[1:end-1], p.p[2:end])
amplitude(s::System, p::HilbertPath) = prod(dipole(s, p)) 

# lineshapes
function lineshape(s::System, tag1, tag2)
    if !((tag1, tag2) in keys(s.lineshapes)) && any(isground(s, t) for t=(tag1, tag2))
        lineshape!(s, tag1, tag2, zero)
    end
    s.lineshapes[tag1, tag2]
end
lineshape!(s, tag1, tag2, f) = s.lineshapes[tag1, tag2] = f
function lineshape(s::System, p::HilbertPath{3})
    ref, a, _ = p.p
    return (lineshape(s, a, a),)
end

function lineshape(s::System, p::HilbertPath{5})
    ref, a, b, c, _ = p.p
    return (lineshape(s, a, a), lineshape(s, b, b), lineshape(s, c, c),
            lineshape(s, b, a), lineshape(s, c, a), lineshape(s, c, b))
end

# hilbert paths
isvalid(s::System, p::HilbertPath) = !consecutive(p.p) && allallowed(s, p)

# TODO: test
# We'll need to profile this...
# this is clunky.
# or make an iterator with internal graph representation and adjacency lists...
# or unroll manually with Base.Cartesian-fu
function hilbert_paths(s::System, order::Int)
    prev = Set{HilbertPath}() # keep previously seen
    return Channel() do c
        g=s.ground
        for i in product(fill(states(s), order)...)
            p = HilbertPath((g, i..., g))
            if isvalid(s, p) && !(p in prev)
                push!(prev, p)
                push!(c, p)
            end
        end
    end
end