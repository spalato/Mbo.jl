export g_homo, g_inhomo, g_huang_rhys, g_kubo, LineshapeLUT, g_overdamped

#=                  LINESHAPE FUNCTIONS                                 =#
g_homo(t, γ) = t*γ
g_inhomo(t, σ) = 0.5*σ^2*t^2
function g_huang_rhys(t, omega_0, kt)
    w0t = omega_0*t
    coth(omega_0/kt/2.0)*(1-cos(w0t))+1im*(sin(w0t)-w0t)
end
g_kubo(t, τ, σ) = (σ*τ)^2*(exp(-t/τ)+t/τ-1)
function g_overdamped(t, τ, σ, kt)
    λ = σ^2/2/kt
    g_kubo(t, τ, σ)+1im*λ*τ*(1-exp(-t/τ))
end

#abstract type Lineshape end
# Tried with LazyList. HORRIBLY SLOW (they don't carry type information)
# this is VERY efficient. We should find a way to make this automatic
struct LineshapeLUT{Tx<:Real, Tls<:Number}# <: Lineshape
    lut::Array{Tls,1}
    dx::Tx
    
end
function LineshapeLUT(lut::Array{Tls,1}, dx::Tx) where {Tls, Tx}
    LineshapeLUT{Tx, Tls}(lut, dx)
end
function LineshapeLUT(f, x::Range{Tx}) where {Tx}#, Tls}
    gx = f.(x)
    LineshapeLUT{Tx,eltype(gx)}(gx, step(x)) 
end

(g::LineshapeLUT)(t) = g.lut[Int(t/g.dx)+1]
# we should probably define the addition of multiple lineshape types. Whatever.
