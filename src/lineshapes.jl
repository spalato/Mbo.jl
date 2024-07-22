export g_homo, g_inhomo, g_huang_rhys, g_kubo, LineshapeLUT, g_od, g_od_slow

#=                  LINESHAPE FUNCTIONS                                 =#
g_homo(t, γ) = t*γ
g_inhomo(t, σ) = 0.5*σ^2*t^2
function g_huang_rhys(t, omega_0, kt)
    w0t = omega_0*t
    coth(omega_0/kt/2.0)*(1-cos(w0t))+1im*(sin(w0t)-w0t)
end
g_kubo(t, τ, σ) = (σ*τ)^2*(exp(-t/τ)+t/τ-1)
function g_od_slow(t, τ, σ, kt) # valid if kT >> hbar/tau
    λ = σ^2/2/kt
    g_kubo(t, τ, σ)+1im*λ*τ*(1-exp(-t/τ))
end

function matsu_freq(n, kt, t, τ)
    nun = 2*π*kt*n
    (exp(-nun*t)+nun*t-1)/(nun*(nun^2-1/τ^2))::Real
end

function g_od(t, τ, σ, kt) # overdamped
    t == 0.0 && return 0.0
    λ = σ^2/2/kt
    gi = λ*τ*(1-exp(-t/τ))
    nmin = ceil(Int, 1/(2*π*kt*τ))
    gr = λ*τ*cot(1/2/τ/kt)*(exp(-t/τ)+t/τ-1) + 4*λ*kt/τ * converge_series(k->matsu_freq(k, kt, t, τ), nmin=nmin)
    gr+1im*gi
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
function LineshapeLUT(f, x::AbstractRange{Tx}) where {Tx}#, Tls}
    gx = f.(x)
    LineshapeLUT{Tx,eltype(gx)}(gx, step(x)) 
end

(g::LineshapeLUT)(t) = g.lut[Int(t/g.dx)+1]
# we should probably define the addition of multiple lineshape types. Whatever.
