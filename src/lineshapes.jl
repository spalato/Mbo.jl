export g_homo, g_inhomo, g_huang_rhys, g_kubo, LineshapeLUT, g_od, g_od_slow, g_damped

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

function g_od(t, τ, σ, kt) # overdamped. From eqs 26-29 in https://doi.org/10.1021/ja00103a020
    t == 0.0 && return 0.0+1im*0.0 # type stability!
    λ = σ^2/2/kt
    gi = λ*τ*(1-exp(-t/τ))
    nmin = clamp(ceil(Int, 1/(2*π*kt*τ)), 100, 1000000)
    gr = λ*τ*cot(1/2/τ/kt)*(exp(-t/τ)+t/τ-1) + 4*λ*kt/τ * converge_series(k->matsu_freq(k, kt, t, τ), nmin=nmin)
    gr+1im*gi
end

#= Doesn't work...
function damped_freq(n, kt, t, τ, ω_b)
    nun = 2*π*kt*n
    nunt = nun*t
    nun2 = nun*nun
    numerator = (exp(-nunt)+nunt-1)
    denominator = nun*(nun2-((1+nun2/(ω_b^2))/τ)^2)
    (numerator/denominator)::Real
end

function g_damped(t, τ, σ, kt, ω_B)
    t == 0.0 && return 0.0+1im*0.0 # type stability!
    λ = σ^2/2/kt
    gi = λ*τ*(1-exp(-t/τ))
    nmin = clamp(ceil(Int, 1/(2*π*kt*τ)), 100, 1000000)
    gr = λ*τ*cot(1/2/τ/kt)*(exp(-t/τ)+t/τ-1) + 4*λ*kt/τ * converge_series(k->damped_freq(k, kt, t, τ, ω_B), nmin=nmin)
    gr+1im*gi
end =#

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

(g::LineshapeLUT)(t) = g.lut[Int(round(t/g.dx; digits=4))+1]
# we should probably define the addition of multiple lineshape types. Whatever.
