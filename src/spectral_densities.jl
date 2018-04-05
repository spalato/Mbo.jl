export lorentz_sd, mbo_sd, g_from_sd

"""
g_integrand(c::Function, w::AbstractFloat, kt::AbstractFloat)

Integrand for the calculation of `g`` on a constant grid.
"""
function g_integrand(cw, t, w, kt)
wt = w*t
ingr = (
    (1+coth(w/2/kt))
    *cw/w^2
    *(cis(-wt)+1im*wt -1)
)
ingr
end

"""
g_from_sd(c::Function, kt)

Generate a lineshape function `g(t)` from odd part of spectral density
function ``c(ω)``.
"""
function g_from_sd(c::Function, kt::AbstractFloat, wgrid=linspace(-20,20, 2000))
    cw = c.(wgrid)
    function g(t::AbstractFloat)
        -simps(g_integrand.(cw, t, wgrid, kt), wgrid)/2/pi
    end
    g
end

function lorentz_sd(ω, ω_0, γ)
    n = ω*γ
    l = ω*ω-ω_0*ω_0-γ*γ
    d = l*l+n*n*4
    n/d*ω_0^3*4
end

function mbo_sd(ω, ω_0, γ)
    n = ω*γ
    l = ω*ω-ω_0*ω_0
    d = l*l+n*n
    n/d*ω_0^3*2*sqrt(2)
end