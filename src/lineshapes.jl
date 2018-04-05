export g_homo, g_inhomo, g_huang_rhys, g_kubo

#=                  LINESHAPE FUNCTIONS                                 =#
g_homo(t, γ) = t*γ
g_inhomo(t, σ) = 0.5*σ^2*t^2
function g_huang_rhys(t, omega_0, kt)
    w0t = omega_0*t
    coth(omega_0/kt/2.0)*(1-cos(w0t))+1im*(sin(w0t)-w0t)
end
g_kubo(t, τ, σ) = (σ*τ)^2*(exp(-t/τ)+t/τ-1) # TODO: TRY

