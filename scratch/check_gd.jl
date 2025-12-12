using Mbo
using Plots

function plot_g!(p, t, g, color, name="")
    plot!(p[1], t, real(g), color=color, label="Re "*name)
    plot!(p[2], t, imag(g), color=color, linestyle=:dash, label="Im "*name)
end

function plot_gd()
plotlyjs(size=(800, 600))

kt = ev2angphz(0.025)
#σ = ev2angphz(0.1)
#λ = σ^2/2/kt
λ = ev2angphz(0.55/2)
σ = sqrt(2*λ*kt)
τ = 320.0
ω_B = ev2angphz(0.008)
S = λ/ω_B
@info "" λ S ω_B 2π/ω_B 1/τ
@info clamp(ceil(Int, 1/(2*π*kt*τ)), 100, 1000000)

p1 = plot()
p2 = twinx()
p = [p1, p1]
t = range(0.0, 4000.0, step=1.0)

y = g_od.(t, τ, σ, kt)
plot_g!(p, t, y, :black, "od")

#plot_g!(p, t, g_inhomo.(t, σ), :gray, "inhomo")

plot_g!(p, t, g_kubo.(t, τ, σ), :red, "Kubo")

#plot_g!(p, t, g_damped.(t, τ, σ, kt, ω_B), :blue, "Damped")

plot_g!(p, t, S*(g_huang_rhys.(t, ω_B, kt)), :orange, "HR")

g_mbo = g_from_sd(w -> S*lorentz_sd(w, ω_B, 1/τ), kt, range(-20,20, 200000))
plot_g!(p, t, (g_mbo.(t)), :cyan, "Mbo")
#ylims!(-5, 100)
end