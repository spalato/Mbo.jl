# plot linear data
# usage: julia plot_2d.jl <config>
import YAML
using Plots
using Mbo.TimeGrid
import DSP: fftfreq
using Rsvg
plotlyjs()

function sym_heatmap(x, y, z, fname)
    lim = maximum(abs.(z))
    png(plot(x, y, z, st=:heatmap, c=:balance, clim=(-lim, lim)), fname)
end

function main(args)
cfgf = args[1]
dir = dirname(cfgf)
cfg = open(YAML.load, cfgf)
t1_n = cfg["t1_n"]
t2_n = cfg["t2_n"]
t3_n = cfg["t3_n"]
t1_max = cfg["t1_max"]
t2_max = cfg["t2_max"]
t3_max = cfg["t3_max"]
tg = TimeGrid(
    linspace(0, t1_max, t1_n),
    linspace(0, t2_max, t2_n),
    linspace(0, t3_max, t3_n),
)
rootname = cfg["rootname"]
f1 = fftshift(fftfreq(length(tg.times[1]), 1/diff(tg.times[1])[1]))
f3 = fftshift(fftfreq(length(tg.times[3]), 1/diff(tg.times[3])[1]))

rr = read(joinpath(dir, "$(rootname)_rr.bin"), ComplexF64, size(tg))
sym_heatmap(tg.times[1], tg.times[3], real(rr[:,1,:]), "$(rootname)_rr.png")
rn = read(joinpath(dir, "$(rootname)_rn.bin"), ComplexF64, size(tg))
sym_heatmap(tg.times[1], tg.times[3], real(rn[:,1,:]), "$(rootname)_rn.png")
sr = read(joinpath(dir, "$(rootname)_sr.bin"), ComplexF64, size(tg))
sym_heatmap(f1, f3, real(sr[:,1,:]), "$(rootname)_sr.png")
sn = read(joinpath(dir, "$(rootname)_sn.bin"), ComplexF64, size(tg))
sym_heatmap(f1, f3, real(sn[:,1,:]), "$(rootname)_sn.png")
sa = read(joinpath(dir, "$(rootname)_sa.bin"), ComplexF64, size(tg))
sym_heatmap(f1, f3, real(sa[:,1,:]), "$(rootname)_sa.png")

end

main(ARGS)