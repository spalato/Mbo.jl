# plot linear data
# usage: julia plot_linear.jl <config>
import YAML
using Plots
using Mbo.TimeGrid
using Rsvg
plotlyjs()

function main(args)
cfgf = args[1]
dir = dirname(cfgf)
cfg = open(YAML.load, cfgf)
rootname = cfg["rootname"]

# linear plotsa
dat_t = readdlm(joinpath(dir, "$(rootname)_rlin.txt"))
t = dat_t[:,1]
rlin = dat_t[:,2]+1im*dat_t[:,3]
p = plot(t, [real(rlin), imag(rlin)])
png(p, "$(rootname)_rlin.png")

dat_f = readdlm(joinpath(dir, "$(rootname)_slin.txt"))
f = dat_f[:,1]
slin = dat_f[:,2]+1im*dat_f[:,3]
p = plot(f, [real(slin), imag(slin)])
png(p, "$(rootname)_slin.png")
end # function main

main(ARGS)