# plot linear data
# usage: julia plot_linear.jl <config>
import YAML
using Plots
using Mbo.TimeGrid
using Rsvg
plotlyjs()

function main(args)
cfgf = args[1]
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
);
rootname = cfg["rootname"]

# linear plotsa
dat_t = readdlm("$(rootname)_rlin.txt")
t = dat_t[:,1]
rlin = dat_t[:,2]+1im*dat_t[:,3]
p = plot(t, [real(rlin), imag(rlin)])
png(p, "$(rootname)_rlin.png")

dat_f = readdlm("$(rootname)_slin.txt")
f = dat_f[:,1]
slin = dat_f[:,2]+1im*dat_f[:,3]
p = plot(f, [real(slin), imag(slin)])
png(p, "$(rootname)_slin.png")
end # function main

main(ARGS)