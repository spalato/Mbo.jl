# basic example for a 2 level system coupled to a vibrational mode.

using Mbo
import DSP: fftfreq
import YAML

function run(args)
cfg_f = args[1]
info("Loading parameters from $cfg_f")
cfg = open(YAML.load, cfg_f)
# load parameters
root = cfg["rootname"]
ω_frame = ev2angphz(cfg["e_frame"])
ω_a = ev2angphz(cfg["e_a"])
σ = ev2angphz(cfg["sigma"])
tau = cfg["tau"]
γ = ev2angphz(cfg["gamma"])

t1_n = cfg["t1_n"]
t2_n = cfg["t2_n"]
t3_n = cfg["t3_n"]
t1_max = cfg["t1_max"]
t2_max = cfg["t2_max"]
t3_max = cfg["t3_max"]

info("Setting up system")
s = System("g")
energy!(s, "a", ω_a - ω_frame) 
dipole!(s, "g", "a", 1.0)
lineshape!(s, "a", "a", t->g_homo(t, γ)+g_kubo(t, tau, σ))

# Nothing has changed past that point!

tg = TimeGrid(
    linspace(0, t1_max, t1_n),
    linspace(0, t2_max, t2_n),
    linspace(0, t3_max, t3_n),
)

info("Computing linear response")
r_lin = linear(tg, s)
info("Saving to $(root)_rlin.txt")
writedlm("$(root)_rlin.txt", [tg.times[1] real(r_lin) imag(r_lin)])

r_lin[1] *= 0.5
s_lin = fftshift(ifft(r_lin))
f_lin = fftshift(fftfreq(size(tg)[1], 1/(tg.times[1][2]-tg.times[1][1])))
info("Saving linear spectrum to $(root)_slin.txt")
writedlm("$(root)_slin.txt", [f_lin real(s_lin) imag(s_lin)])

info("Computing third order response")
tic()
rr = R2(tg, s) + R3(tg, s)
rn = R1(tg, s) + R4(tg, s)

info("Saving to $(root)_rr.bin, $(root)_rn.bin")
write("$(root)_rr.bin", rr)
write("$(root)_rn.bin", rn)
rr[1,:,:] *= 0.5
rr[:,:,1] *= 0.5
rn[1,:,:] *= 0.5
rn[:,:,1] *= 0.5
sr = fftshift(ifft(rr, (1,3)), (1,3))
sn = fftshift(ifft(rn, (1,3)), (1,3))

sa = copy(sn)
if iseven(size(sa, 1))
    sa[2:end,:,:] += flipdim(sr[2:end,:,:], 1)
else
    sa += flipdim(sr, 1)
end
dt = toq()
info("Calulation took $(dt) s")
info("Saving rephasing spectrum to $(root)_sr.bin")
write("$(root)_sr.bin", sr)
info("Saving non-rephasing spectrum to $(root)_sn.bin")
write("$(root)_sn.bin", sn)
info("Saving absorptive spectrum to $(root)_sa.bin")
write("$(root)_sa.bin", sa)

end # function main

run(ARGS) # pass in command line arguments to our 'run' function