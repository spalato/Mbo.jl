# basic example for a 3 level system g -> e -> f

using Mbo
import DSP: fftfreq
import YAML
using FFTW
using DelimitedFiles

function run(args)
cfg_f = args[1]
@info("Loading parameters from $cfg_f")
cfg = open(YAML.load, cfg_f)
# load parameters
root = cfg["rootname"]
ω_frame = ev2angphz(cfg["e_frame"])

# build lineshape functions
g = Dict()
for (tag1, tag2) in (["e", "e"], ["e", "f"], ["f","f"])
    γ = ev2angphz(cfg["gamma_$(tag1)$(tag2)"])
    σ = ev2angphz(cfg["sigma_$(tag1)$(tag2)"])
    g[tag1, tag2] = t -> g_homo(t, γ) + g_inhomo(t, σ)
end
g["f", "e"] = g["e", "f"]

t1_n = cfg["t1_n"]
t2_n = cfg["t2_n"]
t3_n = cfg["t3_n"]
t1_max = cfg["t1_max"]
t2_max = cfg["t2_max"]
t3_max = cfg["t3_max"]
e_e = ev2angphz(cfg["e_e"])
e_f = ev2angphz(cfg["e_f"])

@info("Setting up system")
s = System("g")
# rotating frame is a bit artisanal. It shows here.
energy!(s, "e",  e_e - ω_frame)
energy!(s, "f", e_f - 2*ω_frame)
dipole!(s, "g", "e", 1.0)
dipole!(s, "g", "f", 0.0)
dipole!(s, "e", "f", 1.0)
for tag1 in ["e", "f"], tag2 in ["e", "f"]
    lineshape!(s, tag1, tag2, g[tag1, tag2])
end

tg = TimeGrid(
    range(0, stop=t1_max, length=t1_n),
    range(0, stop=t2_max, length=t2_n),
    range(0, stop=t3_max, length=t3_n),
)

@info("Computing linear response")
r_lin = linear(tg, s)
@info("Saving to $(root)_rlin.txt")
writedlm("$(root)_rlin.txt", [tg.times[1] real(r_lin) imag(r_lin)])

r_lin[1] *= 0.5
s_lin = fftshift(ifft(r_lin))
f_lin = fftshift(fftfreq(size(tg)[1], 1/(tg.times[1][2]-tg.times[1][1])))
@info("Saving linear spectrum to $(root)_slin.txt")
writedlm("$(root)_slin.txt", [f_lin real(s_lin) imag(s_lin)])

@info("Computing third order response")
t0 = time()
hpaths = collect(hilbert_paths(s, 3))
rr = zeros(ComplexF64, size(tg))
rn = zeros(ComplexF64, size(tg))
# Rephasing induced absorption is given by R1* 
# Nonrephasing IA is given by R2*
# Should be streamlined...
for p in filter(hp->hp.p[3] == "g", hpaths)
    rr += R2(tg, s, p) + R3(tg, s, p)
    rn += R1(tg, s, p) + R4(tg, s, p)
end
for p in filter(hp->hp.p[3] == "f", hpaths)
    rr += -conj(R1(tg, s, p))
    rn += -conj(R2(tg, s, p))
end

dt = time()-t0
@info("Calulation took $(dt) s")
@info("Saving to $(root)_rr.bin, $(root)_rn.bin")
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
    sa[2:end,:,:] += reverse(sr[2:end,:,:], dims=1)
else
    sa += reverse(sr, dims=1)
end

@info("Saving rephasing spectrum to $(root)_sr.bin")
write("$(root)_sr.bin", sr)
@info("Saving non-rephasing spectrum to $(root)_sn.bin")
write("$(root)_sn.bin", sn)
@info("Saving absorptive spectrum to $(root)_sa.bin")
write("$(root)_sa.bin", sa)

end # function main

run(ARGS) # pass in command line arguments to our 'run' function