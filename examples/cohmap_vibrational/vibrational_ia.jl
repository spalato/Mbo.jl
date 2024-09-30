# 3 level system: ground, excited and doubly excited, with Huang-Rhys coupling
# to a vibrational mode.

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
e_e = ev2angphz(cfg["e_e"])
kt = ev2angphz(cfg["kt"])
e_vib = ev2angphz(cfg["e_vib"])
e_f = ev2angphz(cfg["e_f"])
s_hr = cfg["s_hr"] 
γ = ev2angphz(cfg["gamma"])
σ = ev2angphz(cfg["sigma"])
# setup lineshape function
@inline g_bloch(t) = g_homo(t, γ) + g_inhomo(t, σ)
g=Dict()
g["e", "e"] = t -> g_bloch(t) + s_hr * g_huang_rhys(t, e_vib, kt)
g["e", "f"] = t -> 2*g_bloch(t) + 0 * s_hr * g_huang_rhys(t, e_vib, kt)
g["f", "e"] = g["e","f"]
g["f", "f"] = t -> 4*g_bloch(t) + 1 * s_hr * g_huang_rhys(t, e_vib, kt)


t1_n = cfg["t1_n"]
t2_n = cfg["t2_n"]
t3_n = cfg["t3_n"]
t1_max = cfg["t1_max"]
t2_max = cfg["t2_max"]
t3_max = cfg["t3_max"]

info("Setting up system")
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
# tic()
hpaths = collect(hilbert_paths(s, 3))
rr = zeros(ComplexF64, size(tg))
rn = zeros(ComplexF64, size(tg))
# Rephasing induced absorption is given by R1* 
# Nonrephasing IA is given by R2*
# Should be streamlined...
for p in hpaths
    info("$p")
    if p.p[3] == s.ground
        rr += R2(tg, s, p) + R3(tg, s, p)
        rn += R1(tg, s, p) + R4(tg, s, p)
    else #p.p[3] == "f"
        rr += -conj(R1(tg, s, p))
        rn += -conj(R2(tg, s, p))
    end
end

# dt = toq()
info("Calulation took $(dt) s")
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
    sa[2:end,:,:] += reverse(sr[2:end,:,:], dims=1)
else
    sa += reverse(sr, dims=1)
end

info("Saving rephasing spectrum to $(root)_sr.bin")
write("$(root)_sr.bin", sr)
info("Saving non-rephasing spectrum to $(root)_sn.bin")
write("$(root)_sn.bin", sn)
info("Saving absorptive spectrum to $(root)_sa.bin")
write("$(root)_sa.bin", sa)

end # function main

run(ARGS) # pass in command line arguments to our 'run' function
