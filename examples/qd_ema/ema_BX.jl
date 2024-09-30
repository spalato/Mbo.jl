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

# build lineshape functions
#g_bloch(t, γ, σ) = g_homo(t, γ) + g_inhomo(t, σ)
γ_e = ev2angphz(cfg["gamma_e"])
γ_h1 = ev2angphz(cfg["gamma_h1"])
γ_h2 = ev2angphz(cfg["gamma_h2"])
σ_x1 = ev2angphz(cfg["sigma_x1"])
r_12 = cfg["sigma_r12"]
σ_x2 = r_12*σ_x1

g_e(t) = g_homo(t, γ_e)
g_h1(t) = g_homo(t, γ_h1)
g_h2(t) = g_homo(t, γ_h2)
g_size(t) = 0.5*t^2
g = Dict()
g["x1", "x1"] = t-> g_e(t)+g_h1(t)+σ_x1^2*g_size(t) #g_bloch(t, γ_e+γ_h1, σ_x1)
g["x2", "x2"] = t-> g_e(t)+g_h2(t) + σ_x2^2*g_size(t) #g_bloch(t, γ_e+γ_h2, r_12*σ_x1)
g["x1", "x2"] = t-> g_e(t) + σ_x1*σ_x2*g_size(t) #g_bloch(t, γ_e, sqrt(r_12)*σ_x1)
g["bx12", "bx12"] = t-> g_h1(t)+g_h2(t)+ 4*g_e(t) + (σ_x1+σ_x2)^2*g_size(t) #g_bloch(t, γ_h1+γ_h2+4γ_e, σ_x1*(1+r_12))
g["x1", "bx12"] = t-> g_h1(t)+2g_e(t) + (σ_x1^2+σ_x1*σ_x2)*g_size(t) #g_bloch(t, γ_h1+2γ_e, σ_x1*sqrt(1+r_12))
g["x2", "bx12"] = t-> g_h2(t)+2g_e(t) + (σ_x2^2+σ_x1*σ_x2)*g_size(t) #g_bloch(t, γ_h2+2γ_e, σ_x1*sqrt(r_12*(1+r_12)))
g11 = g["x1", "x1"]
g22 = g["x2", "x2"]
g["x1", "bx11"] = t-> 2*g11(t)
g["bx11", "bx11"] = t -> 4*g11(t)
g["x2", "bx22"] = t -> 2*g22(t)
g["bx22", "bx22"] = t -> 4*g22(t)

# do: g["x2", "x1"] = g["x1", "x2"]
for (t1, t2) in collect(keys(g))
    if t1 != t2
        g[t2, t1] = g[t1, t2]
    end
end

t1_n = cfg["t1_n"]
t2_n = cfg["t2_n"]
t3_n = cfg["t3_n"]
t1_max = cfg["t1_max"]
t2_max = cfg["t2_max"]
t3_max = cfg["t3_max"]

info("Setting up system")
s = System("g")
for tag1 in ["x1", "x2"], tag2 in ["x1", "x2"]
    energy!(s, tag1, ev2angphz(cfg["e_$(tag1)"]) - ω_frame)
    dipole!(s, "g", tag1, cfg["mu_$(tag1)"])
    tag2 != tag1 && dipole!(s, tag1, tag2, 0)
end
# generate 3 exciton states, with the same binding energies and tdms.
energy!(s, "bx12", energy(s, "x1")+energy(s, "x2")-ev2angphz(cfg["binding"]))
energy!(s, "bx11", 2*energy(s, "x1")-ev2angphz(cfg["binding"]))
energy!(s, "bx22", 2*energy(s, "x2")-ev2angphz(cfg["binding"]))
# disallow direct transitions from ground state to biexcitons
dipole!(s, "g", "bx12", 0)
dipole!(s, "g", "bx11", 0)
dipole!(s, "g", "bx22", 0)
# set tdm for biexcitons
dipole!(s, "x1", "bx12", cfg["mu_bx"])
dipole!(s, "x2", "bx12", cfg["mu_bx"])
dipole!(s, "x2", "bx22", cfg["mu_bx"])
dipole!(s, "x1", "bx22", 0)
dipole!(s, "x1", "bx11", cfg["mu_bx"])
dipole!(s, "x2", "bx11", 0)
# disallow transitions between biexcitons
for tag1=["bx11", "bx12", "bx22"], tag2=["bx11", "bx12", "bx22"]
    dipole!(s, tag1, tag2, 0)
end

tg = TimeGrid(
    linspace(0, t1_max, t1_n),
    linspace(0, t2_max, t2_n),
    linspace(0, t3_max, t3_n),
)
steps = map(first ∘ diff, tg.times)
info("Grid steps: $(steps)")
# We cache the value of the lineshape functions on a look-up table.
lut_step = 1.0
@assert all([rem(s, lut_step)==0 for s in steps])
lut_grid = 0.0:lut_step:sum(map(maximum, (t1_max, t2_max, t3_max)))
for (tag1, tag2) in keys(g)
    lineshape!(s, tag1, tag2, LineshapeLUT(g[tag1, tag2].(lut_grid), step(lut_grid)))
end

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
# Rephasing induced absorption is given by R1* 
# Nonrephasing IA is given by R2*
# Should be streamlined...
rr = zeros(Complex128, size(tg))
rn = zeros(Complex128, size(tg))
for hp in hilbert_paths(s, 3)
    # setting mu_bx=0 in the config file will skip the paths automatically
    info("Path: $hp") 
    if hp.p[3] in ["bx12", "bx11", "bx22"] # is ESA
        rr += -conj(R1(tg, s, hp))
        rn += -conj(R2(tg, s, hp))
    else
        #nothing
        rr += R2(tg, s, hp) + R3(tg, s, hp)
        rn += R1(tg, s, hp) + R4(tg, s, hp)
    end
end
dt = toq()
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
    sa[2:end,:,:] += flipdim(sr[2:end,:,:], 1)
else
    sa += flipdim(sr, 1)
end

info("Saving rephasing spectrum to $(root)_sr.bin")
write("$(root)_sr.bin", sr)
info("Saving non-rephasing spectrum to $(root)_sn.bin")
write("$(root)_sn.bin", sn)
info("Saving absorptive spectrum to $(root)_sa.bin")
write("$(root)_sa.bin", sa)

end # function main

run(ARGS) # pass in command line arguments to our 'run' function