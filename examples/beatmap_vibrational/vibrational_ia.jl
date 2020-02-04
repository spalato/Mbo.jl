# 3 level system: ground, excited and doubly excited, with Huang-Rhys coupling
# to a vibrational mode.

using Mbo
import DSP: fftfreq
import YAML

function run(args)
cfg_f = args[1]
info("Loading parameters from $cfg_f")
cfg = open(YAML.load, cfg_f)
# load parameters
ω_frame = ev2angphz(cfg["e_frame"])
e0 = ev2angphz(cfg["e_eg"])
kt = ev2angphz(cfg["kt"])
e_vib = ev2angphz(cfg["e_vib"])
ωf = 2*e0-ev2angphz(cfg["delta"])
s_hr = cfg["s_hr"] 
γ = ev2angphz(cfg["gamma"])
σ = ev2angphz(cfg["sigma"])
# setup lineshape function
@inline g_bloch(t) = g_homo(t, γ) + g_inhomo(t, σ)
@inline g_xx(t) = g_bloch(t) + s_hr * g_huang_rhys(t, e_vib, kt)

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
    lineshape!(s, tag1, tag2, g_xx) # !!! not sure if this will work, I have to think about lineshape funcitons again
end
# !!! keep on working here