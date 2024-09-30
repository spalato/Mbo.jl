# basic example for a 2 level system
# introduces important concepts, both for julia and Mbo
# 

using Mbo     # using <module> loads functions directly, ie: `ev2angphz` or `R1`
import DSP: fftfreq # use a specific function
import YAML   # import <module> requires adding the module name, ie: YAML.load
using IterTools
using FFTW
using DelimitedFiles

function run(args)
# Use a main 'run' function. This is not mandatory, but is required for performance.
cfg_f = args[1] # index by 1, 1:2 translates to [1, 2], the first two elements
@info("Loading parameters from $cfg_f")
cfg = open(YAML.load, cfg_f)
# load parameters
# load parameters
root = cfg["rootname"]
ω_frame = ev2angphz(cfg["e_frame"]) # look carefuly, this is omega, not w. Julia uses unicode
ω_a = ev2angphz(cfg["e_a"])
σ = ev2angphz(cfg["sigma"])
γ = ev2angphz(cfg["gamma"])

t1_n = cfg["t1_n"]
t2_n = cfg["t2_n"]
t3_n = cfg["t3_n"]
t1_max = cfg["t1_max"]
t2_max = cfg["t2_max"]
t3_max = cfg["t3_max"]

# setup system
@info("Setting up system")
s = System("g")
energy!(s, "a", ω_a - ω_frame) # use a rotating frame
dipole!(s, "g", "a", 1.0)
lineshape!(s, "a", "a", t->g_homo(t, γ)+g_inhomo(t, σ)) # use an anonymous function for the lineshape

tg = TimeGrid(
    range(0, stop=t1_max, length=t1_n), # TODO previously linspace
    range(0, stop=t2_max, length=t2_n),
    range(0, stop=t3_max, length=t3_n),
)

@info("Computing linear response")
# compute linear
r_lin = linear(tg, s)
@info("Saving to $(root)_rlin.txt")
writedlm("$(root)_rlin.txt", [tg.times[1] real(r_lin) imag(r_lin)])

r_lin[1] *= 0.5 # first interval is a half-interval
s_lin = fftshift(ifft(r_lin))
f_lin = fftshift(fftfreq(size(tg)[1], 1/(tg.times[1][2]-tg.times[1][1])))
@info("Saving linear spectrum to $(root)_slin.txt")
writedlm("$(root)_slin.txt", [f_lin real(s_lin) imag(s_lin)])

# compute rephasing and non-rephasing separately: the spectra will overlap if
# a fully rotating frame is used (ω_a == ω_frame)
@info("Computing third order response")
rr = R2(tg, s) + R3(tg, s)
rn = R1(tg, s) + R4(tg, s)

# save multidimensional data in binary format.
# Julia uses Fortran-contiguous arrays.
@info("Saving to $(root)_rr.bin, $(root)_rn.bin")
write("$(root)_rr.bin", rr)
write("$(root)_rn.bin", rn)
# first interval is a half-interval
rr[1,:,:] *= 0.5
rr[:,:,1] *= 0.5
rn[1,:,:] *= 0.5
rn[:,:,1] *= 0.5
sr = fftshift(ifft(rr, (1,3)), (1,3))
sn = fftshift(ifft(rn, (1,3)), (1,3))

# # add non-rephasing and a copy of rephasing, flipped along first axis.
# # you can do FFT and IFFT, but be careful with normalizations.
sa = copy(sn)
if iseven(size(sa, 1))
    sa[2:end,:,:] += reverse(sr[2:end,:,:], dims=1)
else
    sa += reverse(sr, dims=1)
end # if-else blocks terminated by 'end'. This is true of all code blocks.

@info("Saving rephasing spectrum to $(root)_sr.bin")
write("$(root)_sr.bin", sr)
@info("Saving non-rephasing spectrum to $(root)_sn.bin")
write("$(root)_sn.bin", sn)
@info("Saving absorptive spectrum to $(root)_sa.bin")
write("$(root)_sa.bin", sa)

end # terminating the `function run` code block

run(ARGS) # pass in command line arguments to our 'run' function