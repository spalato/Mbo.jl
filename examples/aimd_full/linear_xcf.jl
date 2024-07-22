import YAML
using Mbo
using ProgressMeter
import DSP: fftfreq

tag(i) = @sprintf "X%02d" i
function parse_cf(fn)
    dat = readdlm(fn)
    cf = Dict{NTuple{2,Int},Array{<:Real,1}}()
    for i=1:size(dat)[1]
        cf[Int(dat[i,1]), Int(dat[i,2])] = dat[i,3:end]
    end
    cf
end

function run(args)
t0 = now()
cfgf = args[1]
@info("Loading parameters from $(cfgf)")
cfg = open(YAML.load, cfgf)
states_fn = cfg["states_fn"]
cf_fn = cfg["cf_fn"]

t1 = convert.(Float64, cfg["t1"])
t2 = convert.(Float64, cfg["t2"])
t3 = convert.(Float64, cfg["t3"])
σ = cfg["sigma"]

@info("Loading energies from $states_fn")
states = readdlm(states_fn)
energies = states[:,1]
tdm = states[:,2]
@info("States count: $(size(states)[1])")
@info("Loading CF from $cf_fn")
cfs = parse_cf(cf_fn)
@info("CF count: $(length(cfs))")

# build system
@info("Building system...")
t0_sys = now()
s = System("G")
#energy!(s, "G", 0)
#push!(s.grounds, "G")
#lineshape!(s, "G", "G", zero)
frame = ev2angphz(cfg["e_frame"])
dt = 3.0
lut_grid = 0.0:dt:sum(map(maximum, (t1, t2, t3)))
for (i, j) in keys(cfs)
    ti = tag(i)
    tj = tag(j)
    energy!(s, ti, ev2angphz(energies[i])-frame)
    dipole!(s, "G", ti, tdm[i])
    dipole!(s, ti, tj, 0)
    lineshape!(s, "G", ti, zero)
    lineshape!(s, ti, "G", zero)
    # we assume perfect correlation along the inhomogeneous coordinate.
    lut = LineshapeLUT(t->(GriddedCF(cfs[i,j], dt)(t)+g_inhomo(t, ev2angphz(σ))), lut_grid)
    lineshape!(s, ti, tj, lut)
end
@info("    Took $(now()-t0_sys)")
@info("Number of order 1 Hilbert Paths: $(length(collect(hilbert_paths(s, 1))))")
@info("Number of order 3 Hilbert Paths: $(length(collect(hilbert_paths(s, 3))))")

grd_lin = TimeGrid(t1)

t0_calc = now()
out_root = cfg["rootname"]

# compute them separately
totlin = zeros(Complex128, size(grd_lin))
for p in hilbert_paths(s, 1)
    rlin = linear(grd_lin, s, p)
    #writedlm(out_root*"_lin_$(p.p[2]).txt", [grid(grd_lin)[1] real(rlin) imag(rlin)])
    totlin += rlin
end
@info("Saving to $(out_root)_rlin.txt")
writedlm("$(out_root)_rlin.txt", [grid(grd_lin)[1] real(totlin) imag(totlin)])

totlin[1] *= 0.5
s_lin = fftshift(ifft(totlin))
f_lin = fftshift(fftfreq(size(grd_lin)[1], 1/(grd_lin.times[1][2]-grd_lin.times[1][1])))
@info("Saving linear spectrum to $(out_root)_slin.txt")
writedlm("$(out_root)_slin.txt", [f_lin real(s_lin) imag(s_lin)])

@info("    Took $(now()-t0_calc)")
@info("Total runtime: $(now()-t0)")
end # run

run(ARGS)