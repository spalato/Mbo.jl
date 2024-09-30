import YAML
using Mbo
using ProgressMeter

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
root = cfg["rootname"]
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
frame = ev2angphz(cfg["e_frame"])
lut_grid = 0.0:1.0:sum(map(maximum, (t1, t2, t3)))
for (i, j) in keys(cfs)
    ti = tag(i)
    tj = tag(j)
    energy!(s, ti, ev2angphz(energies[i])-frame)
    dipole!(s, "G", ti, tdm[i])
    dipole!(s, ti, tj, 0)
    lineshape!(s, "G", ti, zero)
    lineshape!(s, ti, "G", zero)
    # we assume perfect correlation along the inhomogeneous coordinate.
    lut = LineshapeLUT(t->(GriddedCF(cfs[i,j], 1.0)(t)+g_inhomo(t, ev2angphz(σ))), lut_grid)
    lineshape!(s, ti, tj, lut)
end
@info("    Took $(now()-t0_sys)")
@info("Number of order 1 Hilbert Paths: $(length(collect(hilbert_paths(s, 1))))")
@info("Number of order 3 Hilbert Paths: $(length(collect(hilbert_paths(s, 3))))")

@info("Preparing calculation of third order response")
grd_trd = TimeGrid(t1, t2, t3)
@info("Grid size: $(size(grd_trd))")

rr_fn = root*"_rr.bin"
rn_fn = root*"_rn.bin"
@info("Saving rephasing to: $rr_fn")
@info("Saving nonrephasing to: $rr_fn")
rr_out = open(rr_fn, "w+")
rn_out = open(rn_fn, "w+")
rr = Mmap.mmap(rr_out, Array{ComplexF64, 3}, size(grd_trd))
rn = Mmap.mmap(rn_out, Array{ComplexF64, 3}, size(grd_trd))
pm = Progress(4*length(collect(hilbert_paths(s, 3))))
@info("Computing third order response")
t0_third = now()
rr .+= R2(grd_trd, s, pm)
rr .+= R3(grd_trd, s, pm)
rn .+= R1(grd_trd, s, pm)
rn .+= R4(grd_trd, s, pm)
info("    Took $(now()-t0_third)")
close(rr_out)
close(rn_out)
info("Total runtime: $(now()-t0)")
end # run

run(ARGS)