import YAML
using Mbo
using ProgressMeter

tag(i) = @sprintf "X%02d" i
function parse_corr(fn)
    readdlm(fn)
end

function run(args)
t0 = now()
cfgf = args[1]
info("Loading parameters from $(cfgf)")
cfg = open(YAML.load, cfgf)
states_fn = cfg["states_fn"]
cf_fn = cfg["xcorr_fn"]

t1 = convert.(Float64, cfg["t1"])
t2 = convert.(Float64, cfg["t2"])
t3 = convert.(Float64, cfg["t3"])
sigma = cfg["sigma"]

info("Loading energies from $states_fn")
states = readdlm(states_fn)
energies = states[:,1]
tdm = states[:,2]
info("States count: $(size(states)[1])")
info("Loading correlation matrix from $cf_fn")
corrm = parse_corr(cf_fn)
info("Correlation count: $(length(corrm))")

# build system
info("Building system...")
t0_sys = now()
s = System("G")
frame = ev2angphz(cfg["e_frame"])
lut_grid = 0.0:1.0:sum(map(maximum, (t1, t2, t3)))
for i=1:size(corrm)[1], j=1:size(corrm)[2]
    ti = tag(i)
    tj = tag(j)
    energy!(s, ti, ev2angphz(energies[i])-frame)
    dipole!(s, "G", ti, tdm[i])
    dipole!(s, ti, tj, 0)
    lineshape!(s, "G", ti, zero)
    lineshape!(s, ti, "G", zero)
    #if i == j
        lut = LineshapeLUT(t->(g_inhomo(t, ev2angphz(corrm[i,j]))
                               +g_homo(t, ev2angphz(0.005))
                               +g_inhomo(t, ev2angphz(sigma))), lut_grid)
    #else
    # we assume perfect correlation along the inhomogeneous coordinate.
    #   lut = LineshapeLUT(t->(g_inhomo(t, ev2angphz(corrm[i,j]))
    #                        +g_inhomo(t, ev2angphz(sigma))), lut_grid)
    #end
    lineshape!(s, ti, tj, lut)
end
info("    Took $(now()-t0_sys)")
info("Number of order 1 Hilbert Paths: $(length(collect(hilbert_paths(s, 1))))")
info("Number of order 3 Hilbert Paths: $(length(collect(hilbert_paths(s, 3))))")

grd_lin = TimeGrid(t1)

info("Computing linear responses separately.")
t0_calc = now()
out_root = "$(splitext(basename(states_fn))[1])_$(splitext(basename(cf_fn))[1])"

# compute them separately
totlin = zeros(Complex128, size(grd_lin))
for p in hilbert_paths(s, 1)
    rlin = linear(grd_lin, s, p)
    totlin += rlin
    out_lin = out_root*"_lin_"*join(p.p[2:end-1], "-")*".txt"
    info("Saving linear to: $out_lin")
    writedlm(out_lin, [grid(grd_lin)[1] real(rlin) imag(rlin)])
end
writedlm(out_root*"_lin_tot.txt", [grid(grd_lin)[1] real(totlin) imag(totlin)])
info("    Took $(now()-t0_calc)")
info("Preparing calculation of third order response")
grd_trd = TimeGrid(t1, t2, t3)
info("Grid size: $(size(grd_trd))")

rr_fn = out_root*"_rr.bin"
rn_fn = out_root*"_rn.bin"
info("Saving rephasing to: $rr_fn")
info("Saving nonrephasing to: $rr_fn")
rr_out = open(rr_fn, "w+")
rn_out = open(rn_fn, "w+")
rr = Mmap.mmap(rr_out, Array{Complex128, 3}, size(grd_trd))
rn = Mmap.mmap(rn_out, Array{Complex128, 3}, size(grd_trd))
pm = Progress(4*length(collect(hilbert_paths(s, 3))))
info("Computing third order response")
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