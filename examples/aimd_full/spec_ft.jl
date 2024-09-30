# fourier transform coherences to 2d spectra
# usage: julia spec_ft.jl <config> <in> <out>

import YAML

function run(cfgf, infn, outf)
    cfg = open(YAML.load, cfgf)
    sz = tuple(map(length, cfg["t$i"] for i=1:3)...)
    inf = open(infn, "r+")
    sig = Mmap.mmap(inf, Array{Complex128, 3}, sz)
    sig[1,:,:] *= 0.5
    sig[:,:,1] *= 0.5
    out = open(outf, "w+")
    spec2d = Mmap.mmap(out, Array{Complex128, 3}, sz)
    fill!(spec2d, 0)
    spec2d[:] = fftshift(ifft(sig, (1,3)), (1,3))
    close(out)
    sig[1,:,:] *= 2
    sig[:,:,1] *= 2
    close(inf)
end

cfgf, infn, outf = ARGS[1:3]
run(cfgf, infn, outf)
