# fourier transform coherences to 2d spectra
# usage: julia spec_ft.jl <config> <out>

import YAML

function run(cfgf, root, outf)
    cfg = open(YAML.load, cfgf)
    sz = tuple(map(length, cfg["t$i"] for i=1:3)...)
    inf_r = root*"_sr.bin"
    inf_n = root*"_sn.bin"
    in_r = open(inf_r, "r+")
    in_n = open(inf_n, "r+")
    s_r = Mmap.mmap(in_r, Array{ComplexF64, 3}, sz)
    s_n = Mmap.mmap(in_n, Array{ComplexF64, 3}, sz)
    out = open(outf, "w+")
    try
        spec_abs = Mmap.mmap(out, Array{ComplexF64, 3}, sz)
        spec_abs[:] = copy(s_n)
        spec_abs[2:end,:,:] += reverse(s_r[2:end,:,:], dims=1)
    finally
        close(out)
        close(in_r)
        close(in_n)
    end
end

cfgf, root, outf = ARGS[1:3]
run(cfgf, root, outf)
