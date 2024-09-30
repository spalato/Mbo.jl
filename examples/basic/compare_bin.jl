import YAML
using Mbo

function run(args)
cfgf = args[1]
dir = dirname(cfgf)
cfg = open(YAML.load, cfgf)
rootname = cfg["rootname"]
t1_n = cfg["t1_n"]
t2_n = cfg["t2_n"]
t3_n = cfg["t3_n"]
t1_max = cfg["t1_max"]
t2_max = cfg["t2_max"]
t3_max = cfg["t3_max"]
tg = TimeGrid(
    range(0, t1_max, t1_n),
    range(0, t2_max, t2_n),
    range(0, t3_max, t3_n),
)

for sig in ["rr", "rn"]
    data = Array{ComplexF64}(undef, size(tg))
    read!(joinpath(dir, "$(rootname)_$(sig).bin"), data)

    ref = Array{ComplexF64}(undef, size(tg))
    read!(joinpath(dir, "$(rootname)_$(sig).bin.ref06"), ref)

    if all(isapprox.(data, ref))
        @info "$(cfgf) ok       abs:$(maximum(abs.(data-ref)))"  
    else
        @info "$(cfgf) DIFF!!!  abs:$(maximum(abs.(data-ref)))"
    end
end

end
run(ARGS)