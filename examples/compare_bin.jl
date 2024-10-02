import YAML
using Mbo
using Printf

function run(args)
cfgf = args[1]
@debug "Opening: " cfgf
dir = dirname(cfgf)
cfg = open(YAML.load, cfgf)
rootname = cfg["rootname"]
try
    t1_n = cfg["t1_n"]
    t2_n = cfg["t2_n"]
    t3_n = cfg["t3_n"]
    t1_max = cfg["t1_max"]
    t2_max = cfg["t2_max"]
    t3_max = cfg["t3_max"]
    global tg = TimeGrid(
        range(0, t1_max, t1_n),
        range(0, t2_max, t2_n),
        range(0, t3_max, t3_n),
    )
catch err
    if isa(err, KeyError)
        t1 = convert.(Float64, cfg["t1"])
        t2 = convert.(Float64, cfg["t2"])
        t3 = convert.(Float64, cfg["t3"])
        global tg = TimeGrid(t1, t2, t3)
    else rethrow()
    end
end

for sig in ["rr", "rn"]
    data = Array{ComplexF64}(undef, size(tg))
    read!(joinpath(dir, "$(rootname)_$(sig).bin"), data)

    ref = Array{ComplexF64}(undef, size(tg))
    read!(joinpath(dir, "$(rootname)_$(sig).bin.ref06"), ref)

    if all(isapprox.(data, ref))
        @info "$(cfgf) $(sig) ok       abs="*@sprintf("%.06e", maximum(abs.(data-ref)))
    else
        @info "$(cfgf) $(sig) DIFF!!!  abs="*@sprintf("%.06e", maximum(abs.(data-ref)))
    end
end

end
run(ARGS)