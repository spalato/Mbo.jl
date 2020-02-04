# Use time-correlation functions to obtain lineshape.
export GriddedCF

# externally computed on a constant grid. Allow specific points only.
immutable GriddedCF{Tg<:Real, Tcf<:Number}
    lut::Array{Tcf,1}
    Im::Tcf
    dt::Tg
    tm::Tg # max time
end
function GriddedCF(grid, dt) # performs the double integral on initialization
    grid = copy(grid)
    grid[1] *= 0.5
    I1 = cumsum_kbn(grid).*dt
    # TODO: check next line is ok: cumsum_kbb(I1) vs cumsum_kbb(I1).*dt
    GriddedCF(cumsum_kbn(I1), I1[end], dt, dt*length(grid))
end
# make callable
function (cf::GriddedCF)(t)
    i = Int(t/cf.dt)+1
    i <= length(cf.lut) ? cf.lut[i] : extrapolate(cf, t) # try to inbounds this
end
extrapolate(cf::GriddedCF, t) = (t-cf.tm)*cf.Im + cf.lut[end]
# test! compare to manual zero padding.