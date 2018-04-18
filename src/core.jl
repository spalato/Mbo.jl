
# Allow flexibility in calling the reponse functions.
# Use the @eval macro to operate on R1 to R4.
for R=(:R1, :R2, :R3, :R4)
    @eval begin
    #=unpack 'system' arguments
    It's important to unpack the lineshape functions early: dispatch will then
    ensure each combination of lineshape function gets optimized.
    =#
    $(R)(g::TimeGrid{<:Real,3}, s::System, p::HilbertPath{5}) = $(R)(g, energy(s, p), lineshape(s, p))
    $(R)(g::TimeGrid{T,3}, s::System) where {T} = sum(p->$(R)(g, s, p), hilbert_paths(s, 3))
    # allow unpacking both ways.
    $(R)(t1, t2, t3, ws::NTuple{3,Real}, ls::NTuple{6,Function}) = $(R)(t1, t2, t3, ws..., ls...)
    $(R)(g::TimeGrid{<:Real,3}, ws, ls) = $(R)(grid(g)..., ws..., ls...)
    $(R)(t1::Array{<:Real,3}, t2::Array{<:Real,3}, t3::Array{<:Real,3}, a...) = $(R).(t1, t2, t3, a...)
    #$(R)(g::TimeGrid{<:Real,3}, wa, wb, wc, gaa, gbb, gcc, gba, gca, gcb) = $(R).(grid(g)..., wa, wb, wc, gaa, gbb, gcc, gba, gca, gcb)
    end
end
# TODO: test all of these!

#function Rephasing(g::TimeGrid, ws, ls)
#    if ws[2] > ws[1] # IA
#        return conj(R1(g, ws, ls))
#    else
#        return ()
#end