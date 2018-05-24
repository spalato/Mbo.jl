export TimeGrid, grid#, size

# check StaticArrays
struct TimeGrid{T,N}
    times::NTuple{N,Array{T,1}}
end
TimeGrid(args::NTuple{N,AbstractArray{T,1}}) where {T,N} = TimeGrid{T,N}(args)
TimeGrid(args...) = TimeGrid(args)
function grid(tg::TimeGrid{T,N}) where {T,N}
    Tuple(reshape(t_,
             Tuple(setindex!(ones(Int32, N), length(t_), i)))
          for (i, t_) in enumerate(tg.times))
end
Base.size(tg::TimeGrid) = map(length, tg.times)