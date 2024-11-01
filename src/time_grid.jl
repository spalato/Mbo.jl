export TimeGrid, grid, t_to_f, freqs
using FFTW
import Statistics: mean

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

t_to_f(t::AbstractRange) = fftshift(fftfreq(length(t), 1/step(t)))
# use mean step. Hopefully doesn't matter.
t_to_f(t::AbstractArray{<:Real, 1}) = fftshift(fftfreq(length(t), 1/mean(diff(t))))
"""
    freqs(tg::TimeGrid)

Get FFT frequencies from a `TimeGrid`.

For an `N`-dimensional `TimeGrid`, return `N` arrays.
"""
freqs(tg::TimeGrid) = map(t_to_f, tg.times)