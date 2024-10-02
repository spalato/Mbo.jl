# Migration from julia 0.6 to julia 1.10

`julia` code for versions before v1.0 were never intended to be forward
compatible. Here is a list of updates that should be applied to old calculation
scripts.

## Imports
Quite a few functions were moved out of `Base`. Therefore, the following imports
might need to be added
```julia
using Printf  # for `@sprintf`
using DelimitedFiles # for `writedlm`
using FFTW # for `ifft`, `fft`, etc.
using Mmap # for `Mmap.mmap`
```

## Name changes
The name of the following things has been changed. A simple replace-all should do.
- `Complex128` to `ComplexF64`.
- `info` to `@info`.
- `now` to `time`. Output is s, was ms.
- `linspace(start, stop, length)` to `range(start, stop, length)`. The new `range` accepts various arguments, check it out.


## Changes
The following changes are (very) slightly more complicated:  
- `tic()` and `toc()` have been removed. Instead of `tic(); dt=toq();`, one
should use `t0 = time(); dt = time()-t0;`;
- `flipdim(arr, axis)` to `reverse(arr, dims=axis)`
