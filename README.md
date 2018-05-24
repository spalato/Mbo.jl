<!--This puppy is not ready for the world yet. -->

# Mbo.jl - Multimode Brownian Oscillator [WIP]

Semiclassical modeling of non-linear spectra in the time domain using the
Multimode Brownian Oscillator model, in julia. Designed for the modeling of 2D visible
spectroscopy. Hopefully fast, hopefully easy to use.

## Description

The aim of this code is to faciliate calculations of 2D spectra, putting simple
model systems within reach of the humble experimental spectroscopist. The focus
is on molecular-like system, using the cumulant expansion and classical bath
modes. The centerpiece is a toolkit built around the spectroscopic system: the
user builds a system by specifing state energies, transition dipole moments and
lineshape functions. The linear and non-linear responses can then be
conveniently calculated from the system itself.

For example:
```julia
# build a system. It's still a bit tedious, but that's where the thinking happens.
s = System()
push!(s.grounds, "g") # todo: remove this one ASAP
# Add states "g" and "e" with energies 0 and 1.5.
energy!(s, "g", 0) # you can read this as `set energy for "g" to 0`
energy!(s, "e", 1.5)
# Set the transition dipole moment for g->e
dipole!(s, "g", "e", 1.0)
# Define a lineshape functions (some are supplied).
ls(t) = 0.01*t + 0.5*(0.015*t)^2
lineshape!(s, "g", "e", ls)
# setup a calculation grid
tg = TimeGrid(linspace(0, 100, 100), linspace(0, 100, 100), linspace(0, 100, 100))
# compute
rr = R2(tg, s) + R3(tg, s) # rephasing
rn = R1(tg, s) + R4(tg, s) # non-rephasing
# go ahead and process to your wishes.
```
The example above uses hand-coded constants, but you are free to use any input
or output format you desire. Same goes for lineshape functions: they can be
defined arbitrarily using the programming language (some common forms are
predefined). This allows for flexibility: the system's parameters can equally be
hard coded, read from a file, or computed from more complicated algorithms (for
example, diagonalizing a site basis to an exciton basis). Julia can call
code written in `C`, `Fortran` and `python`.

<!--
The purpose of this code is *not* to build a modeling suite, but to make this
tedious calculation easier to script. It assumes the user is capable of hacking their way through data using `Matlab` or `python`. Modeling suites often fall
into the trap of language design: they slowly grow to support multiple features
and make up file formats and input mini-languages as they go. This is avoided
here by making use of the `julia` programming language to retain the simpler
syntax of `Matlab` and `python` but the computing power of compiled languages.
-->

# Installation
This packages uses the [`julia`](https://julialang.org/) programming language.
You'll need it, but it is easy to install.

`Julia` is an  open-source language  aims to provide clear syntax, similar to 
`Matlab` or `python` while providing performance typical of compiled languages
such as `Fortran` or `C`. Hopefully, the language is easy to pick up for anyone
familiar with any of these. Julia is available for Windows, Mac, Linux. 
Get it from: https://julialang.org/downloads/ (v0.7 or later untested. Please tell me if it breaks!)

Once julia is installed, open `julia` and type:

```julia
julia> Pkg.clone("https://github.com/spalato/Mbo.jl")
```

# Documentation
The main element of this package is the `System` object, which contains all the
physical parameters of the system under study. For the cumulant expansion, this
means state energies, transition dipole moments and lineshape functions. The 
states are represented internally as strings.

When computing a response function, the code builds all the valid paths in
Hilbert space (or Feynmann diagrams) directly from the system. The desired
signals are then summed over all such valid paths (ie: transition dipoles all
greater than 0). The equations for the signals R1 to R4 were taken from
the book by [Cho], eq 5.26-5.30.

## Basics

The attributes are set and read using functions. Functions that modify the
`System` object end with `!`. (This is purely a convention in julia, the `!` has no
special meaning.)

```julia
using Mbo
s = System()
# add a ground state, should get rid of this ASAP
push!(s.grounds, "g")
energy!(s, "g", 0)

# set energies
energy!(s, "x1", 1.2)
energy!(s, "x2", 1.4)
# read energies:
@assert energy(s, "x2")-energy(s, "x1") ≈ 0.2

# set transition dipole momments.
dipole!(s, "g", "x1", 1.0) # this sets both g->x1 and x1->g
dipole!(s, "g", "x2", 1.5)
# IMPORTANT: if nonzero, the x1->x2 transition will be part of the calculation
dipole!(s, "x1", "x2", 0)
# read transition dipole
@assert dipole(s, "g", "x1") == dipole(s, "x1", "g")

# set lineshape functions
f(t) = 0.1*t
lineshape!(s, "g", "g", zero) # todo: remove ASAP.
lineshape!(s, "x1", "g", zero) # todo: remove ASAP
lineshape!(s, "x2", "g", zero) # todo: remove ASAP
lineshape!(s, "x1", "x1", f) # accept any callable.
lineshape!(s, "x2", "x2", t->g_homo(t, 0.1)) # anonymous function in Julia: x->...
lineshape!(s, "x1", "x2", zero) # uncorrelated.
# read lineshape function
@assert lineshape(s, "x1", "x1") === f
```

The spectroscopic response is calculated in the time domain. You can use the
`TimeGrid` object to define the domain of calculation.

```julia
t1=0:20.0 # 0 to 20 included in steps of 1.0
t2=0:10.0:100.0 # 0 to 100 included in steps of 10.0
t3=0:20.0
tg = TimeGrid(t1, t2, t3)
```

The calculation can be performed automatically for all points of the grid and 
for all paths. Alternatively, you can generate the Hilbert Paths and compute
for specific paths manually.
```julia
# compute everyting
rtot = R1(tg, s) + R2(tg, s) + R3(tg, s) + R4(tg, s)
# check hilbert paths
hpaths = collect(hilbert_paths(s, 3)) # third order hilbert paths
@assert length(hpaths) == 4
# compute for a specific path, pick the first.
p = hpaths[1]
r1p = R1(tg, s, p)
...
```

## Features
Currently supported:
- Linear response in the time domain.
- Third order response.
- Computation for all pathways, or a subset.
- Set energies, transition dipoles and lineshape functions.
- Arbitrary callables for lineshape functions 
- Built-in lineshapes: Homogeneous, Inhomogeneous, Huang-Rhys, Kubo Ansatz.
- Lineshape functions from correlation functions.
- Caching lineshape functions for improved performance.
- IO using julia: delimited files, YAML, binary, memory-mapped.

Works, but a bit clumsy:
- Rotating frames.
- Lineshapes from spectral densities.

# Citing
If you found this useful, please consider citing the following paper:
> [WIP!]

# Examples

## Scripts

# Help!
Your best bet is to open an issue using the 'Issues' tab at the top. 

# TODO
- [ ] Handle ground states easily.
- [ ] Compilations to LUT
- [ ] Add rotating frames properly.
- [ ] Add automatic rephasing vs non-rephasing (`rephasing(grid,system)`).
- [ ] Convenience for zeroing transition dipole moments (`forbidden`?)
- [ ] Convenience for correlating (or uncorrelating) lineshape functions?
<!-- 
- [ ] More documentation...
- [ ] Add filtering of Hilbert Paths (eg: IA only).-->
<!--
## May happen...
- Explicit bath modes
-->
[Cho]: Coh, M. *Two-dimensional optical spectroscopy.* (CRC Press, 2009).