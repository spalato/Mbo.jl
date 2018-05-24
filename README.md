<!--This puppy is not ready for the world yet. -->

Mbo.jl - Multimode Brownian Oscillator [WIP]
======================================

Semiclassical modeling of non-linear spectra in the time domain using the
Multimode Brownian Oscillator model, in julia. Designed for the modeling of 2D
visible spectroscopy. Hopefully fast, hopefully easy to use.

Description
-----------

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
using Mbo
s = System("g") # use "g" as a ground state
# Set energy for "e" to 1.5 eV, converted to angular PHz
energy!(s, "e", ev2angphz(1.5)) 
# Set the transition dipole moment for g<->e
dipole!(s, "g", "e", 1.0)
# Define a lineshape functions (some common cases are supplied).
ls(t) = 0.002*t + 0.5*(0.002*t)^2
lineshape!(s, "e", "e", ls)
# setup a calculation grid
tg = TimeGrid(0.0:100, 0.0:100, 0.0:100)
# compute. tedious bookkeeping done by the code.
lin = linear(tg, s) # linear response
rr = R2(tg, s) + R3(tg, s) # rephasing
rn = R1(tg, s) + R4(tg, s) # non-rephasing
# This gives you a 3D array of complex numbers.
# Go ahead and process to your wishes.
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

Installation
============
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

Documentation
============
The main element of this package is the `System` object, which contains all the
physical parameters of the system under study. For the cumulant expansion, this
means state energies, transition dipole moments and lineshape functions.

Physics
-------
The calculation of the third order response is tedious: a given third order 
response function requires 3 time arguments, 3 frequencies, 6 lineshape
functions and 4 transition dipoles. This has to be repeated for every possible
path the system can take in Hilbert space. For `N` states (+ a single ground 
state), there are between `N^2` and approx `N^3` such paths (depending on 
bandwidth constraints and transition dipoles). This tedium is taken care of
using the `System` and `HilbertPath` objects. Handling of the time arguments
is simplified using the `TimeGrid` object.

When computing a response function, the code builds all the valid paths in
Hilbert space directly from the system. The transition energies, dipoles and
lineshape functions are taken from the system directly. The desired
signals are then summed over all such valid paths (ie: transition dipoles all
greater than 0).

The equations for third order responses R1 to R4 were taken from
the book by [Cho](#book_cho), eq 5.26-5.31. These equations take the fluctuating
ground state as a reference. The lineshape functions are for the energies of
the state (not the energy gaps). As a consequence, fluctuations of the ground
state are always 0.

Basics
------
The attributes of the system are set and read using functions. Functions that
modify the `System` object end with `!`. (This is purely a convention in Julia,
the `!` has no special meaning.) The energies and transition dipoles are rather
straightforward.
```julia
using Mbo
s = System("g")
# set energies
energy!(s, "x1", 1.2)
energy!(s, "x2", 1.4)
# read energies:
@assert energy(s, "x2")-energy(s, "x1") == 1.4-1.2
@assert length(states(s)) == 3
# set transition dipole momments.
dipole!(s, "g", "x1", 1.0) # this sets both g->x1 and x1->g
dipole!(s, "g", "x2", 1.5)
# IMPORTANT: if nonzero, the x1->x2 transition will be part of the calculation
dipole!(s, "x1", "x2", 0)
# read transition dipole
@assert dipole(s, "g", "x1") == dipole(s, "x1", "g")
```
The lineshape functions are set in a similar manner. The calculation should work
with any object callable with a single argument. This can be a function (such as
`f` below), an anonymous function (in Julia: `t -> 0.3*t`) or a callable
object (interpolation, integrator, ...). In order to improve performance, you
can use a look-up table (`LineshapeLUT`) object to cache the value of the
lineshape function.
```julia
# set lineshape functions
f(t) = 0.1*t
lineshape!(s, "x1", "x1", f) # accept any callable.
lineshape!(s, "x2", "x2", t->g_homo(t, 0.1)) # anonymous function in Julia: x->...
lineshape!(s, "x1", "x2", zero) # uncorrelated.
lineshape!(s, "x2", "x1", zero) # uncorrelated.
# read lineshape function
@assert lineshape(s, "x1", "x1") === f
```

The spectroscopic response is calculated in the time domain. You can use the
`TimeGrid` object to define the domain of calculation. The calculation can be performed automatically for all points of the grid and 
for all paths. 
```julia
t1=0:20.0 # 0 to 20 included in steps of 1.0
t2=0:10.0:100.0 # 0 to 100 included in steps of 10.0
t3=0:20.0
tg = TimeGrid(t1, t2, t3)
# compute everyting
rtot = R1(tg, s) + R2(tg, s) + R3(tg, s) + R4(tg, s)
```
Alternatively, you can generate Hilbert paths and compute
for specific paths manually.
```julia
hpaths = collect(hilbert_paths(s, 3)) # third order hilbert paths
@assert length(hpaths) == 4
# compute for a specific path, pick the first.
p = hpaths[1]
r1 = R1(tg, s, p)
```
This also works for the linear reponse:
```julia
rlin = linear(tg, s) # use the first time argument.
lin_paths = collect(hilbert_paths(s, 1))
lin_resp = [linear(tg, s, p) for p=lin_paths] # individual contributions.
@assert all(rlin .== sum(lin_resp))
```

Features
--------
Currently supported:
- Linear response in the time domain.
- Third order response.
- Set energies, transition dipoles (real only) and lineshape functions.
- Computation for all pathways, or a subset.
- Arbitrary callables for lineshape functions.
- Built-in lineshapes: Homogeneous, Inhomogeneous, Huang-Rhys, Kubo ansatz.
- Lineshape functions from correlation functions.
- Caching lineshape functions for improved performance.
- Unit conversion: eV -> angular PHz (not automatic.)
- IO using julia: delimited files, YAML, binary, memory-mapped.

Works, but deserve improvement:
- Rotating frames.
- Lineshapes from spectral densities.
- Filtering of pathways.

Citing
======
If you found this useful, please consider citing the following paper:
> [WIP!]

Examples
========
WIP
<!--TODO: do this. 
- two level system (YAML)
- two level system with phonon
- three level system, g e1 e2
- three level system, g e f (ladder)
- Example coherence beatmap (g e1 e2 f)
- Example vibrational beatmap (g e f) + phonon
- Example: correlation functions and many body.
-->

Scripts
-------
Some examples use external scripts for pre or post processing (eg: plotting in
Julia still requires work). As such, some examples are not fully self contained
and can involve scripts written in `python`, `julia` and `powershell`.
<!--TODO: COMPLETE THIS -->

Help!
=====
Please open an issue using the 'Issues' tab at the top. Anything works,
including if you're lost due to missing documentation. 

TODO
====
- [x] Handle ground states easily. (Default lineshape functions in case they're missing)
- [ ] Compilations to LUT
- [ ] Add rotating frames properly.
- [ ] Add automatic rephasing vs non-rephasing (`rephasing(grid,system)`).
- [ ] Convenience for zeroing transition dipole moments (`forbidden`?)
- [ ] Convenience for uncorrelating states (lineshape function=zero)?
<!-- 
- [ ] More documentation...
- [ ] Add filtering of Hilbert Paths (eg: IA only).-->
<!--
## May happen...
- Explicit bath modes
-->

References
==========

1. <a name="book_cho"></a>Cho, M. *Two-dimensional optical spectroscopy.* (CRC Press, 2009).
2. <a name="book_mukamel"></a>Mukamel, S. *Principles of nonlinear optical spectroscopy.* (Oxford University Press, 1995).