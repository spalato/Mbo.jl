Mbo.jl - Multimode Brownian Oscillator
======================================

Semiclassical modeling of non-linear spectra in the time domain using the
Multimode Brownian Oscillator model. Designed for the modeling of 2D
visible spectra. Hopefully fast.

Description
-----------

This code aims to put the modeling of the 2D spectra of simple
model systems within reach of the humble experimental spectroscopist. 
It is a toolkit built around the spectroscopic system: the
user specifies state energies, transition dipole moments and
lineshape functions. The linear and non-linear responses can then be
conveniently calculated from the system object. 

For example:
```julia
# build a system. It's still a bit tedious, but that's where the thinking happens.
using Mbo
s = System("g") # use "g" as a ground state
# Add a state of energy of state "e" is 1.5 eV, which we need to convert
# to angular PHz
energy!(s, "e", ev2angphz(1.5))  
# Set transition dipole moment of the g->e transition
dipole!(s, "g", "e", 1.0)

# Define a lineshape functions. Here, 0.02 t + (0.02 t)^2/2
# Most common cases are supplied. 
ls(t) = ev2angphz(0.02)*t + 0.5*(ev2angphz(0.02)*t)^2
# Set the function `ls(t)` as the lineshape for the e-e response.
lineshape!(s, "e", "e", ls)

# calculation grid, in fs. tedious bookkeeping done by the code.
tg = TimeGrid(0.0:100, 0.0:100, 0.0:100)
lin = linear(tg, s) # linear response
rr = R2(tg, s) + R3(tg, s) # rephasing
rn = R1(tg, s) + R4(tg, s) # non-rephasing
# Go ahead and process to your wishes.
```
The example above uses hand-coded constants, but you are free to use any input
or output format you desire. Same goes for lineshape functions: they can be
defined arbitrarily using the programming language (some common forms are
predefined). This allows for flexibility: the system's parameters can equally be
hard coded, read from a file, or computed from ancillary algorithms (for
example, diagonalizing a site basis to an exciton basis). `Julia` can call
code written in `C`, `Fortran` and `python`.

The best way to get started using the code for modeling is to look at supplied
examples. See the `/examples/` folder and the [Examples section](#Examples)
below.

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
This packages uses the [julia](https://julialang.org/) programming language.
You'll need it, but it is easy to install.

`Julia` is an  open-source language. It aims to provide clear syntax, similar to 
`Matlab` or `python` while providing performance typical of compiled languages
such as `Fortran` or `C`. Hopefully, the language is easy to pick up for anyone
familiar with any of these. Julia is available for Windows, Mac, Linux. Get it
from [here](https://julialang.org/downloads/). 

Once julia is installed, open `julia` and type:

```julia
julia> ]
```
Then 
```julia
@v1.10 pkg> add https://github.com/spalato/Mbo.jl
```
If you get an error regarding version of `git`, try to run `julia` after
deactivating all your `conda` envs. If you run into issues, please file an issue
in `github` using the `Issues` tab above.

<!--
```julia
julia> ]
```

Then 
```julia
@v1.10 pkg> activate /path/to this repository/
```
Your terminal should then look like

```julia
Mbo.jl pkg>
```
Finally precompile the package with 

```julia
Mbo.jl pkg> precompile
```
-->

The initial version of this code (up to `Mbo.jl v1.0`) was written for 
`julia v0.6`, which was never intended for long-term support. The current 
version of `Mbo.jl` has been updated to work with `julia v1.10`. For upgrading
your scripts, see the [migration guide](migration_guide.md). There really isn't
much to do.

Documentation
============
The main element of this package is the `System` object, which contains all the
physical parameters of the system under study. Using the cumulant expansion,
this means state energies, transition dipole moments and lineshape functions.

Physics
-------
The third order response is calculated using the Cumulant expansion to second
order. This method is briefly presented here. For details, refer to the books
by [Mukamel](#book_mukamel)<!-- !! add chapters!! --> and [Cho](#book_cho).

The calculation of 2D spectra is carried out using the 4-point correlation
function of the transition dipole moments. The 4-point correlation function
has 4 time arguments, whose permutations give rise to the more commonly
discussed double-sided Feynman diagrams.

When multiple states are involved, different transition dipole operators can be
involved in the third order response, and thus in the 4-point correlation
function. A given combination of 4 transition dipole operator
take the system through up to 4 spectroscopically coupled states. The path the 
system takes through its manifold of states is called here a Hilbert Path.
As previously mentionned, each of these Hilbert paths gives rise to 4 
double-sided Feynmann diagrams.

The calculation of the third order response is tedious: a given third order 
response function requires 3 time arguments, 3 frequencies, 6 lineshape
functions and 4 transition dipoles. This has to be repeated for every possible
path the system can take in Hilbert space. For `N` states (+ a single ground 
state), there are around `N^2` such paths. This tedium is taken care of
using the `System` and `HilbertPath` objects. Handling of the time arguments
is simplified using the `TimeGrid` object.

The equations for third order responses R<sub>1</sub> to R<sub>4</sub> were
taken from
the book by [Cho](#book_cho), eq 5.26-5.31. These equations take the fluctuating
ground state as a reference energy. The lineshape functions are from the
fluctuations of the energies of
the states (not the energy gaps). As a consequence, fluctuations of the ground
state are always 0. The code is currently limited to a single ground state and
real transition dipoles (there are likely simple fixes for this).

Basics
------
The attributes of the system are set and read using functions. Functions that
modify the `System` object end with `!` (this is purely a `julia` convention,
`!` has no special meaning).  States are indexed by case-sensitive strings.
You can use "G", "X1", "S+3/2L", as you wish.

Setting and reading state energies and transition dipoles is rather
straightforward. No implicit unit conversions are made, it's up to the user
to ensure the proper units are supplied. The code expects energies
to be supplied as angular frequencies, inverse of your time units. If you use a
time grid in fs, you should supply energies in angular PHz.
```julia
using Mbo
# use "g" as a ground state with energy 0
s = System("g")
@assert "g" in states(s) # use the @assert macro for simple tests
@assert energy(s, "g") == 0
# set energies in angular frequencies, inverse of your time axis (PHz for fs)
energy!(s, "x1", 1.2) # x1 has ω = 1.2 ang PHz
energy!(s, "x2", 1.4) # x2 has ω = 1.4 ang PHz
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
Lineshape functions are handled similarly. You can supply
 any object callable with a single argument. This can be a function (`f` 
below), an anonymous function (in Julia: `t -> 0.3*t`) or a callable
object (interpolation, integrator, ...). In order to improve performance, you
can use a look-up table (`LineshapeLUT`) object to cache the value of the
lineshape function. See examples.
```julia
# set lineshape functions
f(t) = 0.1*t
lineshape!(s, "x1", "x1", f) # accept any callable.
lineshape!(s, "x2", "x2", t->g_homo(t, 0.1)) # anonymous function
lineshape!(s, "x1", "x2", zero) # uncorrelated.
lineshape!(s, "x2", "x1", zero) # uncorrelated.
# read lineshape function
@assert lineshape(s, "x1", "x1") === f
```

The
`TimeGrid` object defines the domain of calculation. By default, the 
calculation is performed automatically for all points of the grid and 
for all paths.
```julia
t1=0:20.0 # 0 to 20 included in steps of 1.0
t2=0:10.0:100.0 # 0 to 100 included in steps of 10.0
t3=0:20.0
tg = TimeGrid(t1, t2, t3)
# compute everyting
rtot = R1(tg, s) + R2(tg, s) + R3(tg, s) + R4(tg, s)
```
Alternatively, you can generate Hilbert paths and limit the calculation
to a subset of paths.
```julia
hpaths = collect(hilbert_paths(s, 3)) # third order hilbert paths
@assert length(hpaths) == 4 # 4 paths x 4 reponses -> 16 Feynmann diagrams.
# compute for a specific path, pick the first.
p = hpaths[1]
r1 = R1(tg, s, p)
```
The same `System` and `TimeGrid` can be used to compute the linear reponse.
```julia
rlin = linear(tg, s) # uses the first time delay
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
- Built-in lineshapes: Homogeneous, Inhomogeneous, Huang-Rhys, Kubo ansatz and variations thereof.
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
> S. Palato, H. Seiler, P. Nijjar, O. Prezhdo, P. Kambhampati, Atomic fluctuations in electronic materials revealed by dephasing, *PNAS* (2020) https://doi.org/10.1073/pnas.1916792117

Examples
========

Examples are located in the `examples/` directory. The parameters of the
calculations (energies, lineshape parameters) are stored in a separate `.yaml`
config file. Plotting is carried out separately from the main calculation, see
[Scripts](#Scripts).

- basic: Simple two-level system. Basic functionnality, plotting in julia and python.
- vibrational: Two-level system coupled to a coherent phonon without damping.
- electronic: Three-level system (two singly excited states) demonstrating electronic coherences.
- induced absorption: Three-level system featuring a doubly excited state.
<!--TODO: do this. 
- Example coherence beatmap (g e1 e2 f)
- Example vibrational beatmap (g e f) + phonon
- Example: correlation functions and many body
- Make sure to demonstrate LineshapeLUT
-->

Scripts
-------
Some examples use external scripts for pre or post processing (eg: plotting in
Julia still requires work). As such, some examples are not fully self contained
and can involve scripts written in `python`, `julia` and `powershell`. These
scripts can be found in the `examples/` directory.

- `plot_linear.jl`: plot linear response (time and frequency domain) in julia;
- `plot_linear.py`: plot linear response in python;
- `plot_2d.jl`: plot 3rd order response (t,t,t and f,t,f) in julia;
- `plot_2d.py`: plot 3rd order response (t,t,t and f,t,f) in python.

Example plotting scripts are supplied for `julia` and `python`.
Plotting in julia requires a few extra packages, listed in `examples/julia_requires`.
Plotting in python requires `numpy` and `matplotlib`, see `examples/python_requirements.txt`.

Help!
=====
Please open an issue using the 'Issues' tab at the top. You found a bug? You
can't achieve what you want? Lost and confused? Maybe we can help (maybe not). 

<!-->
TODO
====
- [ ] Update python scripts to remove deprecation warnings
- [ ] Tidy up python dependencies
- [ ] Tidy up julia dependencies

<!--
- [ ] Switch to TOML
- [ ] Add rotating frames properly.
- [ ] Add automatic rephasing vs non-rephasing (`rephasing(grid,system)`).
- [ ] Do we want to add filters? Maybe just a `cookbook` is ok.
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
3. <a name="book_hammzanni"></a>Hamm, P. & Zanni, M. *Concepts and methods of 2D infrared spectroscopy.* (Cambridge University Press, 2011).
