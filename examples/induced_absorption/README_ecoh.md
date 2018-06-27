# Simple electronic coherence example

The calculation is made for a simple three-level system (two singly excited 
states A,B), the simplest system which can yield an electronic coherence.

This simple example demonstrates the properties of electronic coherences.

The cumulant expansion as used here requires 4 lineshape functions: `g_AA`,
`g_BB`, `g_AB` and `g_BA`. We are assuming linear coupling to uncorrelated
classical bath modes, such that `g_BA(t) = g_AB(t)`. This reduces the number of
independant lineshape functions to 3. Similarly to the `basic` example, each
lineshape function is taken in the same limit of pure homogeneous and pure
inhomogeneous broadening and thus requires 2 parameters. In the description
used here, perfect correlation is achieved when `g_AB(t) = sqrt(g_AA(t) g_BB(t))`

To run the calculation, run:
```
julia ecoh.jl <config>
```
For a description of the output files, see the simple example.

There is an extra plotting script: `plot_dyn.py` which plots the dynamics of
regions around the peaks, revealing the electronic coherence.

A utility script to run the calculation and plotting is supplied as
`run_cfg.ps1`. It should be trivial to translate to bash.
