# Simple induced absorption

The calculation is made for a simple three-level system: ground state `g`,
excited state `e` and doubly excited state `f`. This is the simplest system to
yield induced absorption. Again, pure homogeneous and pure inhomogeneous
dephasing are used. In most physical models, the energy of states `e` and `f`
should show some degree of correlation.

The example supplied here is for a system with a constant binding energy `d`
(`E_f = 2 E_e - d`) and perfectly correlated fluctuations of the energy levels.
I strongly suggest trying out different lineshape parameters and 
making sense of the results. Welcome to the rabbit hole.

To run the calculation, run:
```
julia basic_ia.jl <config>
```
For a description of the output files, see the `basic` example. This model has
no dynamics, so it has a single point along `t_2`.