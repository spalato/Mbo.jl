# Simple vibrational coupling example

The calculation is made for a simple two-level system coupled to a vibrational
mode, including pure homogeneous and inhomogeneous broadening.

This simple example demonstrates the effect of vibrational modes. It has 3
extra parameters compared to the `simple` example: Huang-Rhys parameter S,
vibrational mode energy and temperature. Two examples parameters are supplied:
- `vibrational.yaml` for small value of S (0.2). This looks like coupled transitions.
- `wavepacket.yaml` for large value of S (2.0). This looks like wavepacket dynamics.
Have fun playing with the parameters.

To run the calculation, run:
```
julia vibrational.jl <config>
```
The time grid is much larger than in the simple case. For a description of the
output files, see the simple example.
