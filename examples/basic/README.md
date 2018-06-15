# Simple example

The calculation is made for a simple two-level system with pure homogeneous and
inhomogeneous broadening.

This simple example demonstrates basic functionality of the `Mbo` module.
It also serves to introduce julia syntax and python interoperability.
Information about the system is contained in a parameter file called
`basic.yaml`. See comments in the file for explanations.

To run the calculation, run:
```
julia basic.jl basic.yaml
```

This will produce the following output files:
- `basic_rlin.txt`: linear response in the time domain. 3 columns text: time, real part, imaginary part;
- `basic_slin.txt`: linear spectrum in the frequency domain. 3 columns text: frequency (PHz), real part, imaginary part;
- `basic_rr.bin`: rephasing reponse in the time domain, binary;
- `basic_rn.bin`: non-rephasing reponse in the time domain, binary;
- `basic_sr.bin`: rephasing spectrum (f1, t2, f3), binary;
- `basic_sn.bin`: non-rephasing spectrum (f1, t2, f3), binary;
- `basic_sa.bin`: fully absorptive spectrum (f1, t2, f3), binary.

Binary files contain a three dimensional array of complex-128 numbers in Fortran order.
They are produced using julia's [`write`](https://docs.julialang.org/en/stable/stdlib/io-network/#Base.write) function. For how to load in python, see `plot_2d.py`

Plotting scripts are available in both julia and python. Using julia plotting
doesn't require a working python install, but I find it horribly slow (it was
probably made to be web-scale...). Python scripts also demonstrate how to load
binary data in python.