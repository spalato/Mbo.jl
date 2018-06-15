#!python
# make plots, in python
# usage: python plot_2d.py

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import yaml
from numpy.fft import fft, ifft, fftshift, fftfreq, ifftshift


def expand_axis(x):
    """Expand axis by first value, for use with pcolormesh"""
    new = np.zeros(x.shape)
    new[1:] = 0.5*(x[1:]+x[:-1])
    first = 2*new[1]-new[2]
    new[0] = first
    return new

cfgf = sys.argv[1]

with open(cfgf, "r") as f:
    cfg = yaml.load(f)

locals().update(**cfg) # that's how we roll.

t_1 = np.linspace(0, cfg["t1_max"], cfg["t1_n"])
t_2 = np.linspace(0, cfg["t2_max"], cfg["t2_n"])
t_3 = np.linspace(0, cfg["t3_max"], cfg["t3_n"])

for tresp in ["rr", "rn"]:
    fn = "{}_{}.bin".format(rootname, tresp)
    sig = np.memmap(fn, dtype=np.complex128,
                    shape=(t1_n, t2_n, t3_n), order="F")
    plt.figure(figsize=(4,4), dpi=300)
    z = np.real(sig[:,0,:])
    lim = np.max(np.abs(z))
    plt.pcolormesh(t_1, t_3, z.T, cmap=plt.cm.seismic, vmin=-lim, vmax=lim)
    levels = np.linspace(-lim, lim, 11)
    plt.contour(t_1, t_3, z.T, levels=levels, colors='0.8', linewidths=0.5)
    outname = fn[:-4]+".png"
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.savefig(outname)

f_1 = fftshift(fftfreq(t_1.size, t_1[1]-t_1[0]))
#f_2 = phz2ev(fftshift(fftfreq(t_2.size, t_2[1]-t_2[0])))
f_3 = fftshift(fftfreq(t_3.size, t_3[1]-t_3[0]))

for fresp in ["sr", "sn", "sa"]:
    fn = "{}_{}.bin".format(rootname, fresp)
    sig = np.memmap(fn, dtype=np.complex128,
                    shape=(t1_n, t2_n, t3_n), order="F")
    plt.figure(figsize=(4,4), dpi=300)
    z = np.real(sig[:,0,:])
    lim = np.max(np.abs(z))
    plt.pcolormesh(expand_axis(f_1), expand_axis(f_3), z.T, cmap=plt.cm.seismic, vmin=-lim, vmax=lim)
    levels = np.linspace(-lim, lim, 11)
    plt.contour(expand_axis(f_1), expand_axis(f_3), z.T, levels=levels, colors='0.8', linewidths=0.5)
    outname = fn[:-4]+".png"
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.savefig(outname)
