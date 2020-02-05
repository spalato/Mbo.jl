#!python
# make plots, in python
# usage: python plot_2d.py <config>

import sys
import os.path as pth
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import yaml
from scipy.constants import eV, h, c, femto
from numpy.fft import fft, ifft, fftshift, fftfreq, ifftshift
try:
    from styles import presentation
except ImportError:
    pass
else:
    presentation()

def phz2ev(x):
    scale = h/eV/femto
    return x*scale

def expand_axis(x):
    """Expand axis by first value, for use with pcolormesh"""
    new = np.zeros(x.shape)
    new[1:] = 0.5*(x[1:]+x[:-1])
    first = 2*new[1]-new[2]
    new[0] = first
    return new

def setup_time_axis(idx, cfg):
    try:
        ax = np.linspace(0, cfg[f"t{idx}_max"], cfg[f"t{idx}_n"])
    except KeyError:
        ax = np.array(cfg[f"t{idx}"])
    return ax

cfgf = sys.argv[1]
fld = pth.dirname(cfgf)
with open(cfgf, "r") as f:
    cfg = yaml.load(f)

locals().update(**cfg) # that's how we roll.

t_1 = setup_time_axis(1, cfg)
t_2 = setup_time_axis(2, cfg)
t_3 = setup_time_axis(3, cfg)
t1_n = t_1.size
t2_n = t_2.size
t3_n = t_3.size

for tresp in ["rr", "rn"]:
    fn = "{}_{}.bin".format(rootname, tresp)
    sig = np.memmap(pth.join(fld, fn), dtype=np.complex128,
                    shape=(t1_n, t2_n, t3_n), order="F")
    plt.figure(figsize=(4,4), dpi=300)
    z = np.real(sig[:,0,:])
    lim = np.max(np.abs(z))
    plt.pcolormesh(expand_axis(t_1), expand_axis(t_3), z.T, cmap=plt.cm.seismic, vmin=-lim, vmax=lim)
    levels = np.linspace(-lim, lim, 11)
    plt.contour(t_1, t_3, z.T, levels=levels, colors='0.8', linewidths=0.5)
    outname = fn[:-4]+".png"
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.savefig(outname)

f_1 = phz2ev(fftshift(fftfreq(t_1.size, t_1[1]-t_1[0])))
#f_2 = phz2ev(fftshift(fftfreq(t_2.size, t_2[1]-t_2[0])))
f_3 = phz2ev(fftshift(fftfreq(t_3.size, t_3[1]-t_3[0])))

for fresp in ["sr", "sn", "sa"]:
    fn = "{}_{}.bin".format(rootname, fresp)
    sig = np.memmap(pth.join(fld, fn), dtype=np.complex128,
                    shape=(t1_n, t2_n, t3_n), order="F")
    plt.figure(figsize=(4,4), dpi=300)
    z = np.real(sig[:,0,:])
    lim = np.max(np.abs(z))
    plt.pcolormesh(expand_axis(f_1), expand_axis(f_3), z.T, cmap=plt.cm.seismic, vmin=-lim, vmax=lim)
    levels = np.linspace(-lim, lim, 11)
    plt.contour(f_1, f_3, z.T, levels=levels, colors='0.8', linewidths=0.5)
    outname = fn[:-4]+".png"
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.savefig(outname)
