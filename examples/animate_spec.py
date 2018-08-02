#!python3
#animate_spec: make movie of a 2d spectra vs popt
# usage: python animate_spec.py <config> <in_file> <out_file>

import sys
import os.path as pth
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.animation import FFMpegWriter
import yaml
from scipy.constants import eV, h, c, femto
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(it):
        print("Rendering...")
        return it
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

cfgf, in_f, out_f = sys.argv[1:4]

metadata = {
    "title": in_f,
    "artist": "Mbo.jl",
}
writer = FFMpegWriter(fps=10, metadata=metadata)

fld = pth.dirname(cfgf)
with open(cfgf) as f:
    cfg = yaml.load(f)

t_1 = np.linspace(0, cfg["t1_max"], cfg["t1_n"])
t_2 = np.linspace(0, cfg["t2_max"], cfg["t2_n"])
t_3 = np.linspace(0, cfg["t3_max"], cfg["t3_n"])

f_1 = phz2ev(fftshift(fftfreq(t_1.size, t_1[1]-t_1[0])))
f_3 = phz2ev(fftshift(fftfreq(t_3.size, t_3[1]-t_3[0])))

spec = np.memmap(in_f,
                  dtype=np.complex128,
                  shape=(cfg["t1_n"], cfg["t2_n"], cfg["t3_n"]),
                  order="F")

fig = plt.figure(figsize=(4.8, 4.0), dpi=300)

z = np.real(spec)

lim = np.max(np.abs(z))
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

u = max([np.max(f) for f in (f_3, f_1)])
l = min([np.min(f) for f in (f_3, f_1)])
levels = np.linspace(-lim, lim, 11)

with writer.saving(fig, out_f, 300):
    for i, t2_ in tqdm(enumerate(t_2), total=t_2.size):
        #print(x.size, y.size, z[:,i,:].shape)
        fig.clear()
        plt.pcolormesh(expand_axis(f_1), expand_axis(f_3), np.real(z[:,i,:]).T, 
                cmap="seismic", vmin=-lim, vmax=lim)
        plt.contour(f_1, f_3, np.real(z[:,i,:]).T,
            levels=levels, colors='0.8', linewidths=1)

        plt.xlabel(r"$\Delta E_1$ (meV)")
        plt.ylabel(r"$\Delta E_3$ (meV)")
        plt.text(0.05, 0.95, r"$t_2$"+": {:4.0f} fs".format(t2_),
            horizontalalignment="left",
            verticalalignment="top",
            transform=plt.gca().transAxes)
        plt.gca().xaxis.set_major_locator(plt.MultipleLocator(0.2))
        plt.gca().yaxis.set_major_locator(plt.MultipleLocator(0.2))
        plt.xlim(np.min(f_1), np.max(f_1))
        plt.ylim(np.min(f_3), np.max(f_3))
        plt.tight_layout()
        writer.grab_frame()
        #if i >= 10: break
