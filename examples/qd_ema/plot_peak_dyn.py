#!python
# compute and plot the beatmap
# usage: python plot_beatmap.py <config>

import sys
import sys
import os.path as pth
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import yaml
from scipy.constants import eV, h, c, femto
from numpy.fft import fft, ifft, fftshift, fftfreq, ifftshift
# try:
#     from styles import presentation
# except ImportError:
#     pass
# else:
#     presentation()

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

def between(a, l, u):
    """Boolean mask for values in `a` larger than `l`, smaller than `u`"""
    return (a > l) & (a < u)

def scale_to_last(arr, n=1, func=np.median):
    return arr/np.abs(func(arr[-n:]))


cfgf = sys.argv[1]
fld = pth.dirname(cfgf)
with open(cfgf, "r") as f:
    cfg = yaml.load(f)

locals().update(**cfg) # yep

t_1 = np.linspace(0, cfg["t1_max"], cfg["t1_n"])
t_2 = np.linspace(0, cfg["t2_max"], cfg["t2_n"])
t_3 = np.linspace(0, cfg["t3_max"], cfg["t3_n"])

f_1 = phz2ev(fftshift(fftfreq(t_1.size, t_1[1]-t_1[0])))+cfg["e_frame"]
f_2 = phz2ev(fftshift(fftfreq(t_2.size, t_2[1]-t_2[0])))+cfg["e_frame"]
f_3 = phz2ev(fftshift(fftfreq(t_3.size, t_3[1]-t_3[0])))+cfg["e_frame"]


fn = "{}_sa.bin".format(rootname)
sig = np.memmap(pth.join(fld, fn), dtype=np.complex128,
                shape=(t1_n, t2_n, t3_n), order="F", mode="c")

roi_hw = 0.03
roi = (cfg["e_x2"], cfg["e_x1"]-cfg["binding"]/2)
e1_mask = between(f_1, roi[0]-roi_hw, roi[0]+roi_hw)
e3_mask = between(f_3, roi[1]-roi_hw, roi[1]+roi_hw)
xpeak = np.real(sig).compress(e1_mask, axis=0).compress(e3_mask, axis=2).sum(axis=(0,2))

fig = plt.figure(dpi=150, figsize=(6, 2.5))
gs = plt.GridSpec(1, 2, width_ratios=[2,3])
axes = [
    plt.subplot(gs[0]),
    plt.subplot(gs[1])
]
plt.sca(axes[0])
z = np.real(sig[:,0,:])
lim = np.max(np.abs(z))
z = z / lim
plt.pcolormesh(expand_axis(f_1), expand_axis(f_3), z.T, vmin=-1, vmax=1, cmap="seismic")
plt.contour(f_1, f_3, z.T, levels=np.linspace(-1, 1, 11), colors="0.9", linewidths=0.5)
plt.xlim(1.85, 2.1)
plt.ylim(1.85, 2.1)
plt.gca().add_patch(plt.Rectangle(
    [c-roi_hw for c in roi], 2*roi_hw, 2*roi_hw, fill=False, zorder=10
))
plt.xlabel("E$_1$ (eV)")
plt.ylabel("E$_3$ (eV)")

plt.sca(axes[1])
plt.plot(t_2, -scale_to_last(xpeak))
plt.xlabel("t$_2$ (fs)")
plt.ylabel("Signal (arb.)")
plt.tight_layout()
plt.savefig(f"{rootname}_xpeak_dyn.png")
