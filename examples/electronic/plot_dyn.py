# plot peak dynamics
# usage: python plot_dyn.py <config>

import sys
import os.path as pth
from itertools import product, cycle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import yaml
from numpy.fft import fft, ifft, fftshift, fftfreq, ifftshift
from scipy.constants import eV, h, c, femto
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from styles import presentation
presentation()

def between(x, lo, hi):
    return (x>lo) & (x<hi)

def slicemask(x, loc, width):
    hw = width/2
    return between(x, loc-hw, loc+hw)

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

cfgf = sys.argv[1]
fld = pth.dirname(cfgf)
with open(cfgf, "r") as f:
    cfg = yaml.load(f)

locals().update(**cfg) # that's how we roll.

t_1 = np.linspace(0, cfg["t1_max"], cfg["t1_n"])
t_2 = np.linspace(0, cfg["t2_max"], cfg["t2_n"])
t_3 = np.linspace(0, cfg["t3_max"], cfg["t3_n"])
f_1 = phz2ev(fftshift(fftfreq(t_1.size, t_1[1]-t_1[0])))
f_3 = phz2ev(fftshift(fftfreq(t_3.size, t_3[1]-t_3[0])))

fn = "{}_sa.bin".format(rootname)
sig = np.memmap(pth.join(fld, fn), dtype=np.complex128,
                shape=(t1_n, t2_n, t3_n), order="F")

tags = "AB"
locs = [e_A-e_frame, e_B-e_frame]
sl_w = 0.08
sl_hw = sl_w/2
fig = plt.figure()
ax = plt.gca()
ax_inset = inset_axes(ax, "25%", "30%", 4)
lstyles = cycle(["-", "--"])
for (t1, t3), (l1, l3) in zip(product(tags, tags), product(locs, locs)):
    sl = sig.compress(slicemask(f_1, l1, sl_w), axis=0).compress(slicemask(f_3, l3, sl_w), axis=2)
    sl = sl.sum(axis=(0,2))
    #sl = sig[argnear(f_1, l1),:,argnear(f_3,l3)]
    l = ax.plot(t_2, np.real(sl), label="{},{}".format(t1, t3),
                alpha=0.7, ls=next(lstyles))
    color = l[0].get_color()
    #ax_inset.axvline(l1, color=color, lw=0.5, alpha=0.5)
    #ax_inset.axhline(l3, color=color, lw=0.5, alpha=0.5)
    ax_inset.add_patch(plt.Rectangle(
        (l1-sl_hw, l3-sl_hw), sl_w, sl_w, facecolor=color, alpha=0.5, edgecolor='k'
        ))
#ax.legend()

z = np.real(sig[:,0,:])
vlim = np.max(np.abs(z))
ax_inset.pcolormesh(expand_axis(f_1), expand_axis(f_3), z.T, cmap="seismic", vmin=-vlim, vmax=vlim, zorder=-10)
ax_inset.tick_params(labelleft=False, labelbottom=False)
margin = 0.05
bounds = [locs[0]-margin, locs[1]+margin]
ax_inset.set_xlim(*bounds)
ax_inset.set_ylim(*bounds)
outname = fn[:-4]+"_peakdyn.png"
plt.savefig(outname, dpi=300)
    

