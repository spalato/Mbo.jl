#!python3
# beatmap.py: extract beatmaps
# usage: python beatmap.py <config>

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import yaml
from styles import presentation
from spectro.utils import phz2ev
from numpy.fft import fft, ifft, fftshift, fftfreq, ifftshift
from analysis2d.plot_utils import plot_2d

g_kw = {
    "lw": 0.5,
    "c":"k",
    "alpha": 1,
    "ls":":",
}


cfgf = sys.argv[1]

with open(cfgf, "r") as f:
    cfg = yaml.load(f)

locals().update(**cfg)
inf = rootname+"_sa.bin"

outroot = inf[:-4]
s_r = np.memmap(inf, dtype=np.complex128,
                shape=(t1_n, t2_n, t3_n), order="F")
t_1 = np.linspace(0, cfg["t1_max"], cfg["t1_n"])
t_2 = np.linspace(0, cfg["t2_max"], cfg["t2_n"])
t_3 = np.linspace(0, cfg["t3_max"], cfg["t3_n"])
f_1 = phz2ev(fftshift(fftfreq(t_1.size, t_1[1]-t_1[0])))
f_2 = phz2ev(fftshift(fftfreq(t_2.size, t_2[1]-t_2[0])))
f_3 = phz2ev(fftshift(fftfreq(t_3.size, t_3[1]-t_3[0])))

resid = s_r - np.mean(s_r, axis=1)[:,np.newaxis,:]
bm_real = fftshift(ifft(np.real(resid), axis=1), axes=1)

n_tones = 1
tones = [i for i in range(1, n_tones+1) if i != 0]
locs = [np.argmin(np.abs(f_2-e_vib*i)) for i in tones]
# total intensity spectrum
z = np.sum(np.abs(bm_real), axis=(0,2))
#m_ = f_2>0
fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(4.6, 4.6))
for ax in axes:
    plt.sca(ax)
    plt.plot(f_2[f_2!=0], z[f_2!=0], 'k.-')
    for i, l in zip(tones, locs):
        plt.plot(f_2[l], z[l], 'r+', mew=1)
    plt.xlim(0, 0.155)
plt.sca(axes[0])
plt.ylabel("Tot. int. (abs)")
plt.gca().yaxis.set_label_coords(-0.12, 0)
plt.sca(axes[1])
plt.xlabel(r"$E_{coh}$")
plt.yscale("log")
plt.tight_layout()
plt.subplots_adjust(hspace=0.1, top=0.95, bottom=0.12)
plt.savefig(outroot+"_cohmap_totint.png", dpi=300)
plt.close()
# maps


selected = [bm_real[:,l,:] for l in locs]
lim = max([np.max(np.abs(s)) for s in selected])
levels = np.linspace(0, lim, 11)
for i, sel in zip(tones, selected):
    plt.figure(figsize=(4.8, 4.0))
    plot_2d(f_1, f_3, np.abs(sel).T,
        cmap='YlGn',
        symz=False,
        colors='k',
        diag="k",
        levels=levels,
        pcolor_kw={"vmin":0,},
        )
    e0 = e_e
    e_0 = e0-e_frame
    for l in [e_0-e_vib, e_0, e_0+e_vib]:
        plt.axhline(l, **g_kw)
        plt.axvline(l, **g_kw)
    txt = plt.text(0.05, 0.95,
        "Max: {:.03g}".format(np.max(np.abs(sel))),
        transform=plt.gca().transAxes, va='top',
        fontsize='small', color='k', zorder=10)
    #txt.set_path_effects([pe.withStroke(linewidth=2, foreground='k')])
    plt.xlim(-0.15, 0.15)
    plt.ylim(-0.15, 0.15)
    plt.xlabel(r"$E_1$ (eV)")
    plt.ylabel(r"$E_3$ (eV)")
    plt.tight_layout()
    plt.savefig(outroot+"_cohmap_amp_{}.png".format(i), dpi=300)

# phase maps
for i, sel in zip(tones, selected):
    fig = plt.figure(figsize=(5.4, 4.0))
    phase = 0.5*(np.angle(sel)/np.pi+1)
    colors = plt.cm.hsv(phase.T)
    alphas = np.abs(sel.T)
    alphas /= np.max(alphas)
    alphas /= 0.8
    alphas = np.clip(alphas, 0, 1)
    colors[:,:,-1] = alphas
    plt.contour(f_1, f_3, np.abs(sel.T),
                levels=levels, colors='k', linewidths=1)
    im = plt.imshow(colors,
                    extent=(np.min(f_1), np.max(f_1), np.min(f_3), np.max(f_3)),
                    origin="bottom",
                    cmap='hsv', vmin=-1, vmax=1,
                    aspect='auto',
            )
    txt = plt.text(0.05, 0.95,
        "Max: {:.03g}".format(np.max(np.abs(sel))),
        transform=plt.gca().transAxes, va='top',
        fontsize='small', color='k', zorder=10)
    #cb.ax.set_yticklabels([r"$\pi$", ])
    plt.plot([-1, 1], [-1,1], 'k:', lw=1)
    for feat in [e_0+e_vib*i for i in [-1,0,1]]:
        plt.axhline(feat, lw=0.5, color='0.5')
        plt.axvline(feat, lw=0.5, color='0.5')
    #for feat in [-delta, -delta-de, -delta+de]:
    #    plt.axhline(feat, lw=0.5, ls="--", color='0.5')

    plt.xlim(-0.15, 0.15)
    plt.ylim(-0.15, 0.15)
    plt.xlabel(r"$E_1$ (eV)")
    plt.ylabel(r"$E_3$ (eV)")

    ax = plt.gca()
    cb = fig.colorbar(im, ticks=np.linspace(-1,1,9), ax=ax)
    cb.ax.set_yticklabels([-1, "", "", "", 0, "", "", "", 1])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize="small")
    #plt.subplots_adjust(left=0.2)
    plt.tight_layout()
    plt.savefig(outroot+"_cohmap_phi_{}.png".format(i), dpi=300)
