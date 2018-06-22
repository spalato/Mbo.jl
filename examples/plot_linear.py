#!python
# plot_linear.py: plot linear response in the time and frequency domains.
# usage: python plot_linear.py <config>

import sys
import os.path as pth
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import yaml

cfgf = sys.argv[1]

with open(cfgf, "r") as f:
    cfg = yaml.load(f)

fld = pth.dirname(cfgf)
rootname = cfg["rootname"]

t, rr, ri = np.loadtxt(pth.join(fld, "{}_rlin.txt".format(rootname)), unpack=True)
plt.figure()
plt.plot(t, rr)
plt.plot(t, ri)
plt.tight_layout()
plt.savefig("{}_rlin.png".format(rootname), dpi=300)
plt.close()

f, sr, si = np.loadtxt(pth.join(fld, "{}_slin.txt".format(rootname)), unpack=True)
plt.figure()
plt.plot(f, sr)
plt.plot(f, si)
plt.tight_layout()
plt.savefig("{}_slin.png".format(rootname), dpi=300)
