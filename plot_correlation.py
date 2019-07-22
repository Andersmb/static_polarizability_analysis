import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import os
from pprint import pprint

import functions

import yaml

with open("mw_data_0001v2.yaml") as f, \
     open("hg_data.yaml") as g, \
     open("data_gto_response.yaml") as h:
    mw = yaml.load(f)
    hg = yaml.load(g)
    orca = yaml.load(h)

skip = []
molecules = functions.common_species()
molecules = [mol for mol in molecules if mol not in skip]
print("Number of species: ", len(molecules))

re_hgmw = [100 * (hg[mol]["pbe"]["mean"] / mw[mol]["pbe"]["mean"] - 1) for mol in molecules]
re_hgorca = [100 * (hg[mol]["pbe"]["mean"] / orca[mol]["mean"] - 1) for mol in molecules]

FS = 14
fig = plt.Figure(figsize=(10, 10), dpi=100)
ax = plt.gca()
ax.set_ylabel("Relative Error [%]", fontsize=FS)
ax.set_xlabel("Relative Error [%]", fontsize=FS)
ax.tick_params("x", labelsize=FS)
ax.tick_params("y", labelsize=FS)

slope, intercept, r, pval, stderr = stats.linregress(re_hgmw, re_hgorca)
regress_x = np.arange(min(re_hgorca), max(re_hgorca), 0.0005)
regress_y = [slope * i + intercept for i in regress_x]

ax.scatter(re_hgmw, re_hgorca, edgecolors="black", s=2)
ax.plot(regress_x, regress_y, linewidth=0.5, linestyle="--", color="red")

print("Least-square fit parameter r^2: ", r**2)

plt.grid()
plt.tight_layout()
plt.savefig("fig_{}.png".format(__file__.split(".")[0]), dpi=100)
plt.show()
