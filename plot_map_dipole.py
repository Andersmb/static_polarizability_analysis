"""
This data has been collected in the following way:

1) From the zero-field calculations, determine which of the dipole components that is very close to zero
   --> Found that both x and y directions were virtually zero
2) Perform calculations with applied electric fields in the x-direction of varying strengths:
   0.001, 0.002, 0.003, 0.004, 0.005, 0.006
3) Using the external field calculations, plot the x component of the dipole as a function of the applied field

Motivation for this: What is the shape of the function u(e)? This may explain the reason why certain molecules give
positive errors and others give negative errors.

"""


import matplotlib.pyplot as plt
import os, sys, yaml
from glob import glob
from pprint import pprint
sys.path.append("/Users/abr121/Documents/github/computational_chemistry")
from MRChem import MrchemOut
import numpy as np

datadir = os.path.join(os.getcwd(), "datafiles_map_u")

outputfiles = glob("{}/*.out".format(datadir))
molecules = map(lambda x: os.path.basename(x).split("_")[0], outputfiles)
fields = map(lambda x: os.path.basename(x).split("_")[2], outputfiles)

data = {mol:{"x": [], "y": []} for mol in molecules}
for output, mol, field in zip(outputfiles, molecules, fields):
    f = MrchemOut(output)
    if not f.normaltermination():
        continue
    else:
        data[mol]["y"].append(f.dipole_vector()[0])
        data[mol]["x"].append(float(field) / 1000)

# Now sort data based on increasing fields
for mol in data.keys():
    if len(data[mol]["x"]) == 0:
        continue
    data[mol]["x"], data[mol]["y"] = zip(*sorted(zip(data[mol]["x"], data[mol]["y"])))

FS = 16

fig = plt.Figure(figsize=(10,10))
ax = plt.gca()
ax.set_xlabel("Field strength along x [a.u.]", fontsize=FS)
ax.set_ylabel("x-component of dipole moment [a.u.]", fontsize=FS)
plt.grid()
for mol in data.keys():
    if len(data[mol]["y"]) > 1:
        # Extrapolate the line formed by the two first points: e=0 and e=0.001
        slope = (data[mol]["y"][1] - data[mol]["y"][0]) / (data[mol]["x"][1] - data[mol]["x"][0])
        xs = np.linspace(0, 0.006)
        ys = slope * xs + data[mol]["y"][0]

        # Then plot
        ax.scatter(data[mol]["x"], data[mol]["y"], marker="o", edgecolor="black", s=20, label=mol)
        ax.plot(xs, ys, linewidth=0.5, color="black", linestyle="--")
        ax.set_xlim(-0.0001, 0.0061)
        ax.set_ylim(-0.015, 1)

ax.legend()
fig.tight_layout()
plt.savefig("fig_{}.png".format(__file__.split(".")[0]), dpi=100)
plt.show()