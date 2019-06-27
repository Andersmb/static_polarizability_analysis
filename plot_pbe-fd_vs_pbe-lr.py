"""
Validation plot. Compare FD to LR to assess whether hyperpolarizabilities has contaminated FD.
"""

import functions
from pprint import pprint
import matplotlib.pyplot as plt
import sys, operator, yaml
from pprint import pprint

# Get data structures
with open("mw_data_0001v2.yaml") as f, \
        open("mw_data_response_lda.yaml") as g, \
        open("mw_data_response_pbe.yaml") as h, \
        open("datafiles_orca_response_V2.yaml") as i:
    data1 = yaml.load(f)
    data2 = yaml.load(g)
    data3 = yaml.load(h)
    data4 = yaml.load(i)

molecules = functions.incommon(data1.keys(), data2.keys(), data4.keys())
skip = []
spin_filter = False

# Make sure none of the diagonal elements are exactly zero
# which is a consequence of the SCF not being converged
for mol in molecules:
    for c in data2[mol]["pbe"]["diagonal"]:
        if float(c) == 0.0:
            skip.append(mol)

# Filter out molecules matching the spin polarization
if spin_filter:
    molecules = filter(lambda mol: data1[mol]["multiplicity"] > 1, molecules)

# Skip desired molecules
molecules = [mol for mol in molecules if mol not in skip]

# Define the xticks for the plots
xticks = range(len(molecules))

# Now extract the data we want: relative errors for the mean polarizability for each molecule
rel_err = [100 * (data2[mol]["pbe"]["mean"] / data1[mol]["pbe"]["mean"] - 1) for mol in molecules]

# Sort data based on the PBE relative error results
molecules_sorted, rel_err_sorted = zip(*sorted(zip(molecules, rel_err), reverse=True, key=operator.itemgetter(1)))

print("Number of species: ", len(molecules_sorted))
print("Mean relative error (MRE): ", sum(rel_err_mw_gto) / len(rel_err_mw_gto))
print("RMSRE: ", math.sqrt(sum(map(lambda x: x**2, rel_err_mw_gto)) / len(molecules_sorted)))

# Define edge colors based on spin polarizability
spin_colors = ["deepskyblue" if data1[mol]["multiplicity"] == 1 else "crimson" for mol in molecules_sorted]

# Set up the figure with subplots
FS = 14
width = 0.7

fig = plt.figure(figsize=(15, 5), dpi=100)
ax = plt.gca()
ax.tick_params(axis="y", labelsize=FS)
ax.set_ylabel("Relative Error [%]", fontsize=FS)

# Plot data
for i in range(len(molecules_sorted)):
    ax.bar(xticks[i], rel_err_sorted[i], color=spin_colors[i], edgecolor="black", width=width)

ax.grid(True, linestyle="--")
ax.set_xlim(-1, len(molecules))

# Place the molecule names on the xtick positions, rotation by 90 degrees
plt.xticks(xticks, [mol.upper() for mol in molecules_sorted], rotation=90, fontsize=FS)

fig.tight_layout()
plt.savefig("fig_{}.png".format(__file__.split(".")[0]), dpi=100)
plt.show()
