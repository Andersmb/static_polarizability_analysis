"""
Validation plot. Compare FD to LR to assess whether hyperpolarizabilities has contaminated FD.
"""
from pprint import pprint
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys, operator, yaml
from pprint import pprint
import math

import functions

# Get data structures
with open("mw_data_0001v2.yaml") as f, open("mw_data_response_lda.yaml") as g:
    data1 = yaml.load(f)
    data2 = yaml.load(g)

molecules = functions.incommon(data1.keys(), data2.keys())
skip = []
spin_filter = False

# Make sure none of the diagonal elements are exactly zero
# which is a consequence of the SCF not being converged
for mol in molecules:
    for c in data2[mol]["lda"]["diagonal"]:
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
rel_err = [100 * (data2[mol]["lda"]["mean"] / data1[mol]["lda"]["mean"] - 1) for mol in molecules]

# Sort data based on the PBE relative error results
molecules_sorted, rel_err_sorted = zip(*sorted(zip(molecules, rel_err), reverse=True, key=operator.itemgetter(1)))

print("Number of species: ", len(molecules_sorted))
print("Mean relative error (MRE): ", sum(rel_err) / len(rel_err))
print("RMSRE: ", math.sqrt(sum(map(lambda x: x**2, rel_err)) / len(molecules_sorted)))

# Define edge colors based on spin polarizability
spin_colors = ["deepskyblue" if data1[mol]["multiplicity"] == 1 else "crimson" for mol in molecules_sorted]

# Set up the figure with subplots
FS = 20
width = 0.7

# Define custom lines that we will use for making a custom legend
# that explains the red and blue colors
lines = [Line2D([0], [0], color="deepskyblue"),
         Line2D([0], [0], color="crimson")]

fig = plt.figure(figsize=(15, 5), dpi=100)
ax = plt.gca()
ax.tick_params(axis="y", labelsize=FS)
ax.set_ylabel("Relative Error [%]", fontsize=FS)

# Plot data
for i in range(len(molecules_sorted)):
    ax.bar(xticks[i], rel_err_sorted[i], color=spin_colors[i], edgecolor="black", width=width)

ax.grid(True, linestyle="--", linewidth=0.3)
ax.set_xlim(-1, len(molecules))

# Place the molecule names on the xtick positions, rotation by 90 degrees
plt.xticks(xticks, [mol.upper() for mol in molecules_sorted], rotation=90, fontsize=9)
ax.grid(True, linestyle="--", linewidth=0.3)

# Title
plt.title("LDA Finite Differences   vs   LDA Response", fontsize=FS+2)

fig.tight_layout()
plt.legend(lines, ["Closed-shell", "Open-shell"], fontsize=FS)
plt.savefig("fig_{}.png".format(__file__.split(".")[0]), dpi=700, transparent=True)
plt.show()
