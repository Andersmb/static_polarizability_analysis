import functions as f
from pprint import pprint
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys, operator
import yaml

#data1 = f.get_mw_pol_fdu("mw_rawdata_0001_v2.csv", fieldstrength=0.001)
with open("hg_data.yaml") as f, open("hg_data_notextrapolated.yaml") as g:
    data1 = yaml.load(f)
    data2 = yaml.load(g)

# Only use molecules common in both data sets
skip = []
spin_filter = "ALL" # either "ALL", "NSP" or "SP"
molecules = []
for mol in data1.keys() + data2.keys():
    if mol in data1.keys() and mol in data2.keys() and mol not in skip and mol not in molecules:
        molecules.append(mol)

# Filter based on the spin information
if spin_filter != "ALL":
    molecules = filter(lambda mol: data1[mol]["spin"] == spin_filter, molecules)

# Define the xticks for the plots
xticks = range(len(molecules))

# Now extract the data we want: relative errors for the mean polarizability for each molecule
rel_err = [100 * (data2[mol]["ccsd(t)"]["mean"] / data1[mol]["mean"] - 1) for mol in molecules]

# Sort data based on the PBE relative error results
molecules_sorted, rel_err_sorted = zip(*sorted(zip(molecules, rel_err), reverse=True, key=operator.itemgetter(1)))

# Define edge colors based on spin polarizability
spin_colors = ["deepskyblue" if data1[mol]["spin"] == "NSP" else "crimson" for mol in molecules_sorted]

# Set up the figure with subplots
fontsize = 14
width=0.88

# Define custom lines that we will use for making a custom legend
# that explains the red and blue colors
lines = [Line2D([0], [0], color="deepskyblue"),
         Line2D([0], [0], color="crimson")]

fig = plt.figure(figsize=(15, 5))
ax = plt.gca()

ax.set_ylabel("Relative Error [%]", fontsize=fontsize)
#ax.set_yticks(range(9))
ax.tick_params("y", labelsize=fontsize)
ax.set_xlim(-1, len(molecules))

# Plot data
for i in range(len(molecules_sorted)):
    ax.bar(xticks[i], rel_err_sorted[i], color=spin_colors[i], edgecolor="black", width=width)

ax.grid(True, linestyle="--", linewidth=0.3)

# Place the molecule names on the xtick positions, rotation by 90 degrees
plt.xticks(xticks, [mol.upper() for mol in molecules_sorted], rotation=90, fontsize=10)

plt.legend(lines, ["Closed-shell", "Open-shell"], fontsize=fontsize)
plt.tight_layout()
plt.savefig("fig_ccsdt_extrap_vs_noextrap.png", dpi=300)
