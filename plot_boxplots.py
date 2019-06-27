import functions as f
from pprint import pprint
import matplotlib.pyplot as plt
import sys, operator

# Get data structures
data1 = f.get_HG_data("hg_data.csv", ["pbe", "spw92", "ccsd(t)"])
data2 = f.get_mw_pol_fdu("mw_rawdata_0001_v2.csv", fieldstrength=0.001)

# Only use molecules common in both data sets
skip = []
molecules = []
for mol in data1.keys() + data2.keys():
    if mol in data1.keys() and mol in data2.keys() and mol not in skip and mol not in molecules:
        molecules.append(mol)

rel_err_mw_cc = [100 * (data2[mol]["pbe"]["mean"] / data1[mol]["ccsd(t)"]["mean"] - 1) for mol in molecules]
rel_err_gto_cc = [100 * (data1[mol]["pbe"]["mean"] / data1[mol]["ccsd(t)"]["mean"] - 1) for mol in molecules]

fig = plt.figure(figsize=(10, 10))
ax = plt.gca()
FS = 16
props = {"linewidth": 3,}

# Plot data
ax.boxplot([rel_err_gto_cc, rel_err_mw_cc], 
            notch=False,
            showfliers=True,
            positions=[1, 1.2],
            showmeans=True,
            boxprops=props,
            medianprops={"linewidth": 4, "color": "gray"},
            flierprops={"marker": "o", "markersize": 10})

ax.set_xticklabels(["PBE/GTO\nvs\nCCSD(T)", "PBE/MW\nvs\nCCSD(T)"], fontsize=FS)
ax.set_ylabel("Relative Error [%]", fontsize=FS)
ax.tick_params(axis="y", labelsize=FS)
ax.set_xlim(0.9, 1.3)

plt.tight_layout()
plt.savefig("fig_boxplots", dpi=200)
