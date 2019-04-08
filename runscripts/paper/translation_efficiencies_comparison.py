import os
import cPickle
import csv
import numpy as np
import matplotlib.pyplot as plt

from wholecell.utils.fitting import normalize

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
path_to_sim_data = os.path.join(root_dir, "out", "manual", "kb", "simData_Fit_1.cPickle")
path_to_derived_TEs = os.path.join(root_dir, "reconstruction", "ecoli", "flat", "derived_translation_efficiencies.tsv")
path_to_output = os.path.join(root_dir, "runscripts", "paper", "translation_efficiencies_comparison.png")

sim_data = cPickle.load(open(path_to_sim_data, "rb"))
original_TEs = normalize(sim_data.process.translation.translationEfficienciesByMonomer)
with open(path_to_derived_TEs, "r") as f:
	reader = csv.reader(f, delimiter="\t")
	derived_TEs = np.array([line for line in reader])

ids_r_proteins = sim_data.moleculeGroups.rProteins
ids_rnap = sim_data.moleculeGroups.rnapIds
ids_monomers = sim_data.process.translation.monomerData["id"]

fig, ax = plt.subplots(1, 1, figsize = (11, 8.5))
xticklabels = []
for i, molecule in enumerate(ids_r_proteins + ids_rnap):
	original_TE = original_TEs[np.where(ids_monomers == molecule)[0][0]]
	derived_TE = derived_TEs[np.where(derived_TEs[:, 2] == molecule)[0][0], 3].astype(float)
	fold_change = derived_TE / original_TE
	ax.bar(i, fold_change)
	xticklabels.append(derived_TEs[np.where(derived_TEs[:, 2] == molecule)[0][0], 1])
ax.set_xticks(range(len(xticklabels)))
ax.set_xticklabels(xticklabels, rotation=45, ha="right")
ax.set_ylabel("Required TE Adjustment")
plt.savefig(path_to_output)
plt.close("all")