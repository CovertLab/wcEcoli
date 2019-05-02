import os
import numpy as np
import matplotlib.pyplot as plt

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.simulation_data import SimulationDataEcoli

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__))))
PATH_TO_OUTPUT = os.path.join(root_dir, "paper", "investigation",
                              "{}.png".format(
                                  os.path.split(__file__)[-1].split(".")[0]))

# Generate pre-fit sim_data
raw_data = KnowledgeBaseEcoli()
sim_data = SimulationDataEcoli()
sim_data.initialize(
    raw_data=raw_data,
    basal_expression_condition="M9 Glucose minus AAs")

# Get IDs
r_protein_ids = [x.split("[c]")[0] for x in sim_data.moleculeGroups.rProteins]
rnap_monomer_ids = [x.split("[c]")[0] for x in sim_data.moleculeGroups.rnapIds]

# Map gene names
monomer_to_genes = {str(x['monomerId']): str(x['symbol'])
                    for x in raw_data.genes if x['type'] == 'mRNA'}

# Get translation efficiencies from Li et al. 2014, Table S4
genes = [str(x['name']) for x in raw_data.translationEfficiency]
genes_is_other = [x not in r_protein_ids or x not in rnap_monomer_ids
                  for x in genes]
TEs = np.array([float(x['translationEfficiency'])
                if x['translationEfficiency'] != 'NA' else -1.
                for x in raw_data.translationEfficiency])

# Plot
fig, axesList = plt.subplots(1, 2, figsize=(11, 8.5))

ax = axesList[0]
ax.scatter(np.zeros(len(genes_is_other)), TEs[genes_is_other], alpha=0.7)
ax.scatter(np.ones(len(r_protein_ids)),
           [TEs[genes.index(monomer_to_genes[x])] for x in r_protein_ids],
           alpha=0.7)
ax.scatter(2 * np.ones(len(rnap_monomer_ids)),
           [TEs[genes.index(monomer_to_genes[x])] for x in rnap_monomer_ids],
           alpha=0.7)
xticklabels = ["Other ({})".format(sum(genes_is_other)),
               "Ribosome ({})".format(len(r_protein_ids)),
               "RNA polymerase ({})".format(len(rnap_monomer_ids))]
ax.set_xticks(range(len(xticklabels)))
ax.set_xticklabels(xticklabels)
ax.set_ylabel("translation efficiency")
ax.set_xlabel("-1 indicates 'NA' in Li 2014")
ax.set_title("Translation efficiencies (Li 2014, Table S4)")

ax = axesList[1]
_, bins, _ = ax.hist(TEs[genes_is_other], orientation="horizontal", log=True,
                     alpha=0.2)
ax.hist([TEs[genes.index(monomer_to_genes[x])] for x in r_protein_ids],
        orientation="horizontal", log=True, alpha=0.6, bins=bins)
ax.hist([TEs[genes.index(monomer_to_genes[x])] for x in rnap_monomer_ids],
        orientation="horizontal", log=True, alpha=0.6, bins=bins)

plt.subplots_adjust(bottom=0.3, top=0.7)
plt.savefig(PATH_TO_OUTPUT)
plt.close("all")
