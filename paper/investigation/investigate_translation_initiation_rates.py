"""
Compares ribosome initiation rates in:

- Gorochowski et al. 2019
    - E. coli K12, harbors pBR322-derived plasmid containing inducible lacZ
    - MOPS minimal, 0.4% glycerol, 50 ug/mL arginine
    - Table S1, Sheet "LacZ", Column B "lacZ uninduced"
    - Ribosome initiation rate
    - Units: ribosome / second

- Li et al. 2014
    - E. coli K12 MG1655
    - MOPS minimal, 0.2% glucose
    - Table S1, Column C "MOPS complete"
    - Absolute protein synthesis rate
    - Units: molecules per generation
    - 1 cell cycle:
        21.5 minutes in MOPS complete
        26.5 minutes in MOPS methionine dropout
        56.3 minutes in MOPS minimal

"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.simulation_data import SimulationDataEcoli


root_dir = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__))))
this_dir = os.path.join(root_dir, "paper", "investigation")
path_to_gorochowski2019 = os.path.join(this_dir, "Gorochowski_2019.tsv")
path_to_li2014 = os.path.join(this_dir, "Li_2014.tsv")
path_to_output = os.path.join(this_dir, "{}.pdf".format(os.path.split(__file__)[-1].split(".")[0]))


def openfile(path_to_file):
    with open(path_to_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        data = np.array([x for x in reader])
    return data


# Load data
go2019_raw = openfile(path_to_gorochowski2019)
go2019 = []
for row in go2019_raw[1:]:
    go2019.append({
        "gene": row[0],
        "ribosome_per_sec": float(row[1]) if row[1] != "NA" else -1,
        "mrna_RPKM": float(row[2])})
go2019 = np.array(go2019)

li2014_raw = openfile(path_to_li2014)
li2014 = []
doubling_time = 56.3 # min
doubling_time__complete = 21.5 # min
for row in li2014_raw[1:]:
    ribosome_per_generation = float(row[2]) if "[" not in row[2] else float(row[2][1:-1])
    ribosome_per_generation__complete = float(row[1]) if "[" not in row[1] else float(row[1][1:-1])
    li2014.append({
        "gene": row[0],
        "ribosome_per_sec": ribosome_per_generation /doubling_time / 60.,
        "ribosome_per_sec__complete": ribosome_per_generation__complete / doubling_time__complete / 60.,
        "translation_efficiency": float(row[4]) if row[4] != "NA" else -1,
        "mrna_RPKM": float(row[3]) if "[" not in row[3] else float(row[3][1:-1])})
li2014 = np.array(li2014)

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
r_protein_genes = [monomer_to_genes[x] for x in r_protein_ids]
rnap_monomer_genes = [monomer_to_genes[x] for x in rnap_monomer_ids]

# Plot
fig, axes_list = plt.subplots(4, 3, figsize=(8.5, 11))
for axes_row, dataset, name in zip(
        axes_list[:2, :],
        [li2014, go2019],
        ["Li et al. 2014", "Gorochowski et al. 2019"]):

    ax = axes_row[0]
    genes = [x["gene"] for x in dataset]
    indexes_r_protein = [genes.index(x) for x in r_protein_genes if x in genes]
    indexes_rnap_monomers = [genes.index(x) for x in rnap_monomer_genes if x in genes]
    indexes_others = list(set(range(len(dataset))) - set([0]) - set(indexes_r_protein) - set(indexes_rnap_monomers))

    ax.scatter(np.zeros(len(indexes_others)), [x["ribosome_per_sec"] for x in dataset[indexes_others]])
    ax.scatter(np.ones(len(indexes_r_protein)), [x["ribosome_per_sec"] for x in dataset[indexes_r_protein]])
    ax.scatter(2 * np.ones(len(indexes_rnap_monomers)), [x["ribosome_per_sec"] for x in dataset[indexes_rnap_monomers]])

    xticklabels = ["Other\n({})".format(len(indexes_others)),
                   "Ribosome\n({})".format(len(indexes_r_protein)),
                   "RNA pol\n({})".format(len(indexes_rnap_monomers))]
    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels)
    ax.set_title("{}".format(name))
    ax.set_ylabel("Ribosome / s")


    ax = axes_row[1]
    _, bins, _ = ax.hist([x["ribosome_per_sec"] for x in dataset[indexes_others]],
                         orientation="horizontal", log=True, alpha=0.2)
    ax.hist([x["ribosome_per_sec"] for x in dataset[indexes_r_protein]],
            orientation="horizontal", log=True, alpha=0.6, bins=bins)
    ax.hist([x["ribosome_per_sec"] for x in dataset[indexes_rnap_monomers]],
            orientation="horizontal", log=True, alpha=0.6, bins=bins)


def plot_comparisons(axes_list, gene, highlight=False, text=False):
    facecolor = "r" if highlight else "none"
    edgecolor = "r" if highlight else "k"
    fontsize = 4

    if gene in ["rpmD", "rplX", "rpsH", "rplR", "rplN"]:
        facecolor = "b"

    ax = axes_list[0]
    if go2019[genes_go.index(gene)]["ribosome_per_sec"] != -1:
        ax.scatter(li2014[genes_li.index(gene)]["ribosome_per_sec"],
                   go2019[genes_go.index(gene)]["ribosome_per_sec"],
                   facecolors=facecolor,
                   edgecolors=edgecolor)
        if highlight and text:
            ax.text(li2014[genes_li.index(gene)]["ribosome_per_sec"],
                    go2019[genes_go.index(gene)]["ribosome_per_sec"],
                    gene,
                    fontsize=fontsize)

    ax = axes_list[1]
    if go2019[genes_go.index(gene)]["ribosome_per_sec"] != -1 and li2014[genes_li.index(gene)]["translation_efficiency"] != -1:
        ax.scatter(li2014[genes_li.index(gene)]["translation_efficiency"],
                   go2019[genes_go.index(gene)]["ribosome_per_sec"],
                   facecolors=facecolor,
                   edgecolors=edgecolor)
        if highlight and text:
            ax.text(li2014[genes_li.index(gene)]["translation_efficiency"],
                    go2019[genes_go.index(gene)]["ribosome_per_sec"],
                    gene,
                    fontsize=fontsize)

    ax = axes_list[2]
    test_zero_division = np.logical_and(li2014[genes_li.index(gene)]["mrna_RPKM"] != 0, go2019[genes_go.index(gene)]["mrna_RPKM"] != 0)
    test_data_exists = np.logical_and(li2014[genes_li.index(gene)]["ribosome_per_sec__complete"] != -1, go2019[genes_go.index(gene)]["ribosome_per_sec"] != -1)
    if test_zero_division and test_data_exists:
        ax.scatter(li2014[genes_li.index(gene)]["ribosome_per_sec__complete"] / li2014[genes_li.index(gene)]["mrna_RPKM"],
                   go2019[genes_go.index(gene)]["ribosome_per_sec"] / go2019[genes_go.index(gene)]["mrna_RPKM"],
                   facecolors=facecolor,
                   edgecolors=edgecolor)
        if highlight and text:
            ax.text(li2014[genes_li.index(gene)]["ribosome_per_sec__complete"] / li2014[genes_li.index(gene)]["mrna_RPKM"],
                    go2019[genes_go.index(gene)]["ribosome_per_sec"] / go2019[genes_go.index(gene)]["mrna_RPKM"],
                    gene,
                    fontsize=fontsize)

    ax = axes_list[3]
    if li2014[genes_li.index(gene)]["mrna_RPKM"] != 0 and go2019[genes_go.index(gene)]["ribosome_per_sec"] != -1:
        ax.scatter(li2014[genes_li.index(gene)]["ribosome_per_sec__complete"] / li2014[genes_li.index(gene)]["mrna_RPKM"],
                   go2019[genes_go.index(gene)]["ribosome_per_sec"],
                   facecolors=facecolor,
                   edgecolors=edgecolor)
        if highlight and text:
            ax.text(li2014[genes_li.index(gene)]["ribosome_per_sec__complete"] / li2014[genes_li.index(gene)]["mrna_RPKM"],
                    go2019[genes_go.index(gene)]["ribosome_per_sec"],
                    gene,
                    fontsize=fontsize)

    ax = axes_list[4]
    ax.scatter(li2014[genes_li.index(gene)]["mrna_RPKM"],
               go2019[genes_go.index(gene)]["mrna_RPKM"],
               facecolors=facecolor,
               edgecolors=edgecolor)
    if highlight and text:
        ax.text(li2014[genes_li.index(gene)]["mrna_RPKM"],
                go2019[genes_go.index(gene)]["mrna_RPKM"],
                gene,
                fontsize=fontsize)

    ax = axes_list[5]
    if li2014[genes_li.index(gene)]["mrna_RPKM"] != 0 and li2014[genes_li.index(gene)]["ribosome_per_sec__complete"] != -1:
        ax.scatter(li2014[genes_li.index(gene)]["translation_efficiency"],
                   li2014[genes_li.index(gene)]["ribosome_per_sec__complete"] / li2014[genes_li.index(gene)]["mrna_RPKM"],
                   facecolors=facecolor,
                   edgecolors=edgecolor)


    return

genes = list(set([x["gene"] for x in li2014]).intersection(set([x["gene"] for x in go2019])))
genes_li = [x["gene"] for x in li2014]
genes_go = [x["gene"] for x in go2019]
genes_indicate_on_top = []
axes_list_comparisons = [axes_list[2, 0], axes_list[2, 1], axes_list[3, 0], axes_list[3, 1], axes_list[2, 2], axes_list[3, 2]]
for gene in genes:
    if gene in r_protein_genes:
        genes_indicate_on_top.append(gene)
        continue

    plot_comparisons(axes_list_comparisons, gene, highlight=False, text=False)

for gene in genes_indicate_on_top:
    plot_comparisons(axes_list_comparisons, gene, highlight=True, text=True)


for ax, xlabel, ylabel, li_condition in zip(
        axes_list_comparisons[:-1],
        ["Ribosome initiation rate", "Translation efficiency", "Ribosome initiation rate / RPKM", "Ribosome initiation rate / RPKM", "RPKM"],
        ["Ribosome initiation rate", "Ribosome initiation rate", "Ribosome initiation rate / RPKM", "Ribosome initiation rate", "RPKM"],
        ["minimal", "complete", "complete", "complete", "complete"]):

    ax.set_xlabel("\nLi et al. 2014\n{}\nMOPS {} + glucose".format(xlabel, li_condition), fontsize=8)
    ax.set_ylabel("Gorochowski et al. 2019\n{}\nMOPS minimal + glycerol".format(ylabel), fontsize=8)
    ax.ticklabel_format(axis="both", style="sci", scilimits=(-1, 1))

ax = axes_list_comparisons[-1]
ax.set_xlabel("\nLi et al. 2014\n{}\nMOPS {} + glucose".format("Translation efficiency", "complete"), fontsize=8)
ax.set_ylabel("\nLi et al. 2014\n{}\nMOPS {} + glucose".format("Ribosome initiation rate / RPKM", "complete"), fontsize=8)
ax.ticklabel_format(axis="both", style="sci", scilimits=(-1, 1))


plt.subplots_adjust(wspace=0.8, hspace=0.8, top=0.95, bottom=0.1)
plt.savefig(path_to_output)
plt.close("all")
