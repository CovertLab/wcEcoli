import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import pdb

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from reconstruction.ecoli.dataclasses.process.transcription import RNA_SEQ_ANALYSIS

RNA_SEQ_ANALYSIS_RPKM = "seal_rpkm"
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
PATH_TO_COVERT_TPM = os.path.join(root_dir, "reconstruction", "ecoli", "flat", "rna_seq_data","rnaseq_{}_mean.tsv".format(RNA_SEQ_ANALYSIS))
PATH_TO_COVERT_RPKM = os.path.join(root_dir, "reconstruction", "ecoli", "flat", "rna_seq_data","rnaseq_{}_mean.tsv".format(RNA_SEQ_ANALYSIS_RPKM))
PATH_TO_LI = os.path.join(root_dir, "prototypes", "operons", "RNAseq_GeneWeiLi.tsv")
PATH_TO_OUTPUT = os.path.join(root_dir, "prototypes", "operons", "investigate_specific_toxins.{}")

#pdb.set_trace()

def openfile(filename):
    with open(filename, "r") as f:
        reader = csv.reader(f, delimiter = "\t")
        data = np.array([line for line in reader])
    return data

# Generate pre-fit sim_data
raw_data = KnowledgeBaseEcoli()
sim_data = SimulationDataEcoli()
sim_data.initialize(
    raw_data=raw_data,
    basal_expression_condition="M9 Glucose minus AAs")

# Get IDs
protein_ids = [
'EG10443-MONOMER',
'EG10442-MONOMER',
'EG12096-MONOMER',
'EG11250-MONOMER',
'G7917-MONOMER'
]

# get colormap
cmap = cm.get_cmap('rainbow',len(protein_ids))

# Map gene names
gene_to_monomer = {str(x['symbol']):str(x['monomerId']) for x in raw_data.genes if x['type'] == 'mRNA'}
gene_to_eg_number = {str(x['symbol']):str(x['id']) for x in raw_data.genes if x['type'] == 'mRNA'}

old_to_new_name = {
    "araH_1": "araH", # joined with araH_2
    "araH_2": "araH", # joined with araH_1
    "csdL": "tcdA",
    "ecoZ": "rbn",
    "gntU_1": "gntU", # joined with gntU_2
    "gntU_2": "gntU", # joined with gntU_1
    "rfaL": "waaL",
    "rfaS": "waaS",
    "rfaZ": "waaZ",
    "sdsP": "mdtP",
    "smf_1": "smf", # joined with smf_2
    "smf_2": "smf", # joined with smf_1
    "ybhO": "clsB",
    "ybjK": "rcdA",
    "ycfQ": "comR",
    "ychM": "dauA",
    "ycjZ": "pgrR",
    "ydcW": "patD",
    "ydiI": "menI",
    "yeiM": "psuT",
    "yhbJ": "rapZ",
    "yhhK": "panZ",
    "yibD": "waaH",
    "yigW_1": "tatD", # joined with yigW_2
    "yigW_2": "tatD", # joined with yigW_1
    "yijP": "eptC",
    "yjeP": "mscM",
    "yjgB": "ahr",
    "yjiE": "hypT",
    "yliL": None, # phantom gene
    "ymdC": "clsC",
}

# Convert in-house (this study) RNA-seq data to relative expression
# MG1655 strain, M9 + 0.4% glucose
# 4623 x 1, ID is EG number
raw_data = openfile(PATH_TO_COVERT_RPKM)
covert_genes = raw_data[1:, 0]
covert_data = raw_data[1:, np.where(raw_data[0] == "M9 Glucose minus AAs")[0][0]]

# Convert Gene-Wei Li 2014 (GSM1300282) RNA-seq data to counts
# MG1655 strain, MOPS media with 0.2% glucose, with full supplement (Neidhardt et al., 1974)
# 3724 x 1, ID is common name
raw_data = openfile(PATH_TO_LI)
li_genes = raw_data[:, 0]
li_data = raw_data[:, 1]

# Arrange data
data = []
for gene, li_expression in zip(li_genes, li_data):
    gene = old_to_new_name.get(gene, gene)
    if gene == None:
        continue
    if gene in ["glyU", "proK", "alaT"]: # Li 2014 measured a few tRNAs
        continue

    monomer = gene_to_monomer[gene]
    eg_number = gene_to_eg_number[gene]
    covert_expression = covert_data[np.where(covert_genes == eg_number)[0][0]]
    data.append([gene, monomer, covert_expression, li_expression])

with open(os.path.join(root_dir, 'prototypes', 'operons', 'covert_gwl_rna_seq.csv'),'w') as f:
    writer = csv.writer(f)
    writer.writerow(['Gene', 'Monomer', 'Covert_Expression', 'Li_Expression'])
    writer.writerows(data)
#import ipdb; ipdb.set_trace()

data = np.array(data)
prot_name_list = list(data[:,1])
protein_indexes = {}
for protein in protein_ids:
    protein_indexes[protein] = prot_name_list.index(protein)


pdb.set_trace()
# Plot
plt.figure(figsize=(5, 5))
ax = plt.subplot(1, 1, 1)

small_marker_size = 10
big_marker_size = 20
text_size = 10

# All data, genes of interested colored
x = np.log10(data[:, 2].astype(float))
y = np.log10(data[:, 3].astype(float))
ax.scatter(x, y, s=small_marker_size, facecolors="None", edgecolors="k", alpha=0.2)
for k, prot in enumerate(protein_ids):
    ax.scatter(x[protein_indexes[prot]], y[protein_indexes[prot]], s=big_marker_size, c=cmap(k))

ax.set_title("log10 mRNA expression")
ax.set_xlabel("data for wcEcoli\n(MG1655, M9 + glucose)")
ax.set_ylabel("Li 2014\n(MG1655, MOPS complete + glucose)")

# Format plot
ax_min = min([ax.get_xlim()[0], ax.get_ylim()[0]])
ax_max = max([ax.get_xlim()[1], ax.get_ylim()[1]])
ax.set_xlim([ax_min, ax_max])
ax.set_ylim([ax_min, ax_max])
ax.plot([ax_min, ax_max], [ax_min, ax_max], c="k", lw=0.5)

# Legend
patches = []
print(protein_ids)
for k, prot in enumerate(protein_ids):
    patches.append(mpatches.Patch(color=cmap(k), label=prot))

ax.legend(handles=patches, fontsize=8) #, green_patch, cyan_patch])


plt.tight_layout()
# Save
plt.savefig(PATH_TO_OUTPUT.format("pdf"))
# plt.savefig(PATH_TO_OUTPUT.format("png"))
plt.close("all")

