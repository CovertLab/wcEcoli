import os
import csv
import numpy as np
import matplotlib.pyplot as plt


scale_by_rna = True

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ribosome_density_file = os.path.join(
	root_dir, "paper", "investigation", "SRR7759814_pass.fastq.gz_genelist.csv")
translation_efficiency_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "translationEfficiency.tsv")
rna_expression_file = os.path.join(
	root_dir, "paper", "investigation", "RNAseq_GeneWeiLi.tsv")

output_file = os.path.join(
	root_dir, "paper", "investigation", "analyze_green.pdf")


def open_file(filename, delimiter):
	with open(filename, "r") as f:
		reader = csv.reader(f, delimiter=delimiter)
		raw_data = [x for x in reader]
	return raw_data

data_green = np.array([x[0].split(",") for x in open_file(ribosome_density_file, ",")])
green_ribosome_density = data_green[1:, data_green[0].tolist().index("RPKM")].astype(float)
green_genes = data_green[1:, data_green[0].tolist().index("Alias")]

data_li = np.array(open_file(translation_efficiency_file, delimiter="\t"))
li_translation_efficiency = data_li[1:, data_li[0].tolist().index("translationEfficiency")]
li_genes = data_li[1:, data_li[0].tolist().index("name")]
li_rna_data = np.array(open_file(rna_expression_file, delimiter="\t"))
# todo: compute li's ribosome density

fig, ax = plt.subplots(1, 1, figsize=(11, 8.5))
for gene, data in zip(li_genes, li_translation_efficiency):
	if gene not in green_genes:
		continue

	x = -1 if data == "NA" else float(data)
	y = green_ribosome_density[green_genes.tolist().index(gene)]

	ax.scatter(x, y)

plt.savefig(output_file)
plt.close("all")