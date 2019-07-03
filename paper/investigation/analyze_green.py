"""
Compare ribosome density: Mohammad vs Li
Compare translation efficiency: Mohammad scaled by Li or this study
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt


def open_file(filename, delimiter):
	with open(filename, "r") as f:
		reader = csv.reader(f, delimiter=delimiter)
		data = np.array([x for x in reader])
	return data


# Describe paths
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mohammad_data_file = os.path.join(
	root_dir, "paper", "investigation", "SRR7759814_pass.fastq.gz_genelist.csv")
li_data_file = os.path.join(
	root_dir, "paper", "investigation", "li_table_s4.tsv")
wcm_data_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "rna_seq_data", "rnaseq_seal_rpkm_mean.tsv")
wcm_genes_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "genes.tsv")
ribosome_genes_file = os.path.join(
	root_dir, "paper", "investigation", "ribosome_genes.tsv")
output_plot = os.path.join(
	root_dir, "paper", "investigation", "analyze_green.pdf")
output_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "translationEfficiency_alternate.tsv")
output2_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "translationEfficiency_alternate2.tsv")

# Identify genes of interest to highlight
rnap_genes = ["rpoA", "rpoB", "rpoC"]
ribosome_genes = [x[0] for x in open_file(ribosome_genes_file, "\n")]

# Mohammad et al. 2019 Ribosome Density (RPKM)
mohammad_data = np.array([x[0].split(",") for x in open_file(mohammad_data_file, ",")])
mohammad_data_header = mohammad_data[0].tolist()
mohammad_genes = mohammad_data[1:, mohammad_data_header.index("Alias")]
mohammad_ribosome_density = mohammad_data[1:, mohammad_data_header.index("RPKM")].astype(float)

# Li et al. 2014 Translation Efficiency (AU), mRNA (RPKM)
li_data = np.array(open_file(li_data_file, delimiter="\t"))
li_data_header = li_data[0].tolist()
li_genes = li_data[1:, li_data_header.index("Gene")]
li_mrna = np.array([x[1:-1] if "[" in x else x for x in li_data[1:, li_data_header.index("mRNA level (RPKM)")]], dtype=float)
li_translation_efficiency_raw = li_data[1:, li_data_header.index("Translation efficiency (AU)")]

# Exclude "NA" genes
mask = li_translation_efficiency_raw != "NA"
li_genes = li_genes[mask]
li_mrna = li_mrna[mask]
li_translation_efficiency = li_translation_efficiency_raw[mask].astype(float)

# Compute ribosome density
# Note: Translation Efficiency = Ribosome Density / mRNA
li_ribosome_density = li_translation_efficiency * li_mrna

# This study RNA (RPKM)
wcm_data = open_file(wcm_data_file, "\t")
wcm_data_header = wcm_data[0].tolist()
wcm_genes_egnumber = wcm_data[1:, wcm_data_header.index("Gene")]
wcm_rna = wcm_data[1:, wcm_data_header.index("M9 Glucose minus AAs")].astype(float)
wcm_genes_data = open_file(wcm_genes_file, "\t")
wcm_genes_data_header = wcm_genes_data[0].tolist()
egnumber_to_name = {x[wcm_genes_data_header.index("id")]: x[wcm_genes_data_header.index("symbol")] for x in wcm_genes_data[1:]}
name_to_egnumber = {x[wcm_genes_data_header.index("symbol")]: x[wcm_genes_data_header.index("id")] for x in wcm_genes_data[1:]}
wcm_genes = [egnumber_to_name.get(x) for x in wcm_genes_egnumber]


# Compute translation efficiency by scaling Mohammad et al. 2019 by this study or Li et al. 2014
def derive_translation_efficiency(ribosome_density_genes, ribosome_density, rna_genes, rna):
	derived_genes = []
	derived_translation_efficiency = []
	for gene, ribosome_density_value in zip(ribosome_density_genes, ribosome_density):
		if gene not in rna_genes:
			continue

		rna_value = rna[rna_genes.index(gene)]
		if rna_value == 0:
			continue
		derived_genes.append(gene)
		derived_translation_efficiency.append(ribosome_density_value / rna_value)
	derived_translation_efficiency = np.array(derived_translation_efficiency)
	return derived_genes, derived_translation_efficiency


derived_with_this_study_genes, derived_with_this_study_translation_efficiency = derive_translation_efficiency(
	mohammad_genes, mohammad_ribosome_density, wcm_genes, wcm_rna)
derived_with_li_genes, derived_with_li_translation_efficiency = derive_translation_efficiency(
	mohammad_genes, mohammad_ribosome_density, li_genes.tolist(), li_mrna)


# Analyze
def plot(ax, x_genes, x_data, y_genes, y_data):
	for x_gene, x_value in zip(x_genes, x_data):
		if x_gene not in y_genes:
			continue

		y_value = y_data[y_genes.index(x_gene)]

		color = "tab:orange" if x_gene in rnap_genes else "tab:green" if x_gene in ribosome_genes else "none"
		ax.scatter(x_value, y_value, c=color, edgecolors="tab:blue" if color is "none" else color)
	ax.plot([0, min(ax.get_xlim()[1], ax.get_ylim()[1])], [0, min(ax.get_xlim()[1], ax.get_ylim()[1])], "k")
	return


fig, axes_list = plt.subplots(3, 3, figsize=(11, 8.5))

ax = axes_list[0, 0]
plot(ax, li_genes, li_ribosome_density, mohammad_genes.tolist(), mohammad_ribosome_density)
ax.set_xlabel("Li ribosome density")
ax.set_ylabel("Mohammad ribosome density")

ax = axes_list[0, 1]
plot(ax, li_genes, li_translation_efficiency, derived_with_this_study_genes, derived_with_this_study_translation_efficiency)
ax.set_xlabel("Li translation efficiency")
ax.set_ylabel("Derived translation efficiency\n(Mohammad scaled by this study)")

# Replace missing data with average (mimics procedure in ParCa), Normalize
# n_genes = 4353
# data = derived_with_this_study_translation_efficiency.copy()
# to_plot = np.hstack((data, np.ones(n_genes - data.shape[0]) * data.mean()))
# to_plot /= to_plot.sum()
# axes_list[2].hist(to_plot, bins=np.arange(0, .03, .001))
# is_r_protein = [x in ribosome_genes for x in derived_with_this_study_genes]
# for value in to_plot[np.where(is_r_protein)[0]]:
# 	axes_list[2].axvline(value)


ax = axes_list[0, 2]
n, bins, patches = ax.hist(derived_with_this_study_translation_efficiency, bins=np.arange(0, 5.5, 0.25))
is_r_protein = [x in ribosome_genes for x in derived_with_this_study_genes]
ax.hist(derived_with_this_study_translation_efficiency[is_r_protein], bins=bins)

ax = axes_list[1, 0]
plot(ax, li_genes, li_mrna, wcm_genes, wcm_rna)
ax.set_xlabel("Li RNA (RPKM)")
ax.set_ylabel("This study RNA (RPKM)")

ax = axes_list[1, 1]
plot(ax, li_genes, li_translation_efficiency, derived_with_li_genes, derived_with_li_translation_efficiency)
ax.set_xlabel("Li translation efficiency")
ax.set_ylabel("Derived translation efficiency\n(Mohammad scaled by Li)")

# # Replace missing data with average (mimics procedure in ParCa), Normalize
# n_genes = 4353
# data = derived_with_li_translation_efficiency.copy()
# # id outliers
# mask = np.where(data > 100)[0]
# data[mask] = np.nan
# to_plot = np.hstack((data, np.ones(n_genes - data.shape[0]) * np.nanmean(data)))
# mask = np.isnan(to_plot)
# to_plot[mask] = np.nanmean(data)
#
# to_plot /= to_plot.sum()
# n, bins, patches = axes_list[5].hist(to_plot, bins=np.arange(0, 1.5e-3, 5e-5))
# is_r_protein = [x in ribosome_genes for x in derived_with_li_genes]
# axes_list[5].hist(to_plot[np.where(is_r_protein)[0]], bins=bins)

ax = axes_list[1, 2]
n, bins, patches = ax.hist(derived_with_li_translation_efficiency, bins=np.arange(0, 5.5, 0.25))
is_r_protein = [x in ribosome_genes for x in derived_with_li_genes]
ax.hist(derived_with_li_translation_efficiency[is_r_protein], bins=bins)

ax = axes_list[2, 2]
n, bins, patches = ax.hist(li_translation_efficiency, bins=np.arange(0, 5.5, 0.25))
is_r_protein = [x in ribosome_genes for x in li_genes]
ax.hist(li_translation_efficiency[is_r_protein], bins=bins)


ax = axes_list[2, 0]
original = li_translation_efficiency / li_translation_efficiency.sum()
# remove outliers
processed = derived_with_li_translation_efficiency.copy()
mask = processed > li_translation_efficiency.max()
processed[mask] = li_translation_efficiency.max()
new = processed / processed.sum()

max_val = 0
for gene in ribosome_genes:
	if gene not in li_genes or gene not in derived_with_li_genes:
		continue
	original_index = np.where(li_genes == gene)[0][0]
	new_index = derived_with_li_genes.index(gene)

	x = original[original_index]
	y = new[new_index]

	if max(x, y) > max_val:
		max_val = max(x, y)

	ax.scatter(x, y, color="none", edgecolors="tab:blue")
	ax.text(x, y, gene)
ax.set_xlim([0, max_val])
ax.set_ylim([0, max_val])
ax.plot([0, 1.25e-3], [0, max_val], "k")
ax.set_xlabel("Original TEs")
ax.set_ylabel("New derived TEs")


ax = axes_list[2, 1]
import cPickle
production_term = cPickle.load(open(os.path.join(root_dir, "paper", "investigation", "production_term.cpickle"), "rb"))
production_term_norm = production_term / production_term.sum()
protein_ids = cPickle.load(open(os.path.join(root_dir, "paper", "investigation", "protein_ids.cpickle"), "rb"))
egnumbers = [x.split("_")[0] for x in protein_ids]
names = [egnumber_to_name[x] for x in egnumbers]
norm_density = mohammad_ribosome_density / mohammad_ribosome_density.sum()

later_x = []
later_y = []
for name, protein_dist in zip(names, production_term_norm):
	if name not in mohammad_genes:
		continue
	y = norm_density[np.where(mohammad_genes == name)[0][0]]

	if name in ribosome_genes:
		later_x.append(protein_dist)
		later_y.append(y)
		continue

	ax.scatter(protein_dist, y, color="none", edgecolors="tab:blue")

ax.scatter(later_x, later_y, color="none", edgecolors="tab:orange")

max_val = max(ax.get_ylim()[1], ax.get_xlim()[1])
ax.set_xlim([0, max_val])
ax.set_ylim([0, max_val])
ax.plot([0, max_val], [0, max_val], "k")
ax.set_xlabel("Protein production")
ax.set_ylabel("Ribosome density")

plt.savefig(output_plot)
plt.close("all")


# Write alternate translation efficiency
def write_alternate(genes, data, output, cutoff):
	out = ['"geneId"\t"name"\t"translationEfficiency"']
	for gene, value in zip(genes, data):
		if value > cutoff:
			continue
		if value == 0:
			continue
		egnumber = name_to_egnumber.get(gene)
		if egnumber is None:
			continue
		# todo: look up synonyms

		out.append('"{}"\t"{}"\t{}'.format(egnumber, gene, value))

	with open(output, "w") as f:
		f.write("\n".join(out))


# write_alternate(derived_with_li_genes, derived_with_li_translation_efficiency, output_file, 10.)
write_alternate(derived_with_this_study_genes, derived_with_this_study_translation_efficiency, output2_file, 25.)
