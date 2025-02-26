"""
Plot a scatterplot of the average monomer counts for each variant index.
"""

import os

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import mplcursors
import numpy as np
import pandas as pd


from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader

VAR_OF_INTEREST = "6"

# Load monomer counts data
filename = 'all_labeled_avg_monomer_counts.csv'
labeled_monomer_counts = pd.read_csv(filename)
monomer_ids = labeled_monomer_counts['monomer_ids'].values
is_essential = labeled_monomer_counts['is_essential'].values
monomer_counts = labeled_monomer_counts.drop(columns='monomer_ids')
monomer_counts = monomer_counts.drop(columns='is_essential')
col_labels = monomer_counts.columns
variants = [int(v) for v in col_labels]

normalized_monomer_counts = monomer_counts / monomer_counts.sum(axis=0)

# add 10e-10 to all values to avoid log(0) errors
normalized_monomer_counts = normalized_monomer_counts + 10e-11

poster_colors = {
    "light_gray": (0.75, 0.75, 0.75),
    "poster_green": (66/255, 170/255, 154/255),
    "poster_blue": (27/255, 132/255, 198/255),
    "poster_purple": (188/255, 140/255, 191/255),
    "poster_gold": (221/255, 203/255, 119/255),
    "poster_light_blue": (136/255, 205/255, 240/255),
    "poster_red": (202/255, 0/255, 32/255),
    }

POSTER_VARIANT_COLORS = [ (136/255, 205/255, 240/255), # light blue
						  (188/255, 140/255, 191/255), # light purple
						  (66/255, 170/255, 154/255), # light green
						  (221/255, 203/255, 119/255), # gold
						  (27/255, 132/255, 198/255)] # dark blue

colors = [poster_colors["light_gray"]] * len(normalized_monomer_counts["0"][:-1])
alphas = [0.5] * len(normalized_monomer_counts["0"][:-1])
markers = ['o' for i in range(len(normalized_monomer_counts["0"][:-1]))]
control_var = normalized_monomer_counts["3"][:-1] # TODO: Yingling, change this to show whatever variant you want
var_to_compare = normalized_monomer_counts[VAR_OF_INTEREST][:-1]
for p in range(len(control_var)):
    if np.log10(var_to_compare[p]) > np.log10(control_var[p]) + 0.25:
        if var_to_compare[p] > 10e-11 and control_var[p] > 10e-11:
            colors[p] = poster_colors['poster_blue']
            alphas[p] = 0.85
    elif np.log10(var_to_compare[p]) < np.log10(control_var[p]) - 0.25:
        if var_to_compare[p] > 10e-11 and control_var[p] > 10e-11:
            colors[p] = poster_colors['poster_gold']
            alphas[p] = 0.85

    # check if gene is essential
    if is_essential[p]:
        markers[p] = '^'

alphas = np.array(alphas)
colors = np.array(colors)
markers = np.array(markers)

# Create a scatterplot
fig, ax = plt.subplots(figsize=(7,7))
# force the plot to be square
plt.axis('equal')
# add the dashed line y = x to the plot
plt.plot([0, 1E-1], [0, 1E-1], color='black', linestyle='dashed')
# Loop over unique marker styles
unique_markers = set(markers)
for marker_style in unique_markers:
    mask = (markers == marker_style)
    plt.scatter(
        normalized_monomer_counts["0"][:-1][mask],
        normalized_monomer_counts[VAR_OF_INTEREST][:-1][mask],
        alpha=alphas[mask], s=35, linewidths=0, color=colors[mask], marker=marker_style
    )

sc = plt.scatter(normalized_monomer_counts["0"][:-1], normalized_monomer_counts[VAR_OF_INTEREST][:-1],
            alpha=1./510, s=35, linewidths=0, color=colors)

# TODO: CHANGE THESE TO SHOW VARIANT NUMBER, TRL EFF, EXP INDEX
plt.xlabel("Proteome Fraction: Wildtype (No GFP)", fontsize = 16)
plt.ylabel("Proteome Fraction: Higher GFP Production", fontsize = 16)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xscale('log')
plt.yscale('log')

plt.savefig("proteome_fraction_scatterplot_copy_" + VAR_OF_INTEREST + ".pdf", format="pdf")
# save as svg as well
plt.savefig("proteome_fraction_scatterplot_copy_" + VAR_OF_INTEREST + ".svg", format="svg")

# Create cursor for interactivity
cursor = mplcursors.cursor(sc, hover = True)

# Define what happens when a point is hovered over
@cursor.connect("add")
def on_add(sel):
    index = sel.index
    sel.annotation.set_text(monomer_ids[index])

plt.show()