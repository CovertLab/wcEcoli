"""
Plots RNA expression, translation efficiency, and degradation rate of ribosomal
proteins. Demonstrates computation that takes place in fit_sim_data_1.py.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from reconstruction.ecoli.fit_sim_data_1 import (
    setInitialRnaExpression, createBulkContainer, expressionConverge)
from wholecell.utils.fitting import normalize
from validation.ecoli.validation_data_raw import ValidationDataRawEcoli
from validation.ecoli.validation_data import ValidationDataEcoli


indicate_candidates = True
n_indicate = 5

root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
this_dir = os.path.dirname(os.path.realpath(__file__))
path_to_output = os.path.join(this_dir, "{}.png".format(
    os.path.split(__file__)[-1].split(".")[0]))


# Generate pre-fit sim_data
raw_data = KnowledgeBaseEcoli()
sim_data = SimulationDataEcoli()
sim_data.initialize(
    raw_data=raw_data,
    basal_expression_condition="M9 Glucose minus AAs")
transcription = sim_data.process.transcription
translation = sim_data.process.translation
translation_effs = translation.translationEfficienciesByMonomer


# Load data
doubling_time = sim_data.conditionToDoublingTime["basal"]
pre_fit_expression = setInitialRnaExpression(
    sim_data,
    transcription.rnaExpression["basal"],
    doubling_time)

pre_fit_bulk_container = createBulkContainer(
    sim_data,
    pre_fit_expression,
    doubling_time)

conc_dict = sim_data.process.metabolism.concentrationUpdates.\
    concentrationsBasedOnNutrients("minimal")
conc_dict.update(sim_data.mass.getBiomassAsConcentrations(doubling_time))

post_fit_expression, _, _, _, post_fit_bulk_container, _, _, _ = expressionConverge(
    sim_data,
    transcription.rnaExpression["basal"],
    conc_dict,
    doubling_time,
    disable_rnapoly_active_fraction_fitting=True,
    disable_ribosome_active_fraction_fitting=True)


# Load validation data
validation_data_raw = ValidationDataRawEcoli()
validation_data = ValidationDataEcoli()
validation_data.initialize(validation_data_raw, raw_data)
validation_protein_to_count = {x[0]: x[1]
                               for x in validation_data.protein.schmidt2015Data}


# Get IDs and indexes
ids_r_proteins = sim_data.moleculeGroups.rProteins
indexes_r_proteins_monomer = [np.where(translation.monomerData["id"] == x)[0][0]
    for x in ids_r_proteins]
indexes_r_proteins_rna = sim_data.relation.rnaIndexToMonomerMapping[
    indexes_r_proteins_monomer]


# Plot
def blank(ax):
    for key in ax.spines.keys():
        ax.spines[key].set_visible(False)
    ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
    return

def scatter(ax, x, y, set_y_lim=False):
    ax.scatter(x, y, facecolors="none", edgecolors="k")
    ax.set_xticks([0, len(x)])
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-1, 1))
    ax.set_xlabel("R-proteins")

    if set_y_lim:
        ax.set_ylim([0, max(y) / 0.9])
    return

pre_fit_rna_times_te = (
        normalize(pre_fit_expression)[indexes_r_proteins_rna]\
        * normalize(translation_effs)
        [indexes_r_proteins_monomer])
ranked_r_proteins = np.argsort(pre_fit_rna_times_te)

fig, axes_list = plt.subplots(3, 6, figsize=(11, 8.5))
title_size = 10

for ax, label in zip(
        axes_list[:, 0],
        ["Pre-Fit", "Post-Fit", "Fold Change\n(Post / Pre)"]):
    blank(ax)
    ax.text(0.2, 0.5, label,
            fontsize=12,
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes)


for i, (expression, bulk_container) in enumerate(zip(
        [pre_fit_expression, post_fit_expression],
        [pre_fit_bulk_container, post_fit_bulk_container])):

    ax = axes_list[i, 1]
    scatter(ax, range(len(ids_r_proteins)), expression[indexes_r_proteins_rna],
            set_y_lim=True)
    ax.set_ylabel("Relative Expression")
    if i == 0:
        ax.set_title("RNA\n(this study)\n", fontsize=title_size)
    if indicate_candidates:
        for j in ranked_r_proteins[:n_indicate]:
            ax.scatter(j, expression[indexes_r_proteins_rna[j]])


    ax = axes_list[i, 2]
    scatter(ax, range(len(ids_r_proteins)),
            translation_effs[indexes_r_proteins_monomer])
    ax.set_ylabel("AU")
    if i == 0:
        ax.set_title("Translation Efficiency\n(Li et al. 2014)\n",
                     fontsize=title_size)
    if indicate_candidates:
        for j in ranked_r_proteins[:n_indicate]:
            ax.scatter(j, translation_effs[indexes_r_proteins_monomer[j]])


    ax = axes_list[i, 3]
    y = normalize(expression)[indexes_r_proteins_rna]\
        * normalize(translation_effs)[indexes_r_proteins_monomer]
    scatter(ax, range(len(ids_r_proteins)), y, set_y_lim=True)
    ax.set_ylabel("AU")
    if i == 0:
        ax.set_title("norm(RNA) * norm(TE)\n", fontsize=title_size)
    if indicate_candidates:
        for j in ranked_r_proteins[:n_indicate]:
            ax.scatter(j, y[j])


    ax = axes_list[i, 4]
    scatter(ax, range(len(ids_r_proteins)), bulk_container.counts(ids_r_proteins))
    ax.set_ylabel("Counts")
    if i == 0:
        ax.set_title("Protein\n(bulk container)\n", fontsize=title_size)
    if indicate_candidates:
        for j in ranked_r_proteins[:n_indicate]:
            ax.scatter(j, bulk_container.count(ids_r_proteins[j]))


    ax = axes_list[i, 5]
    scatter(ax, range(len(ids_r_proteins)),
            [validation_protein_to_count[x] for x in ids_r_proteins])
    ax.set_ylabel("Counts")
    if i == 0:
        ax.set_title("Protein\n(Schmidt et al. 2016)\n", fontsize=title_size)
    if indicate_candidates:
        for j in ranked_r_proteins[:n_indicate]:
            ax.scatter(j, validation_protein_to_count[ids_r_proteins[j]])


ax = axes_list[2, 1]
scatter(ax, range(len(ids_r_proteins)),
        (post_fit_expression / pre_fit_expression)[indexes_r_proteins_rna])
ax.set_ylabel("AU")
if indicate_candidates:
    for j in ranked_r_proteins[:n_indicate]:
        ax.scatter(
            j,
            post_fit_expression[indexes_r_proteins_rna[j]]
            / pre_fit_expression[indexes_r_proteins_rna[j]])

for col in range(axes_list.shape[1]):
    ax0, ax1, _ = axes_list[:, col]
    ax_min = min(ax0.get_ylim()[0], ax1.get_ylim()[0])
    ax_max = max(ax0.get_ylim()[1], ax1.get_ylim()[1])
    ax0.set_ylim([ax_min, ax_max])
    ax1.set_ylim([ax_min, ax_max])

for ax in axes_list[2, 2:]:
    blank(ax)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
handles = [mpatches.Patch(color=colors[x % len(colors)]) for x in range(n_indicate)]
labels = ["R-protein with lowest Pre-Fit (RNA * TE)"]
for i in range(n_indicate - 1):
    labels.append("next lowest")

plt.subplots_adjust(wspace=0.7, hspace=0.6, left=0.05, right=0.97, bottom=0.1, top=0.9)
plt.legend(handles, labels)
plt.savefig(path_to_output)
plt.close("all")
