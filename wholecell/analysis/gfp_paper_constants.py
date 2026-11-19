"""
Shared constants for the GFP burden paper's figure analyses.

Every paper analysis imports its generation window, grid definition and
new-gene identifiers from here, so that changing one of them is a single edit
rather than a search across forty-odd scripts. Nothing in this module reads
simulation output or draws anything; see gfp_paper_style for that.

Companion files:
	runscripts/gfp_paper/panels.tsv    panel label -> script, batch, variants
	runscripts/gfp_paper/batches.tsv   batch id -> path, git hash, sim_data hash
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

# --------------------------------------------------------------------------
# Simulation design
#
# Generations 0-7 are wildtype steady state, 8-15 are the induction transient
# while GFP accumulates, and 16-23 are the new GFP steady state. Averages for
# the paper use the last eight generations only.
NEW_GENE_INDUCTION_GEN = 8
N_GENERATIONS = 24
ANALYSIS_FIRST_GEN = 16
ANALYSIS_LAST_GEN = 23                    # inclusive
# Half-open range, matching the MIN_CELL_INDEX / MAX_CELL_INDEX convention of
# the existing heatmap analyses.
MIN_CELL_INDEX = ANALYSIS_FIRST_GEN
MAX_CELL_INDEX = ANALYSIS_LAST_GEN + 1

# Main grid (batch B1). Expression factors are log10 multipliers; the variant
# file raises 10 ** (factor - 1). The first entry must be 0, which is the
# no-GFP control. 8.33 and 8.66 split the one step over which transcription
# initiation goes from roughly 0.6 of the spacing cap to the cap itself.
NEW_GENE_EXPRESSION_FACTORS: List[float] = [0, 7, 8, 8.33, 8.66, 9, 10]
# Translation efficiency multipliers (RBS strength), descending as the variant
# file requires. 0 is transcribed but not translated.
NEW_GENE_TRANSLATION_EFFICIENCY_VALUES: List[float] = [10, 7, 5, 3, 2, 1, 0.1, 0]

GRID_SEEDS = 5                            # scan depth; a subset is deepened later
FOCUS_SEEDS = 60                          # batch B2, single-design panels

# gfp is inserted immediately after the last base of the lambda attB site. The
# resolved coordinate is recorded per batch in batches.tsv, because it shifts
# if the conflict search has to move it.
INSERTION_SITE_NAME = 'lambda attB'
INSERTION_POS = 807342
NEW_GENE_OPTION = 'egfp'


def n_variants(
		expression_factors: Sequence[float] = None,
		translation_efficiencies: Sequence[float] = None) -> int:
	"""
	Number of variant indices the grid produces.

	The variant file lays the grid out as a strict factorial with a single
	shared control, so this is 1 + (non-zero expression factors) x (RBS
	values), not the product of the two list lengths.
	"""
	factors = (NEW_GENE_EXPRESSION_FACTORS if expression_factors is None
		else expression_factors)
	efficiencies = (NEW_GENE_TRANSLATION_EFFICIENCY_VALUES
		if translation_efficiencies is None else translation_efficiencies)
	return 1 + (len(factors) - 1) * len(efficiencies)


def n_simulations(
		seeds: int = GRID_SEEDS, generations: int = N_GENERATIONS) -> int:
	"""Cell-generations in a full grid run, the unit Sherlock queues."""
	return n_variants() * seeds * generations


# --------------------------------------------------------------------------
# New gene identifiers
#
# Never hard-code 'NG-GFP-MONOMER'. The identifier changes with the new gene
# option (EGFP gives NG-EGFP-MONOMER), so look it up from sim_data instead.
def new_gene_ids(sim_data) -> Tuple[List[str], List[str], List[str]]:
	"""
	(cistron ids, mRNA/RNA ids, monomer ids) of the inserted gene(s).

	Raises:
		ValueError: if the simulation has no new genes, which usually means
			an analysis meant for an induction batch was pointed at a
			wildtype one.
	"""
	cistron_data = sim_data.process.transcription.cistron_data.struct_array
	monomer_data = sim_data.process.translation.monomer_data.struct_array
	cistron_ids = cistron_data[cistron_data['is_new_gene']]['id'].tolist()
	if not cistron_ids:
		raise ValueError(
			'No new genes in this sim_data. This analysis expects a new gene'
			' induction batch.')
	cistron_to_monomer = dict(
		zip(monomer_data['cistron_id'], monomer_data['id']))
	monomer_ids = [cistron_to_monomer[c] for c in cistron_ids
		if c in cistron_to_monomer]
	rna_data = sim_data.process.transcription.rna_data.struct_array
	rna_ids = [r for r in rna_data['id'].tolist()
		if any(r.startswith(c.split('_')[0]) for c in cistron_ids)]
	return cistron_ids, rna_ids, monomer_ids


# --------------------------------------------------------------------------
# Figure panels
#
# Panel labels are provisional: the locked storyboard may still be reordered,
# so scripts are named for their content and the label is attached at save
# time from runscripts/gfp_paper/panels.tsv. This mapping exists only so a
# script can assert that the label it was handed is one the paper knows about.
FIGURE_PANEL_COUNTS: Dict[int, int] = {1: 5, 2: 8, 3: 9, 4: 9, 5: 9, 6: 4}


def panel_labels() -> List[str]:
	"""Every panel label in the locked 2026.11.13 storyboard, e.g. '3F'."""
	return ['%d%s' % (figure, chr(ord('A') + i))
		for figure, count in sorted(FIGURE_PANEL_COUNTS.items())
		for i in range(count)]
