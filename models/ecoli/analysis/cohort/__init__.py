# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"aa_conc.py",
	# Was "cell_health_metrics.py", which does not exist -- the runner swallows
	# ModuleNotFoundError, so this analysis silently never ran.
	"subgen_cell_health_metrics.py",
	"centralCarbonMetabolismCorrelationTimeCourse.py",
	"centralCarbonMetabolismScatter.py",
	"doubling_times_histogram_all.py",
	"ecocyc_table.py",
	"expression_dynamics.py",
	"growthDynamics.py",
	"growth_time_series.py",
	# "proteinFoldChangeVsTranscriptionFrequency.py",
	"histogramDoublingTime.py",
	"histogramFinalMass.py",
	"histogramGrowthRate.py",
	"initialVsFinalMass.py",
	"kinetics_flux_comparison.py",
	"mass_fraction_instantaneous_growth_rates.py",
	"promoter_probabilities.py",
	"proteinCopyNumberDistribution.py",
	"replisome_rnap_collision_per_gene.py",
	"rnaCopyNumberDistribution.py",
	# --- Subgenerational, tier def5_ci: these produce or consume the def5 + CI gene
	# classification. Run subgen_raw_extract.py FIRST; the def5 scripts read its
	# per-cell matrices. Order matters here.
	"subgen_raw_extract.py",
	"subgen_definition5_lineage_ci.py",
	"subgenerational_expression_table_def5.py",
	"transcriptFrequency_def5.py",
	"subgen_monomer_dynamics_def5.py",
	"subgen_per_generation_expression.py",
	# --- Subgenerational, tier descriptive: these classify no genes at all (fixed
	# curated panels, raw counts, subsampling). Just as current and maintained as
	# the tier above -- they simply make no subgen claim. Several previously ran
	# only because the sbatch files invoked them by path.
	"subgen_expression_definitions.py",
	"subgen_expression_definitions_complete.py",
	"subgen_peak_counts.py",
	"protein_distribution.py",
	"transcription_prob_subgen.py",
	"extract_subgen_monomer_counts.py",
	"molar_mass_monomers.py",
	"subgen_monomer_dynamics.py",
	"subsampl_mrna_protein_10k.py",
	"subsample_timepoints_for_scRNAseq_comp.py",
	"subsample_cell_monomer_mass.py",
	# --- Subgenerational, tier legacy -- see the SUBGEN_LEGACY tag. Kept on disk and
	# runnable on purpose, but not by default: they classify genes by a superseded
	# presence-frequency rule rather than def5 + CI.
	# "subgenerational_expression_table.py",
	# "id_subgen_monomers.py",
	# "transcriptFrequency.py",
	#
	# subgen_set_uniqueness.py and the other standalone subgen CLIs are NOT listed
	# here on purpose: they have no `Plot` class, and the parallel runner evaluates
	# `mod.Plot` outside its try/except, so listing one aborts the whole run.
	# See runscripts/sbatch/subgen/04_subgen_standalone_cli.sbatch.
	# "transcriptionGenomeCoverage.py",
	# "transcriptionGenomeCoverageSecondHalf.py",
	]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"proteinCopyNumberDistribution.py",  # TODO(jerry): an empty CORE list could be annoying, so include this?
		],
	'DIVISION': [
		"initialVsFinalMass.py",
		],
	'ECOCYC': [
		"ecocyc_table.py",
		],
	'GROWTH': [
		"aa_conc.py",
		"growth_time_series.py",
		],
	'HETEROGENEITY': [
		"proteinCopyNumberDistribution.py",
		"rnaCopyNumberDistribution.py",
		],
	'METABOLISM': [
		"aa_conc.py",
		"centralCarbonMetabolismCorrelationTimeCourse.py",
		"centralCarbonMetabolismScatter.py",
		"kinetics_flux_comparison.py",
		],
	# Definition 5 always means the CI form: a gene is `subgen` iff the 95% CI of
	# its per-lineage completed-transcript rate lies entirely below 1
	# transcript/generation (subgen_common.classify_def5_ci).
	'SUBGENERATIONAL': [   # run subgen_raw_extract.py FIRST; the rest read it
		"subgen_raw_extract.py",
		"subgen_definition5_lineage_ci.py",
		"subgenerational_expression_table_def5.py",
		"transcriptFrequency_def5.py",
		"subgen_monomer_dynamics_def5.py",
		"subgen_per_generation_expression.py",
		# Descriptive companions (no classification).
		"subgen_expression_definitions.py",
		"subgen_expression_definitions_complete.py",
		"subgen_cell_health_metrics.py",
		"subgen_peak_counts.py",
		"protein_distribution.py",
		"transcription_prob_subgen.py",
		"extract_subgen_monomer_counts.py",
		"molar_mass_monomers.py",
		"subgen_monomer_dynamics.py",
		"subsampl_mrna_protein_10k.py",
		"subsample_timepoints_for_scRNAseq_comp.py",
		"subsample_cell_monomer_mass.py",
		],
	# Superseded definitions, kept for provenance. These classify genes by
	# presence frequency (0 < p < 1), NOT by def5 + CI, so their gene lists are not
	# comparable to the def5_CI tables -- see each script's LEGACY header.
	'SUBGEN_LEGACY': [
		"subgenerational_expression_table.py",
		"id_subgen_monomers.py",
		"transcriptFrequency.py",
		],
	'PAPER': [
		"centralCarbonMetabolismScatter.py",
		"doubling_times_histogram_all.py",
		"expression_dynamics.py",
		"kinetics_flux_comparison.py",
		"histogramDoublingTime.py",
		"histogramFinalMass.py",
		"histogramGrowthRate.py",
		"mass_fraction_instantaneous_growth_rates.py",
		],
	'TRANSCRIPTION': [
		# "transcriptFrequency.py",
		# "transcriptionGenomeCoverage.py",
		# "transcriptionGenomeCoverageSecondHalf.py",
		],
	}
