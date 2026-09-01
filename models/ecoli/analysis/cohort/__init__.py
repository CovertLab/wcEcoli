# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"aa_conc.py",
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
	# Subgenerational-expression suite. THE ORDER IS LOAD-BEARING, and is only
	# honored when the runner is serial (analysisCohort.py -c 1, the default):
	#   subgen_extract.py writes the per-cell cache that subgen_expression_table
	#   and subgen_gen_expression HARD-require, and that subgen_definitions,
	#   subgen_seed_ci and subgen_monomer_dynamics degrade gracefully without.
	#   subgen_definitions.py writes the near-ubiquitous-absence table that
	#   subgen_cell_health_metrics.py optionally joins.
	# The six standalone argparse CLIs (subgen_operon_*, subgen_protein_memory*,
	# subgen_set_uniqueness*) and the subgen_helper_functions.py library have no
	# Plot class and must NOT be listed here: analysisBase.py evaluates mod.Plot
	# outside its try/except and would abort the whole run. Drive those from the
	# `cli` stage of sbatch/subgen/subgen_stage.sbatch instead.
	"subgen_extract.py",
	"subgen_seed_ci.py",
	"subgen_expression_table.py",
	"subgen_gen_expression.py",
	"subgen_definitions.py",
	"subgen_peak_counts.py",
	"subgen_protein_distribution.py",
	"subgen_transcription_prob.py",
	"subgen_monomer_dynamics.py",
	"subgen_subsample_mrna_protein_10k.py",
	"subgen_cell_health_metrics.py",
	# "transcriptFrequency.py",
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
	'SUBGEN': [         # subgenerational-expression suite, in dependency order
		"subgen_extract.py",
		"subgen_seed_ci.py",
		"subgen_expression_table.py",
		"subgen_gen_expression.py",
		"subgen_definitions.py",
		"subgen_peak_counts.py",
		"subgen_protein_distribution.py",
		"subgen_transcription_prob.py",
		"subgen_monomer_dynamics.py",
		"subgen_subsample_mrna_protein_10k.py",
		"subgen_cell_health_metrics.py",
		],
	'TRANSCRIPTION': [
		# "transcriptFrequency.py",
		# "transcriptionGenomeCoverage.py",
		# "transcriptionGenomeCoverageSecondHalf.py",
		],
	}
