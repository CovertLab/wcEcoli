#!/usr/bin/bash
#SBATCH --job-name=subgen_stage
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aj0204@stanford.edu
#SBATCH -p mcovert
#SBATCH -c 1
#SBATCH --mem=128GB
source /home/groups/mcovert/pyenv/versions/wcEcoli3/bin/activate

# The whole subgenerational pipeline, in one job script. Which stage runs is
# chosen by $SUBGEN_STAGE; submit_all.sh sets it and overrides --time/--mem/
# --job-name per stage on the sbatch command line (command-line options beat the
# #SBATCH directives above, which are only the safe worst-case defaults).
#
#   extract      the single heavy simOut pass; writes the per-cell matrices that
#                every classification analysis reads. No dependency.
#   cached       reads those matrices; cheap. Depends on `extract`.
#   independent  reads simOut directly and does NOT need `extract`, so it runs in
#                parallel with the other two: the CI reference implementation,
#                the five-definition comparison, cell health, and the
#                descriptive data-collection scripts.
#   cli          the standalone argparse CLIs. These are not analysis Plot
#                classes, so analysisCohort.py cannot run them and they must NOT
#                be listed in cohort/__init__.py -- the parallel runner evaluates
#                `mod.Plot` outside its try/except and would abort the whole run.
#                Depends on `cached`, which writes its per-gene table and its
#                binary expressed matrix.
#
# Run one stage by hand with:
#   SUBGEN_STAGE=cli sbatch subgen_stage.sbatch
# which is what replaced the old single-purpose rerun_*.sbatch file.

source "$(dirname "$0")/subgen_env.sh"

STAGE="${SUBGEN_STAGE:?set SUBGEN_STAGE to extract|cached|independent|cli}"

case "$STAGE" in

extract)
	# Reads each successful cell's simOut ONCE and writes the per-cell matrices
	# into each variant's plotOut under the fixed `subgen_extract` basename.
	run_analysis "$COH/subgen_extract.py" $VARIANTS "$OUT"
	;;

cached)
	# All of these read the extraction from each variant's plotOut, so they are
	# cheap. subgen_monomer_dynamics.py is the exception: it pulls time-series
	# traces from simOut, and uses the extraction only to pick its top-N subgen
	# gene set (skipping that half, with a message, if the cache is absent).
	for script in \
		subgen_expression_table \
		subgen_gen_expression \
		subgen_monomer_dynamics ; do
		run_analysis "$COH/$script.py" $VARIANTS "$OUT"
	done
	;;

independent)
	for script in \
		subgen_seed_ci \
		subgen_definitions \
		subgen_cell_health_metrics \
		subgen_transcription_prob \
		subgen_peak_counts \
		subgen_protein_distribution \
		subgen_subsample_mrna_protein_10k ; do
		run_analysis "$COH/$script.py" $VARIANTS "$OUT"
	done
	;;

cli)
	# These take explicit input paths, produced by the `extract` and `cached`
	# stages:
	#   <plotOut>/subgen_extract_synth_per_cell.tsv             (extract)
	#   <plotOut>/subgen_extract_frac_protein_zero_per_cell.tsv
	#   <plotOut>/subgen_extract_genes.tsv
	#   <plotOut>/subgen_expression_table.tsv                   (cached) -- canonical
	#   <plotOut>/subgen_gen_expression_expressed.tsv           (cached)
	#
	# `pergene` accepts either canonical per-gene table: the one above, or the
	# `independent` stage's subgen_seed_ci_pergene_successful.tsv. Both carry the
	# same def5_CI categories (verified gene-for-gene) and
	# sc.load_def5_categories() reads either spelling. --category defaults to
	# `subgen`, i.e. CI entirely below 1.
	SIM_DATA="$OUT/kb/simData_Modified.cPickle"
	if [ ! -f "$SIM_DATA" ]; then
		SIM_DATA="$OUT/kb/simData.cPickle"
	fi
	echo "  sim_data = $SIM_DATA"

	# One pass per variant directory. Cohort analyses write into
	# <variant>/plotOut, so each variant carries its own complete input set.
	shopt -s nullglob
	for variant_dir in "$OUT"/*_[0-9][0-9][0-9][0-9][0-9][0-9]; do
		PLOTOUT="$variant_dir/plotOut"
		[ -d "$PLOTOUT" ] || continue

		RAW="$PLOTOUT/subgen_extract"
		PERGENE="$PLOTOUT/subgen_expression_table.tsv"
		EXPRESSED="$PLOTOUT/subgen_gen_expression_expressed.tsv"
		CLI_OUT="$PLOTOUT/standalone"

		echo ""
		echo "########## $(basename "$variant_dir") ##########"
		if [ ! -f "$PERGENE" ]; then
			echo "SKIP: no $PERGENE -- run the 'cached' stage first."
			continue
		fi
		mkdir -p "$CLI_OUT"

		run_analysis "$COH/subgen_set_uniqueness.py" \
			"$EXPRESSED" "$PERGENE" "$CLI_OUT"

		run_analysis "$COH/subgen_set_uniqueness_publication_figures.py" \
			"$EXPRESSED" "$PERGENE" "$CLI_OUT"

		run_analysis "$COH/subgen_operon_membership.py" \
			"$SIM_DATA" "$PERGENE" "$CLI_OUT"

		run_analysis "$COH/subgen_operon_unit_jaccard.py" \
			"$EXPRESSED" "$PERGENE" "$SIM_DATA" "$CLI_OUT"

		run_analysis "$COH/subgen_protein_memory.py" \
			"${RAW}_frac_protein_zero_per_cell.tsv" \
			"${RAW}_synth_per_cell.tsv" \
			"${RAW}_genes.tsv" \
			"$PERGENE" "$CLI_OUT"

		# Chained: consumes subgen_protein_memory.py's own per-gene output.
		run_analysis "$COH/subgen_protein_memory_interactive.py" \
			"$CLI_OUT/subgen_protein_memory_pergene.tsv" "$SIM_DATA" "$CLI_OUT"
	done
	;;

*)
	echo "Unknown SUBGEN_STAGE '$STAGE'. Use: extract | cached | independent | cli" >&2
	exit 1
	;;
esac

subgen_summary
