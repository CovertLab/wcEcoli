#!/usr/bin/bash
#SBATCH --job-name=subgen_launch
#SBATCH --time=00:10:00
#SBATCH -p owners
#SBATCH -c 1
#SBATCH --mem=2GB
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=aj0204@stanford.edu
#SBATCH -o /scratch/users/aj0204/subgen_logs/%x-%j.out


# The whole subgenerational analysis pipeline, in one file.
# 	run by submitting the following in to the command line:
#       sbatch subgen.sbatch <sim_dir> [burn_in] [stage]
#       bash   subgen.sbatch <sim_dir> [burn_in] [stage]

# submits a job that then launches the different stages. Runs through all variants in the sim output

#   sim_dir   sim output directory. 

#   burn_in   generations to skip (default 8). 

#   stage     Subgen scripts are divided into four stages to help with runtime : 
#   			extract | cached | independent | cli. 
# 			Specifying the stage runs only that in case you need to rerun something.
# 			Note if a stage depends on the outputs of another it may fail if the other has not run yet 

# Not specifying stage runs the following default:

#   1. extract: (heavy simOut pass; writes the per-cell matrices)
#   2. cached: requires extract, generates subgen_exp[ression, monomer_dynamics, and gen_expression tables and figures
#   3. cli: requires cached, one job per variant standalone CLIs
#   4. independent: (no dependency; reads simOut directly, runs alongside)

# Optional environment:
#   SUBGEN_CLI_PARTITION / SUBGEN_CLI_TIME  cli placement (default owners / 24:00:00).
#       owners is preemptible and a requeue restarts set_uniqueness from zero, so
#       for a long set: SUBGEN_CLI_PARTITION=mcovert SUBGEN_CLI_TIME=2-00:00:00
#   SUBGEN_REPO     repo to run from (default $HOME/repos/wcEcoli)
#   SUBGEN_DRY_RUN=1  print the sbatch / python commands instead of running them

# All output goes to <variant>/plotOut/

set -euo pipefail

REPO="${SUBGEN_REPO:-$HOME/repos/wcEcoli}"
SELF="$REPO/sbatch/subgen/subgen.sbatch"
COH="$REPO/models/ecoli/analysis/cohort"
LOG_DIR="$SCRATCH/subgen_logs"
DRY_RUN="${SUBGEN_DRY_RUN:-0}"

usage() {
	echo "Usage: sbatch|bash subgen.sbatch <sim_dir> [burn_in] [extract|cached|independent|cli]" >&2
	exit 1
}

[ $# -ge 1 ] || usage
[ -d "$1" ] || { echo "Not a directory: $1" >&2; usage; }
OUT="$(cd "$1" && pwd)"
BURN_IN="${2:-8}"
[[ "$BURN_IN" =~ ^[0-9]+$ ]] || { echo "burn_in must be an integer, got '$BURN_IN'" >&2; usage; }
ONLY_STAGE="${3:-}"
case "$ONLY_STAGE" in
	''|extract|cached|independent|cli) ;;
	*) echo "Unknown stage '$ONLY_STAGE'" >&2; usage ;;
esac

export SUBGEN_IGNORE_FIRST_N_GENS="$BURN_IN"

# Every variant directory, as bare indices (condition_000002 -> 2).
VARIANT_IDX=""
for d in "$OUT"/*_[0-9][0-9][0-9][0-9][0-9][0-9]; do
	[ -d "$d" ] || continue
	n="${d##*_}"
	VARIANT_IDX="$VARIANT_IDX $((10#$n))"
done
VARIANT_IDX="${VARIANT_IDX# }"
[ -n "$VARIANT_IDX" ] || { echo "No *_NNNNNN variant directories in $OUT" >&2; exit 1; }


# job launcher
if [ -z "${_SUBGEN_STAGE_JOB:-}" ]; then
	SIM_NAME="$(basename "$OUT")"
	SIM_NAME="${SIM_NAME#*__}"  # drop the timestamp prefix
	TAG="sg_${SIM_NAME}_i${BURN_IN}"
	CLI_PARTITION="${SUBGEN_CLI_PARTITION:-owners}"
	CLI_TIME="${SUBGEN_CLI_TIME:-24:00:00}"
	mkdir -p "$LOG_DIR"

	echo "Submitting subgen pipeline:"
	echo "  OUT      = $OUT"
	echo "  variants = $VARIANT_IDX"
	echo "  burn_in  = $BURN_IN"
	echo "  stage    = ${ONLY_STAGE:-<all>}"
	echo "  logs     = $LOG_DIR/${TAG}_<stage>-<jobid>.out"

	submit() {
		local suffix="$1" stage="$2" variant="$3" time="$4" mem="$5"; shift 5
		local cmd=(sbatch --parsable --export=ALL --requeue --mail-type=ALL
			--job-name="${TAG}_${suffix}" --time="$time" --mem="$mem"
			-o "$LOG_DIR/%x-%j.out" "$@" "$SELF" "$OUT" "$BURN_IN")
		if [ "$DRY_RUN" = 1 ]; then
			echo "[dry-run] _SUBGEN_STAGE_JOB=$stage SUBGEN_CLI_VARIANT=$variant ${cmd[*]}" >&2
			echo "DRY_$suffix"
		else
			env _SUBGEN_STAGE_JOB="$stage" SUBGEN_CLI_VARIANT="$variant" "${cmd[@]}"
		fi
	}

	submit_cli() {  # submit_cli [extra sbatch args...]
		local idx jid
		for idx in $VARIANT_IDX; do
			jid=$(submit "cli_v$idx" cli "$idx" "$CLI_TIME" 8GB -p "$CLI_PARTITION" "$@")
			echo "cli v$idx     -> job $jid (-p $CLI_PARTITION, $CLI_TIME${1:+, $*})"
		done
	}

	case "$ONLY_STAGE" in
	extract)     echo "extract     -> job $(submit extract extract - 12:00:00 128GB -p owners)" ;;
	independent) echo "independent -> job $(submit independent independent - 12:00:00 128GB -p owners)" ;;
	cached)      echo "cached      -> job $(submit cached cached - 4:00:00 64GB -p owners)" ;;
	cli)         submit_cli ;;
	'')
		jid_extract=$(submit extract extract - 12:00:00 128GB -p owners)
		echo "extract     -> job $jid_extract"
		jid_independent=$(submit independent independent - 12:00:00 128GB -p owners)
		echo "independent -> job $jid_independent"
		jid_cached=$(submit cached cached - 4:00:00 64GB -p owners \
			--dependency=afterok:"$jid_extract")
		echo "cached      -> job $jid_cached (afterok:$jid_extract)"
		submit_cli --dependency=afterany:"$jid_cached"
		;;
	esac
	echo "Check with: squeue --me"
	exit 0
fi

# stages 
STAGE="$_SUBGEN_STAGE_JOB"

set +u  # the venv's activate script references unset variables
source /home/groups/mcovert/pyenv/versions/wcEcoli3/bin/activate
set -u

cd "$REPO"
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1

GIT_DIRTY="$(git -C "$REPO" status --porcelain)"
echo "=== subgen pipeline configuration ==="
echo "  REPO     = $REPO"
echo "  OUT      = $OUT"
echo "  VARIANTS = $VARIANT_IDX"
echo "  STAGE    = $STAGE"
echo "  BURN_IN  = $BURN_IN"
echo "  git      = $(git -C "$REPO" rev-parse --short HEAD)" \
	"$(git -C "$REPO" rev-parse --abbrev-ref HEAD)" \
	"$([ -n "$GIT_DIRTY" ] && echo '(DIRTY)' || echo '(clean)')"
echo "====================================="
if [ -n "$GIT_DIRTY" ]; then
	echo "WARNING: the working tree is DIRTY. Outputs will be tagged git_dirty=true"
	echo "         and will not be reproducible. Commit before a publication run."
fi

# Run one analysis script, recording its exit status without aborting the job:
# the scripts are independent, so one failure should not cost the rest.
# subgen_summary at the end makes failures impossible to miss.
SUBGEN_FAILED=()
LAST_STATUS=0
run_analysis() {
	local script="$1"; shift
	LAST_STATUS=0
	echo "=== Running $(basename "$script") $* ==="
	if [ "$DRY_RUN" = 1 ]; then
		echo "[dry-run] python $script $*"
		return 0
	fi
	# Capture the status directly: inside `if python ...; then`, `$?` is the
	# status of the `if` construct, not of python.
	python "$script" "$@" || LAST_STATUS=$?
	if [ "$LAST_STATUS" -eq 0 ]; then
		echo "--- OK: $(basename "$script")"
	else
		echo "--- FAILED (exit $LAST_STATUS): $(basename "$script")"
		SUBGEN_FAILED+=("$(basename "$script")")
	fi
}

# Every cohort script, once per variant.
run_cohort() {  # run_cohort <script-name>
	local idx
	for idx in $VARIANT_IDX; do
		run_analysis "$COH/$1.py" -v "$idx" "$OUT"
	done
}

run_unless_done() {  # run_unless_done <sentinel> <script> [args...]
	local sentinel="$1"; shift
	if [ -f "$sentinel" ] && [ -f "$sentinel.burn_in" ] \
			&& [ "$(cat "$sentinel.burn_in")" = "$BURN_IN" ]; then
		echo "--- SKIP (already complete for burn_in=$BURN_IN): $(basename "$sentinel")"
		return 0
	fi
	rm -f "$sentinel.burn_in"
	run_analysis "$@"
	if [ "$DRY_RUN" != 1 ] && [ "$LAST_STATUS" -eq 0 ] && [ -f "$sentinel" ]; then
		echo "$BURN_IN" > "$sentinel.burn_in"
	fi
}

subgen_summary() {
	if [ ${#SUBGEN_FAILED[@]} -eq 0 ]; then
		echo "=== All scripts completed successfully. ==="
	else
		echo "=== ${#SUBGEN_FAILED[@]} script(s) FAILED: ${SUBGEN_FAILED[*]} ==="
		return 1
	fi
}

case "$STAGE" in

extract)
	# Reads each successful cell's simOut ONCE and writes the per-cell matrices into each variant's plotOut 
	run_cohort subgen_extract
	;;

cached)
	# All of these read the extraction back out of plotOut
	# subgen_monomer_dynamics.py is the exception: it pulls time-series traces from
	# simOut, and uses the extraction only to pick its top-N subgen gene set
	
	for script in \
		subgen_expression_table \
		subgen_gen_expression \
		subgen_monomer_dynamics ; do
		run_cohort "$script"
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
		run_cohort "$script"
	done
	;;

cli)
	# The standalone argparse CLIs. They are not analysis Plot classes, so
	# analysisCohort.py cannot run them (and they must NOT be listed in
	# cohort/__init__.py), and their input paths are built by hand. <p> below is
	# <variant>/plotOut:
	#   <p>/subgen_extract_synth_per_cell.tsv             (extract)
	#   <p>/subgen_extract_frac_protein_zero_per_cell.tsv
	#   <p>/subgen_extract_genes.tsv
	#   <p>/subgen_expression_table.tsv                   (cached) -- canonical
	#   <p>/subgen_gen_expression_expressed.tsv           (cached)
	#
	# `pergene` accepts either canonical per-gene table: the one above, or the
	# `independent` stage's subgen_seed_ci_pergene_successful.tsv. Both carry the
	# same def5_CI categories and sc.load_def5_categories() reads either spelling.
	SIM_DATA="$OUT/kb/simData_Modified.cPickle"
	if [ ! -f "$SIM_DATA" ]; then
		SIM_DATA="$OUT/kb/simData.cPickle"
	fi
	echo "  sim_data = $SIM_DATA"

	# The launcher sets SUBGEN_CLI_VARIANT to one index per job; "-" or unset means
	# every variant (e.g. running this mode by hand).
	CLI_VARIANTS="${SUBGEN_CLI_VARIANT:--}"
	[ "$CLI_VARIANTS" != "-" ] || CLI_VARIANTS="$VARIANT_IDX"
	echo "  cli variants = $CLI_VARIANTS"

	shopt -s nullglob
	for idx in $CLI_VARIANTS; do
		variant_glob=("$OUT"/*_$(printf '%06d' "$idx"))
		if [ ${#variant_glob[@]} -eq 0 ]; then
			echo "SKIP: no variant directory with index $idx in $OUT"
			continue
		fi
		variant_dir="${variant_glob[0]}"

		PLOTOUT="$variant_dir/plotOut"
		RAW="$PLOTOUT/subgen_extract"
		PERGENE="$PLOTOUT/subgen_expression_table.tsv"
		EXPRESSED="$PLOTOUT/subgen_gen_expression_expressed.tsv"
		CLI_OUT="$PLOTOUT/standalone"

		echo ""
		echo "########## $(basename "$variant_dir") ##########"
		if [ ! -f "$PERGENE" ] && [ "$DRY_RUN" != 1 ]; then
			echo "SKIP: no $PERGENE -- run the 'cached' stage first."
			continue
		fi
		[ "$DRY_RUN" = 1 ] || mkdir -p "$CLI_OUT"

		# The expensive one: hours, scaling with cell count. Its sentinel is the
		# summary JSON, which run() writes after every plot.
		run_unless_done "$CLI_OUT/subgen_set_uniqueness_summary.json" \
			"$COH/subgen_set_uniqueness.py" "$EXPRESSED" "$PERGENE" "$CLI_OUT"

		run_analysis "$COH/subgen_set_uniqueness_publication_figures.py" \
			"$EXPRESSED" "$PERGENE" "$CLI_OUT"

		run_analysis "$COH/subgen_operon_membership.py" \
			"$SIM_DATA" "$PERGENE" "$CLI_OUT"

		# The second-most expensive (tens of minutes). Its sentinel is the run
		# metadata, written after _summary.tsv.
		run_unless_done "$CLI_OUT/subgen_operon_unit_jaccard_run_metadata.json" \
			"$COH/subgen_operon_unit_jaccard.py" \
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
	echo "Unknown stage '$STAGE'. Use: extract | cached | independent | cli" >&2
	exit 1
	;;
esac

subgen_summary
