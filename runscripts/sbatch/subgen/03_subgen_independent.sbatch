#!/usr/bin/bash
#SBATCH --job-name=subgen_independent
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aj0204@stanford.edu
#SBATCH -p mcovert
#SBATCH -c 1
#SBATCH --mem=128GB
source /home/groups/mcovert/pyenv/versions/wcEcoli3/bin/activate

# --- Config -------------------------------------------------------------------
REPO=$HOME/repos/wcEcoli
OUT=$SCRATCH/wcEcoli_out/20260706.223256__subgen_paper_sims_1
# -v 0 = single/no-variant sim; --variant-range 0 2 = variants 0,1,2 inclusive.
VARIANTS="-v 0"

cd $REPO
export PYTHONPATH=$PWD
export OPENBLAS_NUM_THREADS=1

# These read simOut directly and do NOT depend on the raw extraction, so they can
# run in parallel with jobs 01/02. The CI classifier + the data-collection /
# Task-1 scripts (now with the strict lineage filter + bug fixes).
COH=$REPO/models/ecoli/analysis/cohort
for script in \
    subgen_definition5_lineage_ci \
    transcription_prob_subgen \
    subgen_peak_counts \
    extract_subgen_monomer_counts \
    subsampl_mrna_protein_10k \
    subsample_timepoints_for_scRNAseq_comp \
    subsample_cell_monomer_mass \
    subgen_monomer_dynamics ; do
    echo "=== Running $script ==="
    python $COH/$script.py $VARIANTS $OUT
done

# Multigen Fig-5 def5 figures (reads simOut directly; runs on seed 0 per variant).
python $REPO/models/ecoli/analysis/multigen/subgenerationalTranscription_def5.py $VARIANTS $OUT
