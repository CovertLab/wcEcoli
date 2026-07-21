#!/usr/bin/bash
#SBATCH --job-name=subgen_raw_extract
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aj0204@stanford.edu
#SBATCH -p mcovert
#SBATCH -c 1
#SBATCH --mem=128GB
source /home/groups/mcovert/pyenv/versions/wcEcoli3/bin/activate

# --- Config (edit OUT to point at your cohort) --------------------------------
REPO=$HOME/repos/wcEcoli
OUT=$SCRATCH/wcEcoli_out/20260706.223256__subgen_paper_sims_1
# Variant selection (cohort analysis runs one variant at a time):
#   -v 0                -> single variant, index 0 (e.g. wildtype_000000). Use
#                          this for sims with no variant applied.
#   (omit entirely)     -> auto-selects the first/only variant.
#   --variant-range 0 2 -> variants 0,1,2 INCLUSIVE (multi-variant sims).
VARIANTS="-v 0"

cd $REPO
export PYTHONPATH=$PWD
export OPENBLAS_NUM_THREADS=1

# STEP 1: read simOut once per cell and write the per-cell matrices that every
# Definition-5 / Task-4 analysis reads. This is the ONLY heavy sim-reading step
# in the def5 pipeline; run it before job 02. Writes into each variant's plotOut.
python $REPO/models/ecoli/analysis/cohort/subgen_raw_extract.py $VARIANTS $OUT
