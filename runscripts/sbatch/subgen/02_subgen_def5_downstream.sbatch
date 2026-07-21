#!/usr/bin/bash
#SBATCH --job-name=subgen_def5_downstream
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aj0204@stanford.edu
#SBATCH -p mcovert
#SBATCH -c 1
#SBATCH --mem=32GB
source /home/groups/mcovert/pyenv/versions/wcEcoli3/bin/activate

# --- Config (must match job 01) -----------------------------------------------
REPO=$HOME/repos/wcEcoli
OUT=$SCRATCH/wcEcoli_out/20260706.223256__subgen_paper_sims_1
# MUST match job 01's variant selection. -v 0 = single/no-variant sim;
# --variant-range 0 2 = variants 0,1,2 inclusive.
VARIANTS="-v 0"

cd $REPO
export PYTHONPATH=$PWD
export OPENBLAS_NUM_THREADS=1

# STEP 2: these all read the raw extraction from each variant's plotOut, so they
# are cheap (no sim reads except monomer_dynamics_def5, which pulls time-series
# traces). SUBMIT WITH: --dependency=afterok:<job01_id> so they run only after
# the extraction finishes.
COH=$REPO/models/ecoli/analysis/cohort
for script in \
    subgenerational_expression_table_def5 \
    id_subgen_monomers_def5 \
    transcriptFrequency_def5 \
    subgen_monomer_dynamics_def5 \
    subgen_per_generation_expression ; do
    echo "=== Running $script ==="
    python $COH/$script.py $VARIANTS $OUT
done
