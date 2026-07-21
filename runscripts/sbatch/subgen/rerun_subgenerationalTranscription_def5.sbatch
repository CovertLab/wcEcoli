#!/usr/bin/bash
#SBATCH --job-name=rerun_subgenTranscription_def5
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aj0204@stanford.edu
#SBATCH -p mcovert
#SBATCH -c 1
#SBATCH --mem=64GB
source /home/groups/mcovert/pyenv/versions/wcEcoli3/bin/activate

# Re-run ONLY the multigen plot that failed (missing-listener-table crash on an
# incomplete generation, now fixed to skip incomplete generations). It reads
# simOut directly and does NOT depend on the raw extraction, so nothing else
# needs to be re-run. Figures are overwritten in place.
REPO=$HOME/repos/wcEcoli
OUT=$SCRATCH/wcEcoli_out/20260706.223256__subgen_paper_sims_1
VARIANTS="-v 0"

cd $REPO
export PYTHONPATH=$PWD
export OPENBLAS_NUM_THREADS=1

python $REPO/models/ecoli/analysis/multigen/subgenerationalTranscription_def5.py $VARIANTS $OUT
