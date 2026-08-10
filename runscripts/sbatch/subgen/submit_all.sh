#!/usr/bin/bash
# Submit the full subgenerational pipeline with the correct dependency:
# job 02 (reads the raw extraction) waits for job 01 to finish successfully.
# Jobs 01 and 03 are independent and start immediately.
#
# Usage:  bash submit_all.sh
# Edit the OUT / --variant-range values inside each .sbatch first.

cd "$(dirname "$0")"

jid_extract=$(sbatch --parsable 01_subgen_raw_extract.sbatch)
echo "01 subgen_raw_extract       -> job $jid_extract"

jid_downstream=$(sbatch --parsable --dependency=afterok:$jid_extract 02_subgen_def5_downstream.sbatch)
echo "02 subgen_def5_downstream   -> job $jid_downstream (afterok:$jid_extract)"

jid_independent=$(sbatch --parsable 03_subgen_independent.sbatch)
echo "03 subgen_independent       -> job $jid_independent (no dependency)"
