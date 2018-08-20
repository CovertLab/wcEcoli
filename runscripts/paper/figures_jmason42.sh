#!/bin/bash

# Warning!  This file is not committed anywhere!  Please don't delete it
# (unless you are Travis and you know what you're doing).

## Figure 2

# Panel A

# Figure 2A (tentatively) is the four variants with optional fitting of
# ribosome and RNA polymerase expression, corresponding* to sets
# D1 - rib fit, RNAp fit
# D2 - rib unfit, RNAp fit
# D3 - rib fit, RNAp unfit
# D4 - rib unfit, RNAp unfit
# *The exact composition of the panels is TBD

FIG2A_SCRIPT=models/ecoli/analysis/cohort/doubling_times_histogram_all.py
SET_D_ROOT=/scratch/PI/mcovert/wc_ecoli/paper/SET_D
SET_D1=$SET_D_ROOT/20180815.085246.549352__SET_D1_4_gens_256_seeds/
SET_D2=$SET_D_ROOT/20180815.181420.368837__SET_D2_4_gens_256_seeds,_unfit_ribosome_expression/
SET_D3=$SET_D_ROOT/20180815.235510.655007__SET_D3_4_gens_256_seeds,_unfit_rna_poly_expression/
SET_D4=$SET_D_ROOT/20180817.001302.336012__SET_D4_4_gens_256_seeds,_unfit_ribosome_and_rna_poly_expression/

python $FIG2A_SCRIPT $SET_D1
python $FIG2A_SCRIPT $SET_D2
python $FIG2A_SCRIPT $SET_D3
python $FIG2A_SCRIPT $SET_D4


# Supplementary

SUPL_SCRIPT1=models/ecoli/analysis/multigen/massFractionSummary.py

python $SUPL_SCRIPT1 $SET_D1 -s 0
python $SUPL_SCRIPT1 $SET_D2 -s 0
python $SUPL_SCRIPT1 $SET_D3 -s 0
python $SUPL_SCRIPT1 $SET_D4 -s 0
