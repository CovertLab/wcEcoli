# Script to run analysis to generate published figures
# When reproducing, /scratch/PI/mcovert/wc_ecoli/paper will be the out directory
# Simulations in runscripts/paper/paper_runs.sh must be run to generate the data before performing analysis


## Figure 1
#Bottom Simulations Left to Right, Top to bottom:
 python models/ecoli/analysis/multigen/growth_dynamics_grid.py \
 /scratch/PI/mcovert/wc_ecoli/paper/SET_B/20190825.220728__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
 --variant_index 2 --seed 7


## Figure 2
# Panel a
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190821.000038__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20190822.092543__SET_D2_4_gens_256_seeds,_unfit_ribosome_expression/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20190822.230154__SET_D3_4_gens_256_seeds,_unfit_rna_poly_expression/

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_D/20190824.144321__SET_D5_4_gens_256_seeds,_unfit_ribosome_and_rna_poly_expression/

# Panel b
python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190821.000038__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/ \
--variant_index 1

python models/ecoli/analysis/cohort/doubling_times_histogram_all.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190821.000038__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/ \
--variant_index 2

# Panel c
python models/ecoli/analysis/variant/growth_condition_comparison_validation.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190821.000038__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/


## Figure 3
# panel A
# ???

# panel B-D: experimental data

# panel E
# ???

# panel F
python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/ \
--variant 0

# panel G
python models/ecoli/analysis/variant/flux_sensitivity.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_H/20190822.204306__SET_H_1_gen_flux_sensitivity/

# panel H
python models/ecoli/analysis/variant/kinetic_objective_interactions.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/

# panel I
python models/ecoli/analysis/variant/kinetic_objective_comparison.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/

# panel J
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/ \
--variant 167 --seed 0 --gen 0

# panel K
python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/ \
--variant 167

# panel L
python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/groups/mcovert/wc_ecoli/paper/SET_I/20190910.155355__SET_I_kinetic_constraint_factorial_design/ \
--variant 167


## Figure 4
# panel A
python models/ecoli/analysis/multigen/proteinCountsValidation.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190819.181225__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 0

# panel B
python models/ecoli/analysis/cohort/expression_dynamics.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190819.181225__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/

# panel C and D
python models/ecoli/analysis/multigen/subgenerational_transcription.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190819.181225__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 2

# panel E
python models/ecoli/analysis/multigen/functionalUnits.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190819.181225__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 2

# panel F
python models/ecoli/analysis/multigen/pabx_limitations.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_A/20190819.181225__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/ \
--seed 2


## Figure 5
# panel A: made by JC in repo EcoliFoldChanges

# panel B: made by JC from his experimental results


## Supplemental figure 1
python models/ecoli/analysis/multigen/growth_dynamics_panel.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_B/20190825.220728__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
--variant_index 2 --seed 7

python models/ecoli/analysis/multigen/environmental_shift_fluxes.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_B/20190825.220728__SET_B_9_gens_8_seeds_shift_to_plus_AA_without_growth_noise_with_D_period/ \
--seed 7


## Supplemental figure 2
python models/ecoli/analysis/variant/adder_sizer.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190821.000038__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

python models/ecoli/analysis/variant/rna_protein_ratio_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190821.000038__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/

python models/ecoli/analysis/variant/rna_protein_ratio_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_L/20190828.123526__SET_L_4_gens_256_seeds_3_conditions_unfit_ribosome_and_rna_poly_expression/

python models/ecoli/analysis/variant/adder_sizer_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_C/20190821.000038__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period/


## Supplemental figure 3
# panel A
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0 --generation 0 --seed 0

python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0

python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 0

# panel B
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3 --generation 0 --seed 0

python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3

python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 3

# panel C
python models/ecoli/analysis/single/massFractionSummary.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6 --generation 0 --seed 0

python models/ecoli/analysis/cohort/kinetics_flux_comparison.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6

python models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/ \
--variant_index 6

# panels D, E and F
python models/ecoli/analysis/variant/metabolism_kinetic_objective_weight.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_E/20190826.093027__SET_E_8_gens_8_seeds_10_metabolism_weighting_values/


## Supplemental figure 4
python models/ecoli/analysis/variant/meneSensitivity.py \
/scratch/PI/mcovert/wc_ecoli/paper/SET_F/20180814.113407.723416__SET_F_8_gens_8_seeds_9_menE_expression_values


## Supplemental figure 5
# made by JC from his experimental results
