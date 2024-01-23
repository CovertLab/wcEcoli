# Scripts used to run analyses that generate the published figures in the operon
# paper
# You should run the simulations in runscripts/paper/operon_paper/paper_runs.sh
# prior to running these analyses. You will need to edit the timestamps in the
# directory names.


## Figure 1
# Panel B
python models/ecoli/analysis/comparison/polycistronic_transcription.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 2
# Panel B
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20230719.104706__SET_D_8_gens_128_seeds_operons_v1_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media

# Panel C
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20230719.104706__SET_D_8_gens_128_seeds_operons_v1_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media

# Panel E
# See Rend-seq repository

# Panel G
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20230722.192107__SET_F_8_gens_128_seeds_operons_v2_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media

# Panel H
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20230722.192107__SET_F_8_gens_128_seeds_operons_v2_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media


## Figure 3
# Panel A
python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20230722.192107__SET_F_8_gens_128_seeds_operons_v2_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media

# Panel F
python models/ecoli/analysis/comparison/mRNA_copy_numbers_short_genes.py \
out/20230726.175251__SET_G_8_gens_128_seeds_operons_v3_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel H (+ Figure S2A)
python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20230726.175322__SET_H_8_gens_128_seeds_operons_v3_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media

python models/ecoli/analysis/comparison/mRNA_copy_numbers.py \
out/20230822.135902__SET_J_8_gens_128_seeds_operons_on_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media

# Panels I, J
# See Rend-seq repository


## Figure 4
# Panel A
python models/ecoli/analysis/comparison/mRNA_length_histogram.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel B
python models/ecoli/analysis/comparison/mRNA_counts_histogram.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel C
python models/ecoli/analysis/comparison/mRNA_mass_histogram.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure 5
# Panel A
python models/ecoli/analysis/comparison/polycistronic_transcription_extended.py \
out/20230830.144259__SET_L_32_gens_8_seeds_operons_on_with_glucose_minimal_media \
out/20230830.144235__SET_K_32_gens_8_seeds_operons_off_with_glucose_minimal_media

# Panel B (+ Table S2)
python models/ecoli/analysis/comparison/coexpression_probabilities.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panels C, D (+ Table S3)
python models/ecoli/analysis/comparison/protein_stoichiometry.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panels E, F, G (+ Table S4)
python models/ecoli/analysis/comparison/excess_protein_monomers.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure S1
# Panel B
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20230719.104634__SET_C_8_gens_128_seeds_operons_v1_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel C
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20230719.104634__SET_C_8_gens_128_seeds_operons_v1_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel D
python models/ecoli/analysis/comparison/mRNA_copy_numbers_growth_genes.py \
out/20230722.192024__SET_E_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel E
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20230722.192024__SET_E_8_gens_128_seeds_operons_v2_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media


## Figure S2
# Panel A
# See Figure 3H


## Figure S3
# Panel A
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media

# Panel B
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20230822.135902__SET_J_8_gens_128_seeds_operons_on_with_rich_media \
out/20230717.142940__SET_B_8_gens_128_seeds_operons_off_with_rich_media

# Panels C, D
python models/ecoli/analysis/comparison/proteomics_fluxomics_comparison.py \
out/20230822.135846__SET_I_8_gens_128_seeds_operons_on_with_glucose_minimal_media \
out/20230716.121449__SET_A_8_gens_128_seeds_operons_off_with_glucose_minimal_media
