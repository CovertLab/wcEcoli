# Scripts used to run analyses that generate the published figures in the rRNA
# paper
# You should run the simulations in runscripts/paper/rrna_paper/paper_runs.sh
# prior to running these analyses. You will need to edit the timestamps in the
# directory names.


## Figure 1
# Panel B
python models/ecoli/analysis/multigen/rrna_dynamics.py \
out/20240219.172845__SET_A3_24_gens_32_seeds_WT_with_shift_from_minimal_to_rich_media


## Figure 2 (WT vs single-gene-TUs)
# Panel A
# Variant scheme

# Panel B
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20240220.124636__SET_B2_16_gens_32_seeds_monocistronic_rRNA_TUs_with_rich_media \
out/20240219.172825__SET_A2_16_gens_32_seeds_WT_with_rich_media/

# Panel C
python models/ecoli/analysis/comparison/active_ribosome_counts_histogram.py \
out/20240220.124636__SET_B2_16_gens_32_seeds_monocistronic_rRNA_TUs_with_rich_media/ \
out/20240219.172825__SET_A2_16_gens_32_seeds_WT_with_rich_media/

# Panel D
python models/ecoli/analysis/comparison/rrna_to_ribosome_yield.py \
out/20240220.124636__SET_B2_16_gens_32_seeds_monocistronic_rRNA_TUs_with_rich_media/ \
out/20240219.172825__SET_A2_16_gens_32_seeds_WT_with_rich_media/

# Panel E
python models/ecoli/analysis/comparison/instantaneous_doubling_times.py \
out/20240220.124659__SET_B3_24_gens_32_seeds_monocistronic_rRNA_TUs_with_shift_from_minimal_to_rich_media/ \
out/20240219.172845__SET_A3_24_gens_32_seeds_WT_with_shift_from_minimal_to_rich_media/


## Figure S1 (associated with Figure 2)
# Panel A
python models/ecoli/analysis/comparison/doubling_time_histogram.py \
out/20240220.124611__SET_B1_16_gens_32_seeds_monocistronic_rRNA_TUs_with_glucose_minimal_media \
out/20240219.172800__SET_A1_16_gens_32_seeds_WT_with_glucose_minimal_media

# Panel B
python models/ecoli/analysis/comparison/active_ribosome_counts_histogram.py \
out/20240220.124611__SET_B1_16_gens_32_seeds_monocistronic_rRNA_TUs_with_glucose_minimal_media \
out/20240219.172800__SET_A1_16_gens_32_seeds_WT_with_glucose_minimal_media

# Panel C
python models/ecoli/analysis/comparison/rrna_to_ribosome_yield.py \
out/20240220.124611__SET_B1_16_gens_32_seeds_monocistronic_rRNA_TUs_with_glucose_minimal_media \
out/20240219.172800__SET_A1_16_gens_32_seeds_WT_with_glucose_minimal_media


## Figure 3 (knockouts)
# Panel A
# Histogram

# Panel B
# Variant scheme

# Panel C
python models/ecoli/analysis/variant/doubling_time_histogram.py \
out/20240221.112212__SET_C1_16_gens_32_seeds_rRNA_operon_knockouts_with_glucose_minimal_media/

# Panel D
python models/ecoli/analysis/variant/doubling_time_histogram.py \
out/20240221.112243__SET_C2_16_gens_32_seeds_rRNA_operon_knockouts_with_rich_media/

# Panel E
python models/ecoli/analysis/variant/active_ribosome_counts_histogram.py \
out/20240221.112212__SET_C1_16_gens_32_seeds_rRNA_operon_knockouts_with_glucose_minimal_media/

# Panel F
python models/ecoli/analysis/variant/active_ribosome_counts_histogram.py \
out/20240221.112243__SET_C2_16_gens_32_seeds_rRNA_operon_knockouts_with_rich_media/

# Panel G
