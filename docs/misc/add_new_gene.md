Addition of new genes to the chromosome<br>(in progress)
---

<b>Directions for adding a new gene</b><br>

* Add one directory per genome insertion to (`reconstruction/ecoli/flat/new_gene_data`), following file structure in the `template` subdirectory.  
  * Each directory can have multiple new genes (represented as multiple rows in the tables), but these genes must all be inserted in the same place in the genome.
  * <b>ID of new genes, RNAs, and proteins MUST begin with NG</b>
  * Necessary data files:
    * `insertion_location.tsv` - Specify desired genome insertion location
    * `genes.tsv` - List of the genes you would like to insert, with relative coordinates within the sequence to be inserted
    * `gene_sequences.tsv` - Sequences of genes to be inserted
    * `rnas.tsv` - RNAs corresponding to those genes
    * `proteins.tsv` - Proteins corresponding to those genes
    *  `rnaseq_rsem_tpm_mean.tsv` - Best practices: set these each to .01 (so that the Parca does not try to correct for the new genes) and modify basal expression levels later via variant (more details below)
  * Currently supported optional data files:
    * `rna_half_lives.tsv` - Can specify if desired, otherwise will default to average of the other RNAs
    * `protein_half_lives_measured` - Can specify if desired, otherwise will default to average of the other proteins

With `new_genes_option != 'off'`, `KnowledgeBaseEcoli` uses the information in the new gene subdirectory, which is specified in the command line via `new_genes_option`, to make the addition of the gene to the Ecoli chromosome.
For example: `python runscripts/manual/runParca.py --new-genes 'gfp'` will add the genes from the `new_genes_data/gfp` subdirectory to the Ecoli chromosome and then run the Parca.

The following steps occur in incorporating the new genes into the chromosome:
* Read in all data tables from the specified `<new_genes_option>` (e.g. 'gfp'), and join rows with information about the new genes to the corresponding tables with information about the original genes.
  * For example, the rows in `reconstruction/ecoli/flat/new_gene_data/<new_genes_option>/rnas.tsv` are read in and added to the data attribute for `reconstruction/ecoli/flat/rnas.tsv`.
  * If any of the optional files are not provided, then default values will be filled in at a later point in the simulation.
* Based upon the global coordinates in `insertion_location.tsv`, check for conflicts and move insertion location as necessary. For example, if the user specified insertion location would be in the middle of another gene or transcription unit, instead move the insertion location to be immediately after that gene or transcription unit.
* Using the relative coordinates in `genes.tsv` and sequences in `gene_sequences.tsv`, determine the sequence to be inserted.
* Insert the sequence into the reference genome at the insertion location. Update the new gene relative coordinates to global coordinates and add to the orginal genes attribute. Update the coordinates of original genes and transcription units to reflect the state of the genome after the insertion.

---

* Variants
  * `models/ecoli/variants/new_gene_expression.py` - index specifies the factor to multiply the expression level of all new genes by
    * 0: no new gene expression
    * 1: factor = 10^0 = 1
    * x > 1: factor = 10^(x-1)
  * TODO: translational efficiency

---

* Listeners (`models/ecoli/listeners/mRNA_counts.py` and `models/ecoli/listeners/protein_counts.py`) did not have to be modified

---

* Analysis scripts
  * `models/ecoli/analysis/single/newGeneCounts.py` creates two plots - one with the mRNA counts for each new gene in the simulation, and one with the protein counts for each new gene in the simulation
  * `models/ecoli/analysis/multigen/newGeneCounts.py` creates two plots - one with the mRNA counts for each new gene in the simulation, and one with the protein counts for each new gene in the simulation (extention for multiple generations)
  * `models/ecoli/analysis/variant/newGeneCounts.py` creates histograms and scatterplots for each new gene - one with the mRNA counts for that new gene, and one with the protein counts for that new gene, both are colored by variant index
  * `models/ecoli/analysis/variant/doubling_time_histogram.py` creates two plots - one with a histogram of the doubling time, and one with the proportion of seeds that successfully reached the maximum generation in the simulation, both colored by variant index
  * `models/ecoli/analysis/variant/ribosome_counts_histogram.py` creates a histogram of the ribosome counts, colored by variant index
  * `models/ecoli/analysis/variant/rnap_counts_histogram.py` creates a histogram of the RNA polymerase counts, colored by variant index
  * `models/ecoli/analysis/variant/ppgpp_concentration_histogram.py` creates a histogram of the ppGpp concentration, colored by variant index
  * `models/ecoli/analysis/variant/new_gene_protein_mass_fraction_histogram.py` creates a histogram of the proportion of total protein mass that is accounted for by new gene proteins, colored by variant index

Note: for the variant scripts, the average value for each generation is plotted. These variant scripts can be used to analyze the impact that increasing new gene expression level has on the cell. In each of these scripts, you can decide whether to exlcude generations that reached the maximum simulation time. In addition, there is an option for three separate figures to be created to encompass all generations, early generations (0-3), and late generations (4 and onwards). It is recommended to reference the late generation plots in your analysis, as the early generations may be impacted by the initialization process and may not be the most representative.

 ---

Sample Commands

`python runscripts/manual/runParca.py --new-genes gfp`

`python runscripts/manual/runSim.py --variant new_gene_expression 0 10 --generations 8`   
Note: The numbers after new_gene_expression must be integers. A different simulation will run for each value between the first and second specified numbers (inclusive). Note that for higher numbers, you may see `RuntimeError: GLP_EFAIL: Solver failure`. This is a representation of cell death, as this failure usually arises due to strain from insufficient cellular resources.


`python models/ecoli/analysis/multigen/newGeneCounts.py` 

`python models/ecoli/analysis/variant/newGeneCounts.py`

`python models/ecoli/analysis/variant/doubling_time_histogram.py`