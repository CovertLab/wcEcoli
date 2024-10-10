Addition of new genes to the chromosome<br>(in progress)
---

<b>Directions for adding a new gene</b><br>

* Add one directory per genome insertion to (`reconstruction/ecoli/flat/new_gene_data`)
  * Begin by copying the `reconstruction/ecoli/flat/new_gene_data/template` subdirectory.  
  * Each directory can have multiple new genes (represented as multiple rows in the tables), but these genes must all be inserted in the same place in the genome.
  * <b>ID of new genes, RNAs, and proteins MUST begin with NG</b>
  * Data files that must be modified:
    * `insertion_location.tsv` - Specify desired genome insertion location
    * `genes.tsv` - List of the genes you would like to insert, with relative coordinates within the sequence to be inserted. Note that only a single contiguous gene insertion is supported at this time.
    * `gene_sequences.tsv` - Sequences of genes to be inserted
    * `rnas.tsv` - RNAs corresponding to those genes
    * `proteins.tsv` - Proteins corresponding to those genes
  * Currently supported optional data files:
    * `rna_half_lives.tsv` - Can specify if desired, otherwise will default to average of the other RNAs
    * `protein_half_lives_measured` - Can specify if desired, otherwise will default to average of the other proteins
    * `metabolic_reactions_external` - Can specify if an external metabolic pathway uses the genes added
    * `metabolites` - Can specify is an external metabolic pathway uses the genes added

With `new_genes_option != 'off'`, `KnowledgeBaseEcoli` uses the information in the new gene subdirectory, which is specified in the command line via `new_genes_option`, to make the addition of the gene to the Ecoli chromosome.
For example: `python runscripts/manual/runParca.py --new-genes 'gfp'` will add the genes from the `new_genes_data/gfp` subdirectory to the Ecoli chromosome and then run the ParCa. 
Note that even though the gene is added, it will have no expression (i.e. is knocked-out) unless a new
expression level and translation efficiency are set using `models/ecoli/sim/variants/new_gene_internal_shift.py`.

The following steps occur in incorporating the new genes into the chromosome:
* Read in all data tables from the specified `<new_genes_option>` (e.g. 'gfp'), and join rows with information about the new genes to the corresponding tables with information about the original genes.
  * For example, the rows in `reconstruction/ecoli/flat/new_gene_data/<new_genes_option>/rnas.tsv` are read in and added to the data attribute for `reconstruction/ecoli/flat/rnas.tsv`.
  * If any of the optional files are not provided, then default values will be filled in at a later point in the simulation.
* Based upon the global coordinates in `insertion_location.tsv`, check for conflicts and move insertion location as necessary. For example, if the user specified insertion location would be in the middle of another gene or transcription unit, instead move the insertion location to be immediately after that gene or transcription unit.
* Using the relative coordinates in `genes.tsv` and sequences in `gene_sequences.tsv`, determine the sequence to be inserted.
* Insert the sequence into the reference genome at the insertion location. Update the new gene relative coordinates to global coordinates and add to the orginal genes attribute. Update the coordinates of original genes and transcription units to reflect the state of the genome after the insertion.

---
<b>Variants </b><br>

* New Gene (Expression, Translation Efficiency, and Induction/Knockout)
  * `models/ecoli/variants/new_gene_internal_shift.py`
  * Variant index specifies media condition,
    factor to multiply the expression level of all new genes by, and the 
    value to use for translation efficiency for all new genes. User must provide a 
    list of the new gene 
    expression factors and the translation efficiency values, and all pairs 
    from those lists will be run. The mapping between variant indices and 
    the values for the new gene expression and translation are saved to the 
    simulation metadata so that they may be used in analysis.
    * Expected variant indices (int, positive):
      * z > 0: converted to an index for media condition, index for a list of
          new gene expression variant factors, and an index for a list of new gene
          translation efficiency values 
        * condition_index = z div 1000 
        * new_gene_index = z - condition_index * 1000
      * Expected condition_index values
        * (dependent on sim_data.ordered_conditions and should be the same order as rows in `reconstruction/flat/condition/condition_defs.tsv`)
        * 0: control (minimal media)
        * 1: with amino acids 
        * 2: acetate 
        * 3: succinate 
        * 4: minimal media (anaerobic)
      * Expected new_gene_index values
        * 0: control (knockout new gene expression)
        * 1000 > y > 0: converted to an index for a list of
            new gene expression variant factor and an index for a list of new gene
            translation efficiency values 
          * separator = number of translation efficiency values to try 
          * expression_index = y div separator + 1
          * translation_efficiency_index = y mod separator
    * Example: 
      * `NEW_GENE_EXPRESSION_FACTORS = [0, 7, 8]`
        * `NEW_GENE_EXPRESSION_FACTORS[expression_index]` specifies the factor to multiply the expression level of all new genes by
            * 0: no new gene expression
            * 1: factor = 10^0 = 1
            * x > 1: factor = 10^(x-1)
      * `NEW_GENE_TRANSLATION_EFFICIENCY_VALUES = [2.5, 1, 0.5]`
        * `NEW_GENE_TRANSLATION_EFFICIENCY_VALUES[translation_efficiency_index]` will be used as the translation efficiencies of all new genes
      * `SEPARATOR = len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES) = 3`
      * Variant Index: (Media Condition, New Gene Expression Factor, Translation Efficiency)
        * 0: (minimal media, 0, 0), 1: (minimal media, 7, 1), 2: (minimal media, 7, 0.5), 3: (minimal media, 7, 2.5), 4: (minimal media, 8, 1), 5: (minimal media, 8, 
          0.5), 6: (minimal media, 8, 2.5)
        * 1000: (with amino acids, 0, 0), 1001: (with amino acids, 7, 1), 1002: (with amino acids, 7, 0.5), 1003: (with amino acids, 7, 2.5), 1004: (with amino acids, 8, 1), 1005: (with amino acids, 8, 
          0.5), 1006: (with amino acids, 8, 2.5)
  * You can make the changes to 
    the new genes at the beginning of a specified daughter generation. You can 
    specify a `NEW_GENE_INDUTION_GEN` and `NEW_GENE_KNOCKOUT_GEN` in 
    `models/ecoli/variants/new_gene_expression_and_translation_efficiency_internal_shift.py`
    * From generations [0, `NEW_GENE_INDUCTION_GEN`), the new genes will not be 
      transcribed or translated because new genes are knocked out by default in wildtype 
      simulation.
    * From generations [`NEW_GENE_INDUCTION_GEN`, `NEW_GENE_KNOCKOUT_GEN`), 
      the new genes will be transcribed and translated using the expression 
      factor and translation efficiency value from the variant index.
    * From generations `NEW_GENE_KNOCKOUT_GEN` and onwards, new 
      gene expression probabilities will be set to 0 corresponding to new 
      gene knockout.
    * If you don't intend to do new gene induction and/or knockout, you can 
      set `NEW_GENE_INDUTION_GEN = -1` and/or `NEW_GENE_KNOCKOUT_GEN = -1`, 
      respectively.
    * Note: `NEW_GENE_INDUTION_GEN` and `NEW_GENE_KNOCKOUT_GEN` cannot be set to 0.
      * New genes must be induced after the first generation to establish an accurate shift.
      * New genes are knocked out by default, so induction should happen before knockout.
    * Note: if the values you choose for `NEW_GENE_INDUTION_GEN` 
      and `NEW_GENE_KNOCKOUT_GEN` are greater than the number of 
      generations you run, then you won't see their corresponding effects.
  * If you'd like different behavior (e.g. knock in, knock out, then 
    knock in again) you can define your own variant using the 
    internal shift variant framework described in 
    `docs/misc/internal_shift.md`.

---
<b>Listeners</b><br>
* Listeners `models/ecoli/listeners/mRNA_counts.py` and `models/ecoli/listeners/protein_counts.py` did not have to be modified
* Listener `models/ecoli/listeners/ribosome_data.py` and process 
  `models/ecoli/processes/polypeptide_initiation.py` were modified to 
  include a physical limit to the number of ribosome initiation events on 
  an mRNA, based on ribosome footprint size.

---
<b>Analysis scripts</b><br>
* New Gene Expression
  * `models/ecoli/analysis/single/new_gene_counts.py` creates two plots - one 
    with the mRNA counts for each new gene in the simulation, and one with the protein counts for each new gene in the simulation
  * `models/ecoli/analysis/multigen/new_gene_counts.py` creates two plots - 
    one with the mRNA counts for each new gene in the simulation, and one with the protein counts for each new gene in the simulation (extention for multiple generations)
  * `models/ecoli/analysis/variant/new_gene_counts.py` creates histograms and 
    scatterplots for each new gene - one with the mRNA counts for that new gene, and one with the protein counts for that new gene, both are colored by variant index
  * `models/ecoli/analysis/variant/doubling_time_histogram.py` creates two plots - one with a histogram of the doubling time, and one with the proportion of seeds that successfully reached the maximum generation in the simulation, both colored by variant index
  * `models/ecoli/analysis/variant/active_ribosome_counts_histogram.py` creates a histogram of the active ribosome counts, colored by variant index
  * `models/ecoli/analysis/variant/rnap_counts_histogram.py` creates a histogram of the RNA polymerase counts, colored by variant index
  * `models/ecoli/analysis/variant/ppgpp_concentration_histogram.py` creates a histogram of the ppGpp concentration, colored by variant index
  * `models/ecoli/analysis/variant/new_gene_protein_mass_fraction_histogram.py` creates a histogram of the proportion of total protein mass that is accounted for by new gene proteins, colored by variant index

Note: for the variant scripts, the average value for each generation is plotted, except for the plot for active ribosomes, where the initial count from each generation is plotted. These variant scripts can be used to analyze the impact that increasing new gene expression level has on the cell. In some of these scripts, you can decide whether to exlcude generations that reached the maximum simulation time. In addition, there is an option in some scripts for three separate figures to be created to encompass all generations, early generations (0-3), and late generations (4 and onwards). It is recommended to reference the late generations in your analysis, as the early generations may be impacted by the initialization process and may not be the most representative.

* New Gene Expression and Translation Efficiency
  * `models/ecoli/analysis/variant/new_gene_translation_efficiency_heatmaps.py` plots a number of heatmaps, where each square in the heatmap 
    represents an average over all seeds and generations for that variant 
    index (i.e. combination of translation efficiency value and new gene 
    expression factor). The heatmaps are colored by the magnitude of the 
    values. The text number displayed on each square can optionally be 
    colored based on whether some percentage of seeds reached the full 
    number of generations for that variant. Each heatmap displays different 
    data pertaining to new genes or cell stress markers, the options (so far) 
    are:
    * Percent of simulation seeds that successfully reached a given generation 
      number
    * Average doubling time
    * Average cell volume, mass, dry cell mass, mRNA mass, and protein mass
    * Average new gene mRNA count
    * Average new gene mRNA mass fraction
    * Average new gene NTP mass fraction
    * Average new gene protein count
    * Average new gene protein mass fraction
    * Average new gene initialization rate for RNA polymerase (RNAP) and 
      ribosomes
    * Average fraction of time new gene is overcrowded by RNAP and ribosomes
    * Average number of overcrowded genes for RNAP and ribosomes
    * Average number of ribosomes
    * Average number of RNA polymerases
    * Average ppGpp concentration

Note: for the variant scripts, the average value for each generation is 
plotted. These variant scripts can be used to analyze the impact that 
increasing new gene expression level or translation efficiency value has on 
the cell. In addition, there is an option to specify the minimum and 
maximum generation index to be plotted. It is recommended to reference the 
late generation plots in your analysis, as the early generations may be impacted by the initialization process and may not be the most representative.

* New Gene Expression and Translation Efficiency Internal Shift
  * (In progress for a future PR)

 ---

<b>Sample Commands</b><br>

`python runscripts/manual/runParca.py --new-genes gfp`

`python runscripts/manual/runSim.py --variant new_gene_internal_shift 0 10 --generations 8`   
Note: The numbers after new_gene must be integers. A different simulation will run for each value between the first and second specified numbers (inclusive). Note that for higher numbers, you may see `RuntimeError: GLP_EFAIL: Solver failure`. This is a representation of cell death, as this failure usually arises due to strain from insufficient cellular resources.


`python models/ecoli/analysis/multigen/new_gene_counts.py` 

`python models/ecoli/analysis/variant/new_gene_counts.py`

`python models/ecoli/analysis/variant/doubling_time_histogram.py`

`python models/ecoli/analysis/variant/new_gene_translation_efficiency_heatmaps.py`

---
<b>Addition of new metabolic genes<br>(in progress)

* In the genome insertion directory ('reconstruction/ecoli/flat/new_gene_data/gene_names') populate the files: 

  * metabolic_reaction_external.tsv - list the reactions from the additional pathway. Specify for each reaction their:
  
    * id (e.g. ‘rxnpath1’),  

    * Stoichiometry represented as a dictionary where the keys are the reactants and the products, and the values represent the coefficients for each reactant (positive coefficient) or product (negative coefficient) 

    * Direction for reading the reaction (e.g. L2R) 

    * A list containing the name of the enzyme catalyzing the reactions. At the moment each reaction can only be catalyzed by 1 enzyme. 

    * A list containing the name of the substrate of the reaction. At the moment, each reaction can only have one substrate. 

    * The turnover number of an enzyme (kcat) 

    * The Michaelis Menten constant (KM) 

  * metabolites.tsv - list of the metabolites produced by the reaction. Specify for each of them: 

    * Their ID as shown on Ecocyc (e.g. CPD-11890) 

    * Their common name according to Ecocyc. 

    * A list of synonyms. This can be empty. 

    * Chemical formula (e.g. C11H9) 

    * Molecular charge. 

    * Smiles. 

  * If these files are not empty, the WCM will run the metabolic_reactions_external process. 

A new process was created to represent the metabolism corresponding to these reactions. This new process is based on Michaelis Menten kinetics and it is run at each time step before the Metabolism process. To achieve this, two Python files were added: 

  * reconstruction/ecoli/dataclasses/process/metabolism_external_pathway.py an – This script loads the raw data from the files populated above and stores it in the sim_data object. It also makes it accessible as an instance of the new process class. In this file we also define the methods that represent the Michaelis Menten model. This script is called when runParca is run. 

  * models/ecoli/metabolism_external_pathway.py is called when running the sim. The script defines a class corresponding to the process that initializes all the variables necessary, it calls the Michaelis Menten model and updates the state of the cell accordingly. Furthermore, in this sript it is possible to update any listeners with the data produced. 

Analysis scripts  

Four analysis scripts are added when there is data in the metabolic_reaction_external.tsv and metabolites.tsv: 

  * models/ecoli/analysis/single/flux_external_metabolic_pathway.py - creates a plot of the fluxes corresponding to the external reactions added. This will be one plot for each cell. 

  * models/ecoli/analysis/multigen/flux_external_metabolic_pathway.py - creates a plot of the fluxes corresponding to the external reactions added. This will show all these fluxes for all generations simulated in one variant. 

  * models/ecoli/analysis/single/molecules_external_pathway.py - creates plots of molecule counts for each molecule involved in the external reactions, for each cell separately, from every generation. 

  * models/ecoli/analysis/multigen/molecules_external_pathway.py - creates plots of molecule counts for each molecule involved in the external reactions, for all cells from all generations corresponding to a variant. 

The variants and the sample commands for running this code are the same as in the case of gene addition (without any metabolic pathway involved).  