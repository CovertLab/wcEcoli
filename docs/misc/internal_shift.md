Internal Shift Variants
---
Developed by @rjuenemann in collaboration with @ggsun and @1fish2 .

<b>Motivation</b><br>
* As described in Issue #1376, we previously added an option to the ParCa 
  to add new genes to the E. coli chromosome. We ran large 
  batches of simulations on Sherlock to investigate the impact of promoter 
  and ribosome binding site strength of new genes on product production and 
  cell health. To do so in a more automated/scaled way, we developed a 
  variant that will generate the different parameter combinations and 
  modify the relevant sim data accordingly 
  (`models/ecoli/sim/variants
  /new_gene_expression_and_translation_efficiency.py`). We then wanted to
  run simulations where new gene expression is turned off (eg `sim_data.
  process.transcription.rna_synth_prob = 0` and `sim_data.process.
  translation.translation_efficiencies_by_monomer = 0` for the new genes) 
  for the first few generations (say 8) and then will be turned on (eg 
  `sim_data.process.transcription.rna_synth_prob = 0.01` and `sim_data.
  process.translation.translation_efficiencies_by_monomer = 1`) for the 
  next generations. This would allow us to observe what happens in the 
  transition when new gene expression is induced.
* @tahorst has run transition variants before. One example is 
  `models/ecoli/sim/variants/add_one_aa_shift.py`. These shifts were 
  time-based, and modify components of the external_state via timelines, 
  which cannot be easily extended to `sim_data` not in `external_state`.
* A new "internal shift" framework was developed to handle such use cases.

---
<b>How to Create an Internal Shift Variant </b><br>
* Template: `models/ecoli/sim/variants/template_internal_shift.py` 
  (duplicate, rename, and then modify)
* <b> Please ensure `_internal_shift` is included at the end of the variant name.</b>
* Internal shift variants work by defining a dictionary (which is saved 
  to `sim_data`) specifying what changes need to happen at the beginning of 
  which generations.
  * When a daughter cell is initialized, the generation index of the 
    daughter cell will be compared against the generation index 
    keys in the dictionary. The generation index keys in the 
    dictionary correspond to the starting index of the range of 
    generation indices where the change function(s) (the values of the 
    dictionary) will be applied.
  * Example: We are running an 8 generation simulation. We want 
    change 1 to be applied to generation indices 4-5 and change 
    2 to be applied to generation indices 6-7. (Remember we start counting at 
    generation 0, so generation index 7 is really the 8th cell run.)
    * First, define functions change1 and change2 that take as input the 
      variant index (i) and modify sim_data as desired in change 1 and 
      change 2, respectively.
    * Then, fill in `sim_data.internal_shifts.internal_shift_dict` to 
      look like: `{4: [(change1, i)], 6: [(change2, i)]}`
    * If we would also like to make change 3 before change 1 to generation 5 
      only, define a change3 function that takes as input the variant index (i) and 
      modifies sim_data as desired in change 3. Then modify `sim_data.internal_shifts.internal_shift_dict` to 
      look like: `{4: [(change1, i)], 5: [(change3, i), (change1, i)], 6: [(change2, i)]}`.
* Aside from defining this dictionary, internal shift variants work the 
  same as other variants! The change functions you define will look similar to 
  functions you define for normal variants.
* Don't forget to update header and inline comments for your variant.
* No other model changes need to be made outside of the file specifying 
  the variant. The model has been modified such that when initializing each 
  daughter cell it will check `sim_data.internal_shifts.internal_shift_dict` and apply any functions that are there to modify 
  `sim_data` before beginning the simulation.
* You may refer to `models/ecoli/sim/variants/new_gene_expression_and_translation_efficiency_internal_shft.py` as an 
  example implementation.

---
<b>Reasons for Implementation Decisions </b><br>
* Works within existing variant pipelines.
* Generalizable to multiple shifts and for applications outside of new gene 
  chromosome insertion.
* Defining `sim_data.internal_shifts.internal_shift_dict` in `sim_data` 
  means it is accessible to analysis scripts without the need to generate 
  more metadata files, save multiple copies of `sim_data`, etc.
* Chose a generation based shift instead of a time based shift because...
  * Ease of analysis.
  * Did not (currently) have a use case for multiple internal shifts per 
    generation, 
    especially since there is some time delay between shift introduction and 
    visible changes in model output.
* Note: The simulation initialization order needed to be changed in order 
  for this to work as desired. See Pull Request #1395 for details.

