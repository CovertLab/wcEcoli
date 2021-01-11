"""
Removes amino acid inhibition feedback for single amino acids to match data in
Sander et al. Allosteric Feedback Inhibition Enables Robust Amino Acid
Biosynthesis in E. coli by Enforcing Enzyme Overabundance. 2019.

Associated variant analysis:
	remove_aa_inhibition: to produce plot similar to Fig 1B from Sander et al.

Modifies:
	sim_data.process.metabolism.aa_kis

Expected variant indices (dependent on sorted order of sim_data.conditionActiveTfs):
	0: wildtype
	1-7: enzyme mutants that show no feedback control
		Arg (argA)
		His (hisG)
		Ile (ilvA)
		Leu (leuA)
		Pro (proB)
		Thr (thrA)
		Trp (trpE)
"""

import numpy as np


AA_TO_ADJUST = ['ARG[c]', 'HIS[c]', 'ILE[c]', 'LEU[c]', 'PRO[c]', 'THR[c]', 'TRP[c]']


def remove_aa_inhibition(sim_data, index):
	if index > 0:
		metabolism = sim_data.process.metabolism
		aa = AA_TO_ADJUST[index-1]
		aa_index = metabolism.aa_aas.index(aa)
		metabolism.aa_kis[aa_index] *= np.inf
		short = aa
		desc = f'remove {aa} inhibition'
	else:
		short = 'wt'
		desc = 'wildtype with no enzyme adjustment'

	return dict(shortName=short, desc=desc), sim_data
