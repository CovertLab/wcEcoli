"""
Test the units of the ribosome crowding cap in PolypeptideInitiation.

The cap on translation initiation probability is

	max_p = elongation_rate * time_step / footprint / n_ribosomes_activated

The elongation rate is in amino acids per second, so the footprint has to be
expressed in codons. It is tabulated in nucleotides (24 nt for the active
ribosome), and dividing aa/s by a footprint in nucleotides made the cap three
times too strict: it allowed one ribosome per 24 codons instead of one per 24
nucleotides, about 10 ribosomes on a 717 nt mRNA where ~30 fit.

The old expression also multiplied by units.s and then called .asNumber() on a
quantity carrying [s/nucleotide], so the stray unit was silently dropped
instead of raising. These tests pin down both the conversion and the fact that
the result is a plain dimensionless number.

See CovertLab/vEcoli#306.
"""

import unittest

from models.ecoli.processes.polypeptide_initiation import NT_PER_CODON
from wholecell.utils import units

FOOTPRINT_NT = 24           # active ribosome footprint, footprint_sizes.tsv
GFP_LENGTH_NT = 717


class TestRibosomeFootprintUnits(unittest.TestCase):

	def test_nt_per_codon(self):
		self.assertEqual(NT_PER_CODON, 3)

	def test_footprint_converts_to_codons(self):
		"""A 24 nt footprint is 8 codons, so the cap uses 8, not 24."""
		footprint_codons = FOOTPRINT_NT / NT_PER_CODON
		self.assertEqual(footprint_codons, 8)

	def test_ribosomes_that_fit_on_an_mrna(self):
		"""
		Spacing is physical, so it is set by nucleotides: a 717 nt mRNA holds
		about 29 ribosomes at a 24 nt footprint. The pre-fix cap allowed only
		a third of that.
		"""
		fits = GFP_LENGTH_NT // FOOTPRINT_NT
		self.assertEqual(fits, 29)
		pre_fix = GFP_LENGTH_NT // (FOOTPRINT_NT * NT_PER_CODON)
		self.assertEqual(pre_fix, 9)
		self.assertAlmostEqual(fits / pre_fix, 3.2, places=1)

	def test_cap_is_dimensionless(self):
		"""
		A probability must be a plain number. The pre-fix expression produced
		a Unum in [s/nucleotide] and relied on .asNumber() to discard it.
		"""
		elongation_rate = 16.0          # aa / s, as read from the listener
		time_step = 1.0                 # s
		n_activated = 100
		footprint_codons = FOOTPRINT_NT / NT_PER_CODON

		max_p = (elongation_rate * time_step / footprint_codons / n_activated)
		self.assertIsInstance(max_p, float)
		self.assertAlmostEqual(max_p, 16.0 / 8 / 100)

		pre_fix = (elongation_rate / (FOOTPRINT_NT * units.nt)
			* units.s * time_step / n_activated)
		self.assertNotIsInstance(pre_fix, float)
		self.assertEqual(pre_fix.strUnit(), '[s/nucleotide]')
		self.assertAlmostEqual(max_p / pre_fix.asNumber(), float(NT_PER_CODON))


if __name__ == '__main__':
	unittest.main()
