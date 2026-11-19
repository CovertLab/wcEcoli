"""
Tests for the GFP paper's shared data layer and panel naming.

These cover the parts that fail silently rather than loudly: decoding a
variant index into the wrong design, a ratio that divides by zero, a lineage
counted as surviving when it died, and a panel saved under the wrong label.
"""

import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from wholecell.analysis import gfp_paper_data as D
from wholecell.analysis import gfp_paper_style as style

GENE = D.GeneInfo(rna_id='NG001_RNA[c]', monomer_id='NG-GFP-MONOMER[c]',
	monomer_mass_fg=4.5e-5, length_nt=720.0, rnap_footprint_nt=50.0,
	ribosome_footprint_nt=24.0, s30_id='s30', s50_id='s50')

# The 2024-05-02 grid, whose decoding is known from its variant descriptions
OLD_GRID = {'variant': 'new_gene_internal_shift',
	'new_gene_expression_factors': [0, 7, 8, 9, 10],
	'new_gene_translation_efficiency_values': [10, 5, 1, 0.1, 0]}


class TestDesignDecoding(unittest.TestCase):
	def test_control(self):
		self.assertEqual(D.design_of('', OLD_GRID, 0), (0.0, 0.0))

	def test_known_grid_variants(self):
		# 10, 11 and 12 are medium/RBS 10, strong/RBS 5 and strong/RBS 1
		self.assertEqual(D.design_of('', OLD_GRID, 10), (8.0, 10.0))
		self.assertEqual(D.design_of('', OLD_GRID, 11), (9.0, 5.0))
		self.assertEqual(D.design_of('', OLD_GRID, 12), (9.0, 1.0))
		self.assertEqual(D.design_of('', OLD_GRID, 1), (7.0, 5.0))
		self.assertEqual(D.design_of('', OLD_GRID, 5), (7.0, 10.0))

	def test_every_grid_cell_once(self):
		factors = OLD_GRID['new_gene_expression_factors']
		efficiencies = OLD_GRID['new_gene_translation_efficiency_values']
		n = 1 + (len(factors) - 1) * len(efficiencies)
		decoded = [D.design_of('', OLD_GRID, v) for v in range(1, n)]
		expected = {(float(f), float(e)) for f in factors[1:]
			for e in efficiencies}
		self.assertEqual(set(decoded), expected)
		self.assertEqual(len(decoded), len(set(decoded)))

	def test_condition_offset_ignored(self):
		self.assertEqual(D.design_of('', OLD_GRID, 1011),
			D.design_of('', OLD_GRID, 11))

	def test_missing_metadata_raises(self):
		with self.assertRaises(ValueError):
			D.design_of('', {'variant': 'something_else'}, 3)

	def test_design_list_mismatch_raises(self):
		from models.ecoli.sim.variants.new_gene_design_list import (
			NEW_GENE_DESIGNS)
		with tempfile.TemporaryDirectory() as batch:
			meta = os.path.join(batch, 'new_gene_design_list_000003',
				'metadata')
			os.makedirs(meta)
			expression, trl_eff = NEW_GENE_DESIGNS[3]
			with open(os.path.join(meta, 'description'), 'w') as f:
				f.write('New gene design 3: expression factor %s, translation'
					' efficiency %s, induced at generation 8'
					% (expression + 1, trl_eff))
			with self.assertRaises(ValueError):
				D.design_of(batch, {'variant': 'new_gene_design_list'}, 3)
			with open(os.path.join(meta, 'description'), 'w') as f:
				f.write('New gene design 3: expression factor %s, translation'
					' efficiency %s, induced at generation 8'
					% (expression, trl_eff))
			self.assertEqual(
				D.design_of(batch, {'variant': 'new_gene_design_list'}, 3),
				(float(expression), float(trl_eff)))


class TestRatios(unittest.TestCase):
	def test_safe_ratio(self):
		out = D._safe_ratio([1.0, 2.0, 3.0, 4.0], [2.0, 0.0, -1.0, np.nan])
		np.testing.assert_array_equal(out[:1], [0.5])
		self.assertTrue(np.isnan(out[1:]).all())

	def test_derived_columns(self):
		df = pd.DataFrame({'active_ribosomes': [1000.0], 'free_ribosomes':
			[200.0], 'gfp_protein_mass_fg': [10.0], 'proteinMass': [100.0],
			'gfp_mrna': [50.0], 'all_mrna': [5000.0],
			'ribosomes_on_gfp': [400.0], 'rnap_on_gfp': [20.0],
			'gfp_promoter_copies': [2.0], 'active_rnap': [800.0]})
		out = D.add_derived_columns(df, GENE).iloc[0]
		self.assertAlmostEqual(out['ribosomes_per_gfp_mrna'], 8.0)
		self.assertAlmostEqual(out['ribosome_spacing_nt'], 90.0)
		self.assertAlmostEqual(out['ribosome_ceiling_per_gfp_mrna'], 24.0)
		self.assertAlmostEqual(out['rnap_spacing_nt'], 72.0)
		self.assertAlmostEqual(out['gfp_ribosome_share'], 0.4)
		self.assertAlmostEqual(out['gfp_rnap_share'], 0.025)
		self.assertAlmostEqual(out['gfp_rnap_share_at_spacing_limit'],
			720 / 50 * 2 / 800)
		self.assertAlmostEqual(GENE.ribosomes_when_full, 30.0)


def _cycles():
	"""Two seeds, control plus one design; design seed 1 dies at gen 17."""
	rows = []
	for variant, promoter in ((0, 0.0), (5, 8.0)):
		for seed in (0, 1):
			last = 17 if (variant, seed) == (5, 1) else 23
			for gen in range(0, last + 1):
				growth = 1.0 if variant == 0 else 0.5 + 0.1 * seed
				rows.append(dict(variant=variant, seed=seed, generation=gen,
					promoter=promoter, rbs=promoter / 8 * 5, divided=gen < last
					or last == 23, growth_rate_per_h=growth * (1 + seed),
					doubling_time_min=60.0))
	return pd.DataFrame(rows)


class TestLineages(unittest.TestCase):
	def setUp(self):
		self.cycles = _cycles()
		self.lineages = D.lineage_table(self.cycles, GENE)

	def test_survival(self):
		survived = self.lineages.set_index(['variant', 'seed'])['survived']
		self.assertTrue(survived[(0, 0)])
		self.assertTrue(survived[(5, 0)])
		self.assertFalse(survived[(5, 1)])

	def test_window_excludes_the_cycle_that_did_not_divide(self):
		row = self.lineages.set_index(['variant', 'seed']).loc[(5, 1)]
		self.assertEqual(row['n_window_cycles'], 1)        # gen 16 only

	def test_growth_is_relative_to_the_same_seed(self):
		rel = self.lineages.set_index(['variant', 'seed'])[
			'relative_growth_rate']
		# seed 1 control grows at 2.0, its design at 1.2
		self.assertAlmostEqual(rel[(5, 1)], 0.6)
		self.assertAlmostEqual(rel[(5, 0)], 0.5)

	def test_design_uses_survivors_only(self):
		designs = D.design_table(self.lineages, GENE).set_index('variant')
		self.assertTrue(designs.loc[5, 'failed'])
		self.assertEqual(designs.loc[5, 'n_lineages_survived'], 1)
		self.assertAlmostEqual(designs.loc[5, 'relative_growth_rate'], 0.5)
		self.assertFalse(designs.loc[0, 'failed'])


class TestPanelName(unittest.TestCase):
	def test_label_from_manifest(self):
		with tempfile.NamedTemporaryFile('w', suffix='.tsv',
				delete=False) as f:
			f.write('panel\tcontent\tscript\n')
			f.write('3F\tx\tmodels/ecoli/analysis/variant/gfp_paper_foo.py\n')
			f.write('2B\ty\tmodels/ecoli/analysis/variant/gfp_paper_bar.py\n')
			f.write('2C\ty\tmodels/ecoli/analysis/variant/gfp_paper_bar.py\n')
			manifest = f.name
		try:
			self.assertEqual(style.panel_name('/x/gfp_paper_foo.py',
				manifest=manifest), 'Figure_3F_foo')
			self.assertEqual(style.panel_name('gfp_paper_foo.py', 'growth',
				manifest=manifest), 'Figure_3F_foo_growth')
			with self.assertRaises(ValueError):
				style.panel_name('gfp_paper_bar.py', manifest=manifest)
			with self.assertRaises(ValueError):
				style.panel_name('gfp_paper_missing.py', manifest=manifest)
		finally:
			os.remove(manifest)

	def test_every_figure_3_script_is_in_the_real_manifest(self):
		for script in ('gfp_paper_transcriptome_scatter.py',
				'gfp_paper_gfp_mrna_growth_rbs_0.py',
				'gfp_paper_rnap_spacing_on_gfp.py',
				'gfp_paper_ribosome_spacing.py',
				'gfp_paper_ribosomes_per_gfp_mrna_vs_pool_ceiling.py',
				'gfp_paper_share_of_active_rnap_vs_ribosomes.py',
				'gfp_paper_growth_vs_gfp_ribosome_share_per_cycle.py'):
			self.assertTrue(style.panel_name(script).startswith('Figure_3'))


class TestRbsColour(unittest.TestCase):
	def test_fixed_per_value(self):
		self.assertEqual(style.rbs_color(5), style.rbs_color(5.0))
		self.assertNotEqual(style.rbs_color(1), style.rbs_color(10))
		self.assertEqual(style.rbs_color(0), style.RBS_ZERO_COLOR)


if __name__ == '__main__':
	unittest.main()


class TestFigure5Helpers(unittest.TestCase):
	def test_ties_split_credit(self):
		# step 1: protein 0, rRNA 0 (tie); step 2: rRNA alone limits;
		# step 3: protein alone (needs 2 copies, has 2 -> 1 per copy)
		free = np.array([[0, 0], [5, 1], [2, 3]])
		credit, tied = D.split_limiting_credit(free, [2, 1])
		np.testing.assert_allclose(credit, [(0.5 + 0 + 1) / 3,
			(0.5 + 1 + 0) / 3])
		self.assertAlmostEqual(tied, 1 / 3)
		self.assertAlmostEqual(credit.sum(), 1.0)

	def test_replichore_fraction(self):
		class Stub:
			pass
		sim = Stub()
		sim.process = Stub()
		sim.process.replication = Stub()
		sim.process.replication.replichore_lengths = (100, 50)
		np.testing.assert_allclose(D.replichore_fraction(sim,
			[0, 50, 100, -25, -50]), [0, 0.5, 1, 0.5, 1])

	def test_seed_matched_ratio_skips_unmatched_and_dead(self):
		lineages = pd.DataFrame({'variant': [0, 0, 5, 5, 5],
			'seed': [0, 1, 0, 1, 2], 'survived': [True, True, True, False,
			True], 'x': [2.0, 4.0, 1.0, 1.0, 9.0]})
		r = D.seed_matched_ratio(lineages, 5, 'x')
		self.assertEqual(r['n'], 1)           # seed 1 died, seed 2 unmatched
		self.assertAlmostEqual(r['mean'], 0.5)

	def test_nearest_design(self):
		designs = pd.DataFrame({'variant': [0, 1, 2, 3], 'promoter':
			[0, 8, 8, 9], 'n_lineages_survived': [5, 5, 5, 0],
			'gfp_proteome_share': [0.0, 0.10, 0.30, 0.35]})
		self.assertEqual(D.nearest_design(designs, 0.35)['variant'], 2)
		self.assertEqual(D.nearest_design(designs, 0.0)['variant'], 1)

	def test_pooled_response(self):
		rows = []
		for seed in (0, 1):
			for event in (0, 1):
				rows.append(dict(seed=seed, generation=16, event=event,
					window='pre', window_lo=-2, window_hi=0, n_steps=100,
					x=1.0))
				rows.append(dict(seed=seed, generation=16, event=event,
					window='post_0_2', window_lo=0, window_hi=2,
					n_steps=50, x=1.0 + seed))
		r = D.pooled_response(pd.DataFrame(rows), 'x', n_boot=200).iloc[0]
		# pre rate 4/400; post rate (2 + 4)/200 -> ratio 3
		self.assertAlmostEqual(r['ratio'], 3.0)
		self.assertEqual(r['n_events'], 4)
		self.assertLessEqual(r['lo'], r['ratio'])
		self.assertGreaterEqual(r['hi'], r['ratio'])
