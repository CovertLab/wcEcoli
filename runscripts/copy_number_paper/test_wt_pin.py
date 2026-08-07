#!/usr/bin/env python
"""
Verify that PIN_WT_COORDINATE does what it claims, without running a simulation.

The whole mechanism lives in one quantity: the `factor` returned by
`transcription.synth_prob_from_ppgpp`, which is

	factor = (growth + deg_rate) / get_average_copy_number(tau, wt_coordinate)

For the construct, `growth` and `deg_rate` are identical across insertion
positions, so `factor` depends on nothing but the wildtype coordinate. That
makes it an exact test rather than a statistical one.

	unpinned   construct factor at P1 / at P4  ==  0.753    (position cancels)
	pinned     construct factor at P1 / at P4  ==  1.0000   (position survives)

The unpinned figure is tau-dependent and 0.753 is its value at the ParCa's
FITTED doubling time of 44.0 min, which is what this script evaluates at. The
0.787 quoted elsewhere in this study is the same ratio at P4's realised control
tau of 52.1 min -- same Cooper-Helmstetter, same two coordinates, different tau.
The script reports the tau it used so the two are never confused. What it
asserts is that the ratio differs from 1, which holds at any tau.

Assert on `factor`, never on `prob`. `prob` is normalised over the whole genome,
and inserting the construct shifts every downstream coordinate by the insert
length -- about 7.2 kb, 0.15% of the genome -- so `prob` differs between
positions by ~0.1-0.2% for reasons that have nothing to do with the pin.

Usage
-----
	export PYTHONPATH="$(pwd):$PYTHONPATH"

	# One genome: shows pinned vs unpinned side by side.
	python runscripts/copy_number_paper/test_wt_pin.py \
		out/<batch>/kb/simData.cPickle

	# Two or more: adds the cross-position ratios, which are the real test.
	python runscripts/copy_number_paper/test_wt_pin.py \
		out/<p4_batch>/kb/simData.cPickle \
		out/<p1_batch>/kb/simData.cPickle

Pass a raw `kb/simData.cPickle` (the ParCa output). The variant is applied
in-process, so this works before any simulation has been launched.

Exit status is 0 if every check passes, 1 otherwise.
"""

import os
import pickle
import sys

import numpy as np

from models.ecoli.sim.variants import new_gene_burden_ladder as ladder
from models.ecoli.sim.variants.new_gene_internal_shift import (
	determine_new_gene_ids_and_indices)
from wholecell.utils import units

# Any expressing index works; `factor` does not depend on expression at all.
TEST_VARIANT_INDEX = 4

TOLERANCE = 2e-3  # the genome-shift residue, see the module docstring


def _construct_factor(sim_data_path, pin):
	"""
	Apply the variant to a fresh copy of sim_data and return what the
	normaliser hands the construct.

	Returns a dict with the construct's real coordinate, its wildtype
	coordinate, and the `factor` synth_prob_from_ppgpp computes from it.
	"""
	with open(sim_data_path, 'rb') as f:
		sim_data = pickle.load(f)

	ladder.PIN_WT_COORDINATE = pin
	ladder.new_gene_burden_ladder(sim_data, TEST_VARIANT_INDEX)

	transcription = sim_data.process.transcription
	rna_data = transcription.rna_data
	_, new_gene_indices, _, _ = determine_new_gene_ids_and_indices(sim_data)

	# The ParCa's FITTED doubling time, not a realised one. The unpinned factor
	# ratio depends on it -- 0.753 here at 44 min against 0.787 at P4's
	# realised control tau of 52.1 min, same Cooper-Helmstetter at the same two
	# coordinates. Reported below so the number is never mistaken for a miss
	# against an in-simulation figure.
	doubling_time = sim_data.condition_to_doubling_time[sim_data.condition]
	ppgpp = sim_data.growth_rate_parameters.get_ppGpp_conc(doubling_time)
	_, factor = transcription.synth_prob_from_ppgpp(
		ppgpp, sim_data.process.replication.get_average_copy_number)

	# A native gene, to confirm the pin touched only the construct.
	native_index = next(
		i for i in range(len(rna_data['id'])) if i not in new_gene_indices)

	return dict(
		ids=[rna_data['id'][i] for i in new_gene_indices],
		tau=doubling_time.asNumber(units.min),
		coord=float(np.mean(
			[rna_data['replication_coordinate'][i] for i in new_gene_indices])),
		wt_coord=float(np.mean(
			[rna_data['wt_replication_coordinate'][i]
				for i in new_gene_indices])),
		factor=float(np.mean([factor[i] for i in new_gene_indices])),
		native_coord=int(rna_data['replication_coordinate'][native_index]),
		native_wt_coord=int(
			rna_data['wt_replication_coordinate'][native_index]),
		)


def _check(label, passed, detail):
	print('  [%s] %-52s %s' % ('PASS' if passed else 'FAIL', label, detail))
	return passed


def main(paths):
	ok = True
	measured = {}

	for path in paths:
		name = os.path.basename(os.path.dirname(os.path.dirname(path)))
		print('\n%s\n%s\n%s' % ('=' * 78, name or path, '=' * 78))

		pinned = _construct_factor(path, pin=True)
		unpinned = _construct_factor(path, pin=False)
		measured[name] = dict(pinned=pinned, unpinned=unpinned)

		print('  construct TUs: %s' % ', '.join(pinned['ids']))
		print('  ParCa fitted doubling time %9.1f min  '
			'(the factor ratio depends on this)' % pinned['tau'])
		print('  real coordinate            %12.0f' % pinned['coord'])
		print('  wt coordinate, unpinned    %12.0f' % unpinned['wt_coord'])
		print('  wt coordinate, pinned      %12.0f' % pinned['wt_coord'])
		print('  factor, unpinned           %12.6g' % unpinned['factor'])
		print('  factor, pinned             %12.6g' % pinned['factor'])
		print()

		ok &= _check(
			'real coordinate unchanged by the pin',
			pinned['coord'] == unpinned['coord'],
			'%.0f' % pinned['coord'])
		ok &= _check(
			'unpinned wt coordinate == real coordinate',
			unpinned['wt_coord'] == unpinned['coord'],
			'%.0f' % unpinned['wt_coord'])
		ok &= _check(
			'pinned wt coordinate == REFERENCE_WT_COORDINATE',
			pinned['wt_coord'] == ladder.REFERENCE_WT_COORDINATE,
			'%.0f' % pinned['wt_coord'])
		ok &= _check(
			'native gene untouched (wt == real)',
			pinned['native_wt_coord'] == pinned['native_coord'],
			'%d' % pinned['native_coord'])

		# At the reference locus the pin is a no-op by construction. Anywhere
		# else it must move the factor, or it did not take.
		at_reference = (
			round(pinned['coord']) == ladder.REFERENCE_WT_COORDINATE)
		if at_reference:
			ok &= _check(
				'reference locus: pin is a no-op',
				abs(pinned['factor'] / unpinned['factor'] - 1) < TOLERANCE,
				'factor ratio %.6f' % (pinned['factor'] / unpinned['factor']))
		else:
			ok &= _check(
				'off-reference locus: pin moves the factor',
				abs(pinned['factor'] / unpinned['factor'] - 1) > TOLERANCE,
				'factor ratio %.6f' % (pinned['factor'] / unpinned['factor']))

	if len(measured) < 2:
		print('\nOne genome only. The cross-position ratios -- the actual '
			'test -- need at least two.')
		return ok

	print('\n%s\nCross-position construct factor ratios\n%s' % ('=' * 78, '=' * 78))
	names = list(measured)
	base = names[0]
	print('  reference batch: %s\n' % base)
	print('  %-34s %14s %14s' % ('batch', 'unpinned', 'pinned'))

	for name in names[1:]:
		r_un = (measured[name]['unpinned']['factor']
			/ measured[base]['unpinned']['factor'])
		r_pin = (measured[name]['pinned']['factor']
			/ measured[base]['pinned']['factor'])
		print('  %-34s %14.6f %14.6f' % (name, r_un, r_pin))
		ok &= _check(
			'%s: pinned factor ratio == 1' % name,
			abs(r_pin - 1) < TOLERANCE,
			'%.6f' % r_pin)
		ok &= _check(
			'%s: unpinned factor ratio != 1' % name,
			abs(r_un - 1) > TOLERANCE,
			'%.6f  (this is the cancellation the pin removes)' % r_un)

	return ok


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print(__doc__)
		sys.exit(2)
	passed = main(sys.argv[1:])
	print('\nALL CHECKS PASSED' if passed else '\nCHECKS FAILED')
	sys.exit(0 if passed else 1)
