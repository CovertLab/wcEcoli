from __future__ import absolute_import, division, print_function

import copy

from scipy import constants

from vivarium.core.process import Deriver
from vivarium.library.dict_utils import deep_merge
from vivarium.library.units import units
from vivarium_cell.processes.derive_globals import (
	length_from_volume,
	surface_area_from_length,
)


AVOGADRO = constants.N_A * 1 / units.mmol


class WcEcoliDeriveShape(Deriver):

	defaults = {
		'width': 1,  # um
		# Source: Wülfing, C., & Plückthun, A. (1994). Protein folding
		# in the periplasm of Escherichia coli. Molecular Microbiology,
		# 12(5), 685–692.
		# https://doi.org/10.1111/j.1365-2958.1994.tb01056.x
		'periplasm_volume_fraction': 0.3,
	}
	name = "wcEcoliDeriveShape"

	def __init__(self, initial_parameters=None):
		# type: (dict) -> None
		'''Derives cell length and surface area from width and volume.

		Ports:

		* **global**: Should be given the agent's boundary store.

		Arguments:
			initial_parameters (dict): A dictionary that can contain the
				follwing configuration options:

				* **width** (:py:class:`float`): Width of the cell in
				  microns
		'''
		if initial_parameters is None:
			initial_parameters = {}
		parameters = copy.deepcopy(self.defaults)
		deep_merge(parameters, initial_parameters)

		super(WcEcoliDeriveShape, self).__init__(parameters)

	def ports_schema(self):
		default_state = {
			'global': {
				'volume': 0 * units.fL,
				'width': self.parameters['width'],
				'length': 0,
				'surface_area': 0,
				'mmol_to_counts': 0,
				'periplasm_volume': 0,
			}
		}

		schema = {
			'global': {
				variable: {
					'_updater': 'set',
					'_emit': True,
					'_divider': (
						'set' if variable == 'width' else 'split'
					),
					'_default': default_state['global'][variable]
				}
				for variable in default_state['global']
			}
		}
		return schema

	def next_update(self, timestep, states):
		width = states['global']['width']
		volume = states['global']['volume']

		length = length_from_volume(volume.magnitude, width)
		surface_area = surface_area_from_length(length, width)
		mmol_to_counts = (AVOGADRO * volume).to('L/mmol')
		periplasm_volume = volume * self.parameters[
			'periplasm_volume_fraction']

		return {
			'global': {
				'length': length,
				'surface_area': surface_area,
				'mmol_to_counts': mmol_to_counts,
				'periplasm_volume': periplasm_volume,
			},
		}
