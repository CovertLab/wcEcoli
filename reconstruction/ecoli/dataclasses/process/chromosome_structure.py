"""
SimulationData for the ChromosomeStructure process
"""

from __future__ import division, print_function, absolute_import


class ChromosomeStructure(object):
	"""
	SimulationData for the ChromosomeStructure process
	"""

	def __init__(self, raw_data, sim_data):
		self._build_supercoiling_parameters(raw_data, sim_data)


	def _build_supercoiling_parameters(self, raw_data, sim_data):
		"""
		Load parameters used for DNA supercoiling from raw_data.
		"""
		for name, value in raw_data.dna_supercoiling.items():
			self.__setattr__(name, value)

		# Calculate rate constants of gyrase binding/unbinding
		self.gyrase_unbinding_rate_constant = 1 / self.mean_gyrase_dwell_time
		self.gyrase_binding_rate_constant = self.gyrase_unbinding_rate_constant * (
			self.fraction_bound_gyrase / (1 - self.fraction_bound_gyrase)
		)

		# Unique index used for dummy molecule at terC
		self.terC_dummy_molecule_index = -1
