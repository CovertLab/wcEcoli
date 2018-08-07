import time
import random

class SimulationStub(object):
	def __init__(self):
		self.local_time = 0
		self._tracking_environment_change = False

	def time(self):
		return self.local_time

	def initialize_local_environment(self):
		self._tracking_environment_change = True
		self.local_set = False

	def set_local_environment(self, molecule_ids, concentrations):
		self.molecule_ids = molecule_ids
		self.concentrations = concentrations

		if not self.local_set:
			# self.concentrations = np.empty(len(molecule_ids))
			self.environment_change = {}
			for molecule in molecule_ids:
				self.environment_change[molecule] = 0
			self.local_set = True

	def run_incremental(self, run_until):
		time.sleep(1)
		for molecule in self.molecule_ids:
			self.environment_change[molecule] += random.randint(1, 6)
		self.local_time = run_until

	def accumulate_environment_change(self):
		deltas = self.external_states['Environment'].get_deltas()
		for molecule in self.molecule_ids:
			self.environment_change[molecule] += deltas[molecule]

	def get_environment_change(self):
		return self.environment_change

	def finalize(self):
		pass
