import time
import random

class SimulationStub(object):
	def __init__(self):
		self.local_time = 0
		self._tracking_environment_change = False

	def initialize_local_environment(self):
		self._tracking_environment_change = True

	def time(self):
		return self.local_time

	def set_local_environment(self, molecule_ids, concentrations):
		self.molecule_ids = molecule_ids
		self.molecule_changes = {}
		for molecule in self.molecule_ids:
			if not molecule in self.molecule_changes:
				self.molecule_changes[molecule] = 0
		self.concentrations = concentrations

	def run_incremental(self, run_until):
		time.sleep(1)
		for molecule in self.molecule_ids:
			self.molecule_changes[molecule] += random.randint(1, 6)
		self.local_time = run_until

	def get_environment_change(self):
		return self.molecule_changes

	def finalize(self):
		pass
