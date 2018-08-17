import time
import random

class EnvironmentNonSpatial(object):
	def __init__(self, concentrations):
		self._time = 0
		self.run_for = 1
		self.simulations = {}
		self.concentrations = concentrations

	def time(self):
		return self._time

	def add_simulation(self, id):
		state = {}
		self.simulations[id] = state

	def remove_simulation(self, id):
		return self.simulations.pop(id, {})

	def update_concentrations(self, all_changes):
		self._time += self.run_for
		for id, changes in all_changes.iteritems():
			for molecule, change in changes.iteritems():
				self.concentrations[molecule] += change

	def run_until(self):
		until = {}
		for id in self.simulations.keys():
			until[id] = self.time() + self.run_for

		return until

	def molecule_ids(self):
		return self.concentrations.keys()

	def get_concentrations(self):
		concentrations = {}
		for id in self.simulations.keys():
			concentrations[id] = self.concentrations

		return concentrations