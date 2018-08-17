import time
import random

from wholecell.utils import units

from scipy import constants

class EnvironmentNonSpatial(object):
	def __init__(self, concentrations):
		self._time = 0
		self.run_for = 20
		self.simulations = {}
		self.concentrations = concentrations
		self.volume = 1 #TODO (Eran) initialize this value
		self.nAvogadro = constants.N_A #TODO (ERAN) get this from sim_data.constants.nAvogadro


	def counts_to_concentration(self, counts):
		concentrations = [count / (self.volume * self.nAvogadro) for count in counts]
		return concentrations


	def update_counts(self, all_changes):
		self._time += self.run_for
		for sim_id, delta_counts in all_changes.iteritems():
			delta_concentrations = self.counts_to_concentration(delta_counts.values())
			for molecule, delta_conc in zip(delta_counts.keys(), delta_concentrations):
				self.concentrations[molecule] += delta_conc


	def molecule_ids(self):
		return self.concentrations.keys()


	def get_concentrations(self):
		concentrations = {}
		for sim_id in self.simulations.keys():
			concentrations[sim_id] = self.concentrations
		return concentrations


	def time(self):
		return self._time

	def add_simulation(self, id):
		state = {}
		self.simulations[id] = state

	def remove_simulation(self, id):
		return self.simulations.pop(id, {})

	def run_until(self):
		until = {}
		for sim_id in self.simulations.keys():
			until[sim_id] = self.time() + self.run_for

		return until

