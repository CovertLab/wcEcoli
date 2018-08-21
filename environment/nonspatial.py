import time
import random

import numpy as np
from scipy import constants

class EnvironmentNonSpatial(object):
	def __init__(self, concentrations):
		self._time = 0
		self._timestep = 0.2
		self.run_for = 20
		self.size = 1
		self.ndim = 3
		self.id = -1

		self.simulations = {}
		self.locations = {}
		self.concentrations = concentrations
		self.volume = 1000 #TODO (Eran) initialize this value
		self.nAvogadro = constants.N_A #TODO (ERAN) get this from sim_data.constants.nAvogadro


	def evolve(self):
		''' Evolve environment '''
		self.update_locations()


	def update_locations(self):
		''' Update location for all sim_ids '''
		for sim_id, location in self.locations.iteritems():
			location += np.random.normal(0, 0.1, self.ndim)


	def run_incremental(self, run_until):
		''' Simulate until run_until '''
		while self._time < run_until:
			self._time += self._timestep
			self.evolve()


	def counts_to_concentration(self, counts):
		''' Convert an array of counts to concentrations '''
		concentrations = [count / (self.volume * self.nAvogadro) for count in counts]
		return concentrations


	def update_counts(self, all_changes):
		'''
		Use delta counts from all the inner simulations, convert them to concentations,
		and add to the environmental concentrations
		'''

		for sim_id, delta_counts in all_changes.iteritems():
			delta_concentrations = self.counts_to_concentration(delta_counts.values())
			for molecule, delta_conc in zip(delta_counts.keys(), delta_concentrations):
				self.concentrations[molecule] += delta_conc

	def molecule_ids(self):
		''' Return the ids of all molecule species in teh environment '''
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

		# Place cell at a random initial location
		location = np.random.uniform(0,self.size,self.ndim)

		self.simulations[id] = state
		self.locations[id] = location


	def remove_simulation(self, id):
		return self.simulations.pop(id, {})
		return self.locations.pop(id, {})


	def simulations_run_until(self):
		until = {}
		run_until = self.time() + self.run_for
		for sim_id in self.simulations.keys():
			until[sim_id] = run_until

		# pass the environment a run_until
		until[self.id] = run_until

		return until
