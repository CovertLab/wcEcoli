import time
import random

import os

import numpy as np
from scipy import constants

# Constants
N_AVOGADRO = constants.N_A #TODO (ERAN) get this from sim_data.constants.nAvogadro

VOLUME = 1E-10 #(L)

class EnvironmentBatchNonSpatial(object):
	def __init__(self, concentrations):
		self._time = 0
		self._timestep = 0.2
		self.run_for = 10

		self.id = -1

		self.simulations = {}

		self.molecule_ids = concentrations.keys()
		self.concentrations = concentrations.values()

		if os.path.exists("out/manual/batch_environment.txt"):
			os.remove("out/manual/batch_environment.txt")


	def run_incremental(self, run_until):
		''' Simulate until run_until '''
		self.save_environment()

		while self._time < run_until:
			self._time += self._timestep


	def save_environment(self):
		# open in append mode
		env_file = open("out/manual/batch_environment.txt", "a")
		env_file.write("%s\n" % self.concentrations)
		env_file.close()


	def counts_to_concentration(self, counts):
		''' Convert an array of counts to concentrations '''
		concentrations = [count / (VOLUME * N_AVOGADRO) for count in counts]
		return concentrations


	def update_counts(self, all_changes):
		'''
		Use delta counts from all the inner simulations, convert them to concentrations,
		and add to the environmental concentrations of each molecule at each simulation's location
		'''
		for sim_id, delta_counts in all_changes.iteritems():
			delta_concentrations = self.counts_to_concentration(delta_counts.values())

			for molecule, delta_conc in zip(delta_counts.keys(), delta_concentrations):
				self.concentrations[self.molecule_ids.index(molecule)] += delta_conc


	def get_molecule_ids(self):
		''' Return the ids of all molecule species in the environment '''
		return self.molecule_ids


	def get_concentrations(self):
		'''returns a dict with {molecule_id: conc} for each sim give its current location'''
		concentrations = {}
		for sim_id in self.simulations.keys():
			# get concentration
			concentrations[sim_id] = dict(zip(self.molecule_ids, self.concentrations))
		return concentrations


	def time(self):
		return self._time


	def add_simulation(self, id):
		state = {}
		self.simulations[id] = state


	def remove_simulation(self, id):
		self.simulations.pop(id, {})


	def simulations_run_until(self):
		until = {}
		run_until = self.time() + self.run_for
		for sim_id in self.simulations.keys():
			until[sim_id] = run_until

		# pass the environment a run_until
		until[self.id] = run_until

		return until
