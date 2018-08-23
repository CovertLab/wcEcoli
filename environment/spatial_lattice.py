import time
import random

import os

import numpy as np
from scipy import constants

# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt

class EnvironmentSpatialLattice(object):
	def __init__(self, concentrations):
		self._time = 0
		self._timestep = 0.2
		self.run_for = 20
		self.size = 1
		self.ndims = 2
		self.nbins = 20
		self.volume = 1000 #TODO (Eran) initialize this value
		self.nAvogadro = constants.N_A #TODO (ERAN) get this from sim_data.constants.nAvogadro

		self.id = -1

		self.simulations = {}
		self.locations = {}

		self.molecules = concentrations.keys()

		# Create lattice and fill each site with concentrations dictionary
		#TODO (Eran) make this an ndims+1, with +1 for indexing each molecule
		self.lattice = np.empty([self.nbins for dim in xrange(self.ndims)], dtype=dict)
		self.lattice.fill(concentrations)

		if os.path.exists("out/manual/environment.txt"):
			os.remove("out/manual/environment.txt")
		if os.path.exists("out/manual/locations.txt"):
			os.remove("out/manual/locations.txt")


	def save_environment(self):

		glucose_lattice = np.zeros_like(self.lattice, dtype=float)
		for (x, y), value in np.ndenumerate(self.lattice):
			glucose_lattice[x][y] = self.lattice[x][y]['GLC[p]']

		glucose_lattice = glucose_lattice.tolist()

		# open in append mode
		lattice_file = open("out/manual/environment.txt", "a")
		lattice_file.write("%s\n" % glucose_lattice)
		lattice_file.close()


	def save_locations(self):

		locations = self.locations.values()[0].tolist()

		# open in append mode
		locations_file = open("out/manual/locations.txt", "a")
		locations_file.write("%s\n" % locations)
		locations_file.close()


	def evolve(self):
		''' Evolve environment '''

		self.update_locations()

		# diffuse environmental concentrations
		self.diffusion()


	def update_locations(self):
		''' Update location for all sim_ids '''
		for sim_id, location in self.locations.iteritems():
			location += np.random.normal(0, 0.1, self.ndims)

			# lattice cutoff
			location[location < 0] = 0
			location[location >= self.size] = self.size - 0.000001 # minus infinitesimal helps keep within lattice



	def diffusion(self):
		pass


	def run_incremental(self, run_until):
		''' Simulate until run_until '''

		self.save_environment()
		self.save_locations()

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
			location = self.locations[sim_id] * self.nbins
			bin = tuple(np.floor(location).astype(int))

			delta_concentrations = self.counts_to_concentration(delta_counts.values())
			for molecule, delta_conc in zip(delta_counts.keys(), delta_concentrations):
				self.lattice[bin][molecule] += delta_conc


	def molecule_ids(self):
		''' Return the ids of all molecule species in the environment '''
		return self.molecules


	def get_concentrations(self):
		concentrations = {}
		for sim_id in self.simulations.keys():
			# get concentration from cell's given bin
			location = self.locations[sim_id] * self.nbins
			bin = tuple(np.floor(location).astype(int))
			concentrations[sim_id] = self.lattice[bin]
		return concentrations


	def time(self):
		return self._time


	def add_simulation(self, id):
		state = {}

		# Place cell at a random initial location
		location = np.random.uniform(0,self.size,self.ndims)

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
