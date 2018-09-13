from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

import fcl
import numpy as np
from scipy import constants

animating = 'ENVIRONMENT_ANIMATION' in os.environ

import matplotlib

# Turn off interactive plotting when running on sherlock
if animating:
	matplotlib.use('TKAgg')

import matplotlib.pyplot as plt

if animating:
	plt.ion()
	fig = plt.figure()

# Constants
N_AVOGADRO = constants.N_A  # TODO (ERAN) get this from sim_data.constants.nAvogadro
PI = np.pi

# Lattice parameters
N_DIMS = 2
PATCHES_PER_EDGE = 10
TOTAL_VOLUME = 1E-11  # (L)
EDGE_LENGTH = 10.  # (micrometers). for reference: e.coli length is on average 2 micrometers.

# Physical constants
DIFFUSION = 0.001  # diffusion constant. (micrometers^2/s) # TODO (Eran) calculate correct diffusion rate

# Derived environmental constants
PATCH_VOLUME = TOTAL_VOLUME / (PATCHES_PER_EDGE * PATCHES_PER_EDGE)
DX = EDGE_LENGTH / PATCHES_PER_EDGE  # intervals in x- directions (assume y- direction equivalent)
DX2 = DX * DX
# DT = DX2 * DX2 / (2 * DIFFUSION * (DX2 + DX2)) # upper limit on the time scale (go with at least 50% of this)

# Cell constants
CELL_RADIUS = 0.5  # (micrometers)
ORIENTATION_JITTER = PI / 40  # (radians/s)
LOCATION_JITTER = 0.01  # (micrometers/s)


class EnvironmentSpatialLattice(object):
	def __init__(self, concentrations):
		'''FCL Workflow:
		1. Initialize manager to handle collisions between agents.
		2. Initialize a dictionary to hold onto collision geometry.
		'''
		self._time = 0
		self._timestep = 1.0
		self.run_for = 5
		self.agent_id = -1
		self.simulations = {}
		self.locations = {}

		# Initialize a "manager" to handle collisions between "box" CollisionObjects.
		self.manager = fcl.DynamicAABBTreeCollisionManager()
		# Initialize a "boxes" dictionary with "agent_id" keys and "box" values.
		self.boxes = {}
		# Initialize a "collisions" list with "collision anywhere?" boolean and "minimum distance" float
		self.collisions = []

		self._molecule_ids = concentrations.keys()
		self.concentrations = concentrations.values()
		self.molecule_index = {molecule: index for index, molecule in enumerate(self._molecule_ids)}

		# Create lattice and fill each site with concentrations dictionary
		# Molecule identities are defined along the major axis, with spatial dimensions along the other two axes.
		self.lattice = np.empty([len(self._molecule_ids)] + [PATCHES_PER_EDGE for dim in xrange(N_DIMS)],
								dtype=np.float64)
		for idx, molecule in enumerate(self._molecule_ids):
			self.lattice[idx].fill(self.concentrations[idx])

		if os.path.exists("out/manual/environment.txt"):
			os.remove("out/manual/environment.txt")
		if os.path.exists("out/manual/locations.txt"):
			os.remove("out/manual/locations.txt")

		glucose_lattice = self.lattice[self.molecule_index['GLC[p]']]
		plt.imshow(glucose_lattice, vmin=0, vmax=25, cmap='YlGn')
		plt.colorbar()
		plt.axis('off')
		plt.pause(0.0001)

	def evolve(self):
		''' Evolve environment '''

		# TODO (Eran) updating location with the environment causes cells to jump at run_for ts, and then cells deposit
		# all of their deltas at the jumps, rather than along their path laid out during evolve()
		self.update_locations()
		self.run_diffusion()

	def update_locations(self):
		''' Update location for all agent_ids '''
		for agent_id, location in self.locations.iteritems():
			# Move the cell around randomly
			self.locations[agent_id][0:2] = (location[0:2] +
											 np.random.normal(scale=np.sqrt(LOCATION_JITTER * self._timestep),
															  size=N_DIMS)
											 ) % EDGE_LENGTH

			# Orientation jitter
			self.locations[agent_id][2] = (location[2] +
										   np.random.normal(scale=ORIENTATION_JITTER * self._timestep)
										   )

			# Call the CollisionObject "box" from the boxes dictionary
			box = self.boxes[agent_id]

			# Mirror the random translation to box
			box.setTranslation(np.array[self.locations[agent_id][0], self.locations[agent_id][1], 0.0])

			# Mirror the orientation jitter to box CollisionObject
			orientation = (self.locations[agent_id][2])

			# See add_simulation for explanation of this transform
			q_z = np.sin(orientation / 2) * 1.0
			q_w = np.cos(orientation / 2)
			box.setRotation(np.array[0.0, 0.0, q_z, q_w])

			# Undergo collision detection
			self.collision_detection()

	### FCL Requirement ###
	def collision_detection(self):
		'''FCL Workflow:
		1. Setup the collision manager to contain all boxes.
		2. Check for any box-to-box collisions within manager.
		3. If Step 2 = True, check for closest box-to-box distance within manager.
		4. Write collision & distance results to dictionary for visualization.
		'''
		# TODO: Add real volume exclusion using overlapping distance vectors
		# TODO: Feed collisions dictionary to visualization layer

		# Setup CollisionObject manager
		self.manager.setup()

		# Setup data structure to hold collision data called collision_data
		collision_data = fcl.CollisionData()
		# Check box-to-box collisions within manager
		self.manager.collide(collision_data, fcl.defaultCollisionCallback)

		# If there are ANY collisions within manager, then result.is_collision returns True
		if collision_data.result.is_collision == True:
			# Setup data structure to hold distance data called distance_data
			distance_data = fcl.DistanceData()
			# Check box-to-box distances within manager
			self.manager.distance(distance_data, fcl.defaultDistanceCallback)
			# If cells collide, try randomly moving them again
			self.update_locations()  # FIXME: with many cells this may effectively call an infinite loop!

		# Store collision and distance data in "collisions" list
		self.collisions = [collision_data.result.is_collision, distance_data.result.min_distance]

	def run_diffusion(self):
		change_lattice = np.zeros(self.lattice.shape)
		for idx in xrange(len(self.lattice)):
			molecule = self.lattice[idx]

			# run diffusion if molecule field is not uniform
			if (len(set(molecule.flatten())) != 1):
				change_lattice[idx] = self.diffusion_timestep(molecule)

		self.lattice += change_lattice

	def diffusion_timestep(self, lattice):
		''' calculate concentration changes cause by diffusion. Assumes periodic lattice, with wrapping'''

		# TODO (Eran) write this as matrix operation rather than np.roll.
		N = np.roll(lattice, 1, axis=0)
		S = np.roll(lattice, -1, axis=0)
		W = np.roll(lattice, 1, axis=1)
		E = np.roll(lattice, -1, axis=1)

		change_lattice = DIFFUSION * self._timestep * ((N + S + W + E - 4 * lattice) / DX2)

		return change_lattice

	def run_incremental(self, run_until):
		''' Simulate until run_until '''
		self.output_environment()
		self.output_locations()

		while self._time < run_until:
			self._time += self._timestep
			self.evolve()

	def output_environment(self):
		'''plot environment lattice'''
		glucose_lattice = self.lattice[self.molecule_index['GLC[p]']]

		plt.clf()
		plt.imshow(glucose_lattice, cmap='YlGn')
		plt.colorbar()
		plt.axis('off')

	def output_locations(self):
		'''plot cell locations and orientations'''
		for agent_id, location in self.locations.iteritems():
			y = location[0] * PATCHES_PER_EDGE / EDGE_LENGTH
			x = location[1] * PATCHES_PER_EDGE / EDGE_LENGTH
			theta = location[2]
			volume = self.simulations[agent_id]['volume']

			# get length, scaled to lattice resolution
			length = self.volume_to_length(volume)

			dx = length * PATCHES_PER_EDGE / EDGE_LENGTH * np.sin(theta)
			dy = length * PATCHES_PER_EDGE / EDGE_LENGTH * np.cos(theta)

			plt.plot([x - dx / 2, x + dx / 2], [y - dy / 2, y + dy / 2],
					 color='salmon', linewidth=CELL_RADIUS / EDGE_LENGTH * 600, solid_capstyle='round')
			plt.plot([x - dx * 9.5 / 20, x + dx * 9.5 / 20], [y - dy * 9.5 / 20, y + dy * 9.5 / 20],
					 color='slateblue', linewidth=CELL_RADIUS / EDGE_LENGTH * 450, solid_capstyle='round')

		if animating:
			plt.pause(0.0001)

	def volume_to_length(self, volume):
		'''
		get cell length from volume, using the following equation for capsule volume, with V=volume, r=radius,
		a=length of cylinder without rounded caps, l=total length:

		V = (4/3)*PI*r^3 + PI*r^2*a
		l = a + 2*r
		'''

		cylinder_length = (volume - (4 / 3) * PI * CELL_RADIUS ** 3) / (PI * CELL_RADIUS ** 2)
		total_length = cylinder_length + 2 * CELL_RADIUS

		return total_length

	def count_to_concentration(self, count):
		''' Convert count to concentrations '''
		return count / (PATCH_VOLUME * N_AVOGADRO)

	def update_from_simulations(self, all_changes):
		'''
		Use change counts from all the inner simulations, convert them to concentrations,
		and add to the environmental concentrations of each molecule at each simulation's location
		'''
		self.simulations = all_changes
		for agent_id, state in self.simulations.iteritems():
			location = self.locations[agent_id][0:2] * PATCHES_PER_EDGE / EDGE_LENGTH
			patch_site = tuple(np.floor(location).astype(int))

			for molecule, count in state['environment_change'].iteritems():
				concentration = self.count_to_concentration(count)
				index = self.molecule_index[molecule]
				self.lattice[index, patch_site[0], patch_site[1]] += concentration

	def get_molecule_ids(self):
		''' Return the ids of all molecule species in the environment '''
		return self._molecule_ids

	def get_concentrations(self):
		'''returns a dict with {molecule_id: conc} for each sim give its current location'''
		concentrations = {}
		for agent_id in self.simulations.keys():
			# get concentration from cell's given bin
			location = self.locations[agent_id][0:2] * PATCHES_PER_EDGE / EDGE_LENGTH
			patch_site = tuple(np.floor(location).astype(int))

			concentrations[agent_id] = dict(zip(self._molecule_ids, self.lattice[:, patch_site[0], patch_site[1]]))
		return concentrations

	def time(self):
		return self._time

	def add_simulation(self, agent_id, state):
		'''FCL Workflow:
		1. Initialize CollisionObject of type "box" as a bounding box around an agent.
		2. Give box some geometry, a translation, and a rotation to match agent.
		3. Write box to a dictionary value with agent_id as its key.
		'''
		# Place cell at a random initial location
		location = np.random.uniform(0, EDGE_LENGTH, N_DIMS)
		orientation = np.random.uniform(0, 2 * PI)

		# Initialize a CollisionObject with geometry "box" and (x,y,z) = (radius, length, 0.0)
		box = fcl.CollisionObject(fcl.Box(CELL_RADIUS, self.volume_to_length(state['volume']), 0.0))

		# Mirror location onto box CollisionObject using setTranslation()
		# setTranslation reads np.array(x, y, z) = (location(x), location(y), 0.0)
		box.setTranslation(np.array([location.item(0), location.item(1), 0.0]))

		# Mirror orientation onto box CollisionObject using setRotation()
		# setRotation() reads a 4D quaternion np.array (x, y, z, w)

		''' To convert 2D orientation theta to 4D, use axis-angle around x:
		q_x = sin(orientation/2) * axis_x
		q_y = sin(orientation/2) * axis_y
		q_z = sin(orientation/2) * axis_z
		q_w = cos(orientation/2)
		'''

		# In 2D axis.x & axis.y = 0.0 & axis.z is 1.0 since we rotate around z
		q_z = np.sin(orientation / 2) * 1.0
		q_w = np.cos(orientation / 2)
		box.setRotation(np.array[0.0, 0.0, q_z, q_w])

		# Write dictionary "boxes" with "agent_id" as keys and "box" as values
		self.boxes[agent_id] = box
		self.simulations[agent_id] = state
		self.locations[agent_id] = np.hstack((location, orientation))

	def remove_simulation(self, agent_id):
		self.simulations.pop(agent_id, {})
		self.locations.pop(agent_id, {})
		self.boxes.pop(agent_id, {})

	def run_simulations_until(self):
		until = {}
		run_until = self.time() + self.run_for
		for agent_id in self.simulations.keys():
			until[agent_id] = run_until

		# Pass the environment a run_until
		until[self.agent_id] = run_until

		return until
