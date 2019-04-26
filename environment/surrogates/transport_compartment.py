from __future__ import absolute_import, division, print_function

import time
import numpy as np

from agent.inner import CellSimulation

TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [0, 128, 255]]

class TransportCompartment(CellSimulation):
	''''''

	def __init__(self, boot_config, synchronize_config, initialize_processes):

		# initialize transport
		self.processes = {}
		self.processes['transport'] = initialize_processes['transport'](boot_config, synchronize_config)

		# update ecoli's synchronize_config with transport fluxes.
		ecoli_synchronize_config = boot_config
		ecoli_synchronize_config['boundary_reactions'] = self.processes['transport'].transport_reactions_ids

		# initialize ecoli
		self.processes['ecoli'] = initialize_processes['ecoli'](boot_config, ecoli_synchronize_config)

		# resolve_inner_update defines what process each inner update is taken from
		self.resolve_inner_update = {
			'environment_change': 'ecoli',
			'volume': 'ecoli',
			'division': 'ecoli',
			'motile_force': 'transport',
			'transport_fluxes': 'transport'}

		# cross_update defines what inner_updates of one process become outer_update of another
		self.cross_updates = {
			'transport_fluxes': ('transport', 'ecoli')
			}
		self.added_update = {update_id: None for update_id in self.cross_updates.iterkeys()}

		# use one process' functions exclusively
		self.time = self.processes['transport'].time
		self.divide = self.processes['ecoli'].divide
		# self.divide = self.internal.divide  # TODO (eran) what happens when one process declares division? Can it divide the rest?


	def generate_inner_update(self):
		# sends environment a dictionary with relevant state changes

		# get updates from all of the individual processes
		process_updates = {}
		for process_id, process in self.processes.iteritems():
			process_update = process.generate_inner_update()
			process_updates[process_id] = process_update

		# add to inner update by pulling from individual process updates according to resolve_inner_update
		inner_update = {}
		for update_parameter, process_id in self.resolve_inner_update.iteritems():
			inner_update[update_parameter] = process_updates[process_id][update_parameter]

		# apply cross_updates: inner_updates from some processes become outer_updates used by other
		for update, (process_source, process_target) in self.cross_updates.iteritems():
			self.added_update[update] = process_updates[process_source][update]

		# inner updates from compartment
		inner_update['color'] = DEFAULT_COLOR

		return inner_update


	def apply_outer_update(self, outer_update):
		# add to update
		# TODO -- only apply outer update from ecoli to environment, outer update from transport is only for ecoli targets
		outer_update.update(self.added_update)
		for process_id, process in self.processes.iteritems():
			process.apply_outer_update(outer_update)


	def run_incremental(self, run_until):
		for process in self.processes.itervalues():
			process.run_incremental(run_until)

	def finalize(self):
		for process in self.processes.itervalues():
			process.finalize()
