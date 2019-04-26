from __future__ import absolute_import, division, print_function

import time
import numpy as np

from agent.inner import CellSimulation

TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [100, 100, 100]]

class TransportCompartment(CellSimulation):
	''''''

	def __init__(self, boot_config, synchronize_config, initialize_processes):

		# initialize transport
		self.processes = {}
		self.processes['transport'] = initialize_processes['transport'](boot_config, synchronize_config)

		# update ecoli's synchronize_config with transport fluxes.
		ecoli_synchronize_config = boot_config
		ecoli_synchronize_config['target_transport_reactions'] = self.processes['transport'].transport_reactions_ids

		# initialize ecoli
		self.processes['ecoli'] = initialize_processes['ecoli'](boot_config, ecoli_synchronize_config)



		# TODO -- use compartment's color as default


		# # resolve_inner_update defines what process each inner update is taken from
		# self.resolve_inner_update = {
		# 	'environment_change': 'internal',
		# 	'volume': 'internal',
		# 	'division': 'internal',
		# 	'color': 'transport',
		# 	'motile_force': 'transport',
		# 	'transport_fluxes': 'transport'}
		#
		# # cross_update defines what inner_updates of one process become outer_update of another
		# self.cross_updates = {
		# 	'transport_fluxes': ('transport', 'internal')}
		# self.added_update = {update_id: None for update_id in self.cross_updates.iterkeys()}
		#
		# # use one process' functions exclusively


		# self.run_incremental = self.processes['transport'].run_incremental
		self.generate_inner_update = self.processes['transport'].generate_inner_update
		self.apply_outer_update = self.processes['transport'].apply_outer_update
		self.time = self.processes['transport'].time
		self.divide = self.processes['transport'].divide
		# self.divide = self.internal.divide  # TODO (eran) what happens when one process declares division? Can it divide the rest?


	# def generate_inner_update(self):
	# 	# sends environment a dictionary with relevant state changes
	# 	process_updates = {}
	# 	for process_id, process in self.processes.iteritems():
	# 		process_update = process.generate_inner_update()
	# 		process_updates[process_id] = process_update
	#
	# 	# TODO -- cross updates between the processes and pass an aggregated "inner_update"
	# 	# inner_update = {}
	# 	# for update, process_id in self.resolve_inner_update.iteritems():
	# 	# 	inner_update[update] = process_updates[process_id][update]
	# 	#
	# 	# # TODO -- get inner updates that will be passed to outer_update
	# 	# for update, (process_source, process_target) in self.cross_updates.iteritems():
	# 	# 	self.added_update[update] = process_updates[process_source][update]
	#
	#
	# 	inner_update = process_updates['transport']
	# 	inner_update['color'] = DEFAULT_COLOR
	# 	print('inner_update: ' + str(inner_update))
	# 	# import ipdb; ipdb.set_trace()
	#
	# 	return inner_update
	#
	# def apply_outer_update(self, outer_update):
	# 	for process_id, process in self.processes.iteritems():
	# 		process.apply_outer_update(outer_update)
	#
	# 	print('outer_update: ' + str(outer_update))
	# 	# import ipdb; ipdb.set_trace()
	#
	# 	# # add added_update from processes
	# 	# outer_update.update(self.added_update)
	# 	# for process in self.processes.itervalues():
	# 	# 	process.apply_outer_update(outer_update)
	# 	#
	# 	# # TODO -- only apply outer update from ecoli to environment, outer update from transport is only for ecoli targets
	#
	def run_incremental(self, run_until):
		for process in self.processes.itervalues():
			process.run_incremental(run_until)

	def finalize(self):
		for process in self.processes.itervalues():
			process.finalize()
