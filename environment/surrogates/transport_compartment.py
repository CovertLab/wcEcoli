from __future__ import absolute_import, division, print_function

from agent.inner import CellSimulation

TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [0, 128, 255]]


def merge_two_dicts(x, y):
	z = x.copy()  # start with x's keys and values
	z.update(y)  # modifies z with y's keys and values & returns None
	return z

class TransportCompartment(CellSimulation):
	'''
	TODO -- explain the cross wiring configurations that the compartment can take on
	'''

	def __init__(self, boot_config, synchronize_config, network_config):
		initialize_processes = network_config['initilize']
		self.connections = network_config['connections']
		# TODO -- don't hardcode cross_update
		self.cross_update = {
			'transport': {},
			'ecoli': {}}

		# initialize transport
		self.processes = {}
		self.processes['transport'] = initialize_processes['transport'](boot_config, synchronize_config)

		# update ecoli's synchronize_config using transport_reactions_ids from transport
		ecoli_synchronize_config = boot_config
		ecoli_synchronize_config['boundary_reactions'] = self.processes['transport'].transport_reactions_ids

		# initialize ecoli
		self.processes['ecoli'] = initialize_processes['ecoli'](boot_config, ecoli_synchronize_config)

		# use one process' functions exclusively
		self.time = self.processes['transport'].time
		self.divide = self.processes['ecoli'].divide
		# TODO (eran) what happens when one process declares division? Can it divide the rest? How do they share inheritance, etc.


	def generate_inner_update(self):
		'''
		Sends environment a dictionary with relevant state changes.
		Creates dictionary with messages between processes.
		'''

		# get updates from all of the individual processes
		process_updates = {}
		for process_id, process in self.processes.iteritems():
			process_update = process.generate_inner_update()
			process_updates[process_id] = process_update

		# add to inner update by pulling from individual process updates according to resolve_inner_update
		inner_update = {}
		for source, target in self.connections.iteritems():
			source_process, source_message = source.split('.')
			target_process, target_message = target.split('.')

			if target_process == 'compartment':
				inner_update[target_message] = process_updates[source_process][source_message]
			else:
				self.cross_update[target_process][target_message] = process_updates[source_process][source_message]

		# inner updates from compartment
		inner_update['color'] = DEFAULT_COLOR

		return inner_update


	def apply_outer_update(self, outer_update):
		for process_id, process in self.processes.iteritems():
			cross_update = self.cross_update[process_id]
			process_update = merge_two_dicts(outer_update, cross_update)
			process.apply_outer_update(process_update)

	def run_incremental(self, run_until):
		for process in self.processes.itervalues():
			process.run_incremental(run_until)

	def finalize(self):
		for process in self.processes.itervalues():
			process.finalize()
