from __future__ import absolute_import, division, print_function

import csv
import time
import uuid

import numpy as np

from agent.control import AgentControl, AgentCommand

# TODO: read these from CSV files instead of hard-coding them.
MEDIA = ['minimal', 'minimal_plus_amino_acids', 'minimal_minus_oxygen']

TRANSPORT_REACTIONS = ['reaction_1', 'reaction_2', 'reaction_3']

SUBSTRATES = ['GLC[p]', 'GLT[p]', 'CYS[p]', 'GLC[c]', 'GLT[c]', 'CYS[c]']

STOICHIOMETRY = [-1, 2, -1, 1, -1, 1]

INTERNAL_SUBSTRATE_COUNTS = {'GLC[c]': np.random.randint(1000, 100000),
								  'GLT[c]': np.random.randint(1000, 100000),
								  'CYS[c]': np.random.randint(1000, 100000)}
EXTERNAL_SUBSTRATE_COUNTS = {'GLC[p]': 0, 'GLT[p]': 0, 'CYS[p]': 0}

def substrate_counts():
	'''
	Concatenate a list of substrate names with random starting values for each substrate.
	Returns: dictionary of the form {substrate:int}

	'''
	counts = np.random.randint(100,1000, size=len(SUBSTRATES))
	return dict(zip(SUBSTRATES, counts))

def make_flux_distributions():
	'''
	Write flux distributions per media condition, per transport reaction. Used as test data.
	Returns: dictionary of the form {media:{reaction:[array of ints]}}

	'''
	flux_distributions = {}
	for condition in MEDIA:
		flux_distributions[condition] = {}
		for reaction in TRANSPORT_REACTIONS:
			flux_distributions[condition][reaction] = list(np.random.randint(10, size=10))
	return flux_distributions


def make_stoichiometry():
	'''
	Assign stoichiometry to each transport reaction to determine the magnitude and direction of substrates in each reaction
	# TODO: use a for loop to generalize stoichiometry assignment for != 3 TRANSPORT_REACTIONS
	Returns: dictionary of the form {reaction:{substrate:int}}

	'''
	stoichiometry = {
		TRANSPORT_REACTIONS[0]: {SUBSTRATES[0]: STOICHIOMETRY[0], SUBSTRATES[3]: STOICHIOMETRY[1]},
		TRANSPORT_REACTIONS[1]: {SUBSTRATES[1]: STOICHIOMETRY[2], SUBSTRATES[4]: STOICHIOMETRY[3]},
		TRANSPORT_REACTIONS[2]: {SUBSTRATES[2]: STOICHIOMETRY[4], SUBSTRATES[5]: STOICHIOMETRY[5]}}
	return stoichiometry

def write_data_file():
	'''
	Initialize a data file to write test data and visualize in a Jupyter notebook.
	Returns: CSV file 'test_data.csv' initialized with a list of strings as column headers

	'''
	with open('test_data.csv', 'wb') as test_data:
		filewriter = csv.writer(test_data, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		filewriter.writerow(['agent_id', 'timestep', 'volume', 'media', 'GLC[c]', 'GLT[c]', 'CYS[c]'])


class ShepherdControl(AgentControl):

	"""
	Send messages to the other agents in the system to trigger execution and/or shutdown
	the Outer agent (which sends messages to shutdown all the associated Inner agents) or
	shutdown specific Inner agents directly (which then report back to the Outer agent and
	then terminate).
	"""

	def __init__(self, agent_config):
		super(ShepherdControl, self).__init__(str(uuid.uuid1()), agent_config)

	def add_cell(self, agent_type, agent_config):
		# TODO(jerry): Bring back the --variant choice?
		self.add_agent(
			str(uuid.uuid1()),
			agent_type,
			agent_config)

	def lattice_experiment(self, args):
		write_data_file()
		lattice_id = str(uuid.uuid1())
		num_cells = args['number']
		print('Creating lattice agent_id {} and {} cell agents\n'.format(
			lattice_id, num_cells))
		self.add_agent(lattice_id, 'lattice', {})

		time.sleep(15)

		# TODO: change self.add_cell to transport.add_cell, removes need to maintain agent configuration from both control and transport
		for index in range(num_cells):
			self.add_cell(args['type'] or 'ecoli', {
				'outer_id': lattice_id,
				'working_dir': args['working_dir'],
				'seed': index,
				'flux': make_flux_distributions(),
				'stoichiometry': make_stoichiometry(),
				'substrate_counts': substrate_counts(),
				'internal_substrate_counts': INTERNAL_SUBSTRATE_COUNTS,
				'external_substrate_counts': EXTERNAL_SUBSTRATE_COUNTS
			})

	def chemotaxis_experiment(self, args):
		lattice_id = str(uuid.uuid1())
		num_cells = args['number']
		print('Creating lattice agent_id {} and {} cell agents\n'.format(
			lattice_id, num_cells))
		chemotaxis_config = {
			'run_for' : 1.0,
			'static_concentrations': True,
			'gradient': {'seed': True},
			'diffusion': 0.0,
			'translation_jitter': 0.0,
			'rotation_jitter': 0.0,
			'edge_length': 20.0,
			'patches_per_edge': 30}
		self.add_agent(lattice_id, 'lattice', chemotaxis_config)

		for index in range(num_cells):
			self.add_cell(args['type'] or 'chemotaxis', {
				'outer_id': lattice_id,
				'seed': index})

	# TODO: fix this, as it doesn't work correctly
	def transport_experiment(self, args):
		lattice_id = str(uuid.uuid1())
		num_cells = args['number']
		print('Creating lattice agent_id {} and {} cell agents\n'.format(
			lattice_id, num_cells))
		lattice_config = {
			'run_for' : 1.0,
			'static_concentrations': True,
			'gradient': {'seed': True},
			'diffusion': 0.0,
			'translation_jitter': 4.0,	# TODO: implement this movement
			'rotation_jitter': 1.0,		# TODO: implement this movement
			'edge_length': 20.0,
			'patches_per_edge': 30,}
		self.add_agent(lattice_id, 'lattice', lattice_config)

		for index in range(num_cells):
			self.add_cell(args['type'] or 'transport', {
				'outer_id': lattice_id,
				'seed': index,
				'flux': make_flux_distributions(),
				'stoichiometry': make_stoichiometry(),
				'substrate_counts': substrate_counts()})


class EnvironmentCommand(AgentCommand):
	"""
	Extend `AgentCommand` with new commands related to the lattice and ecoli experiments
	"""

	def __init__(self):
		choices = ['chemotaxis-experiment', 'transport-experiment']
		description = '''
Run an agent for the environmental context simulation.

The commands are:
`experiment [--number N] [--type T] [--working-dir D]` ask the Shepherd to run
    a lattice environment with N agents of type T,
`add --id OUTER_ID [--type T] [--config C]` ask the Shepherd to add an agent of
    type T with JSON configuration C to the environment OUTER_ID,
`remove --id AGENT_ID` ask all Shepherds to remove agent AGENT_ID,
`remove --prefix ID` ask all Shepherds to remove agents "ID*",
`trigger [--id OUTER_ID]` start running one or all simulations,
`pause [--id OUTER_ID]` pause one or all simulations,
`divide --id AGENT_ID` ask a cell agent to divide,
`shutdown [--id OUTER_ID]` shut down one or all environment agents and their
     connected agents,
'chemotaxis-experiment [--number N] [--type T]` ask the Shepherd to run a
    chemotaxis environment with N agents of type T
'transport-experiment [--number N] [--type T]` ask the Shepherd to run a 
	transport environment with N agents of type T'''

		super(EnvironmentCommand, self).__init__(choices, description)

	def experiment(self, args):
		self.require(args, 'number', 'working_dir')
		control = ShepherdControl({'kafka_config': self.kafka_config})
		control.lattice_experiment(args)
		control.shutdown()

	def chemotaxis_experiment(self, args):
		self.require(args, 'number')
		control = ShepherdControl({'kafka_config': self.kafka_config})
		control.chemotaxis_experiment(args)
		control.shutdown()

	def transport_experiment(self, args):
		self.require(args, 'number')
		control = ShepherdControl({'kafka_config': self.kafka_config})
		control.transport_experiment(args)
		control.shutdown()

if __name__ == '__main__':
	command = EnvironmentCommand()
	command.execute()
