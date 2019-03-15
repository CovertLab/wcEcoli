from __future__ import absolute_import, division, print_function

import time
import uuid
import numpy as np

from agent.control import AgentControl, AgentCommand

MEDIA = ['minimal', 'minimal_plus_amino_acids', 'minimal_minus_oxygen']

TRANSPORT_REACTIONS = ['reaction_1', 'reaction_2', 'reaction_3']

SUBSTRATES = ['GLC[e]', 'GLT[e]', 'CYS[e]', 'GLC[c]', 'GLT[c]', 'CYS[c]']

def substrate_counts():
	counts = np.random.randint(10, size=len(SUBSTRATES))
	return dict(zip(SUBSTRATES, counts))

def make_flux_distributions():
	flux_distributions = {}
	for condition in MEDIA:
		flux_distributions[condition] = {}
		for reaction in TRANSPORT_REACTIONS:
			flux_distributions[condition][reaction] = list(np.random.randint(10, size=10))
	return flux_distributions


def make_stoichiometry():
	stoichiometry = {
		TRANSPORT_REACTIONS[0]: {SUBSTRATES[0]: -1, SUBSTRATES[3]: 1},
		TRANSPORT_REACTIONS[1]: {SUBSTRATES[1]: -1, SUBSTRATES[4]: 1},
		TRANSPORT_REACTIONS[2]: {SUBSTRATES[2]: -1, SUBSTRATES[5]: 1}}
	return stoichiometry


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
		lattice_id = str(uuid.uuid1())
		num_cells = args['number']
		print('Creating lattice agent_id {} and {} cell agents\n'.format(
			lattice_id, num_cells))
		self.add_agent(lattice_id, 'lattice', {})

		time.sleep(10)

		for index in range(num_cells):
			self.add_cell(args['type'] or 'ecoli', {
				'outer_id': lattice_id,
				'working_dir': args['working_dir'],
				'seed': index,
				'flux': make_flux_distributions(),
				'stoichiometry': make_stoichiometry(),
				'substrate_counts': substrate_counts()})

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

	# TODO: make this work
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
			'translation_jitter': 1.0,
			'rotation_jitter': 0.2,
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
'transport-experiment [--number N] [--type T]` ask the Shepherd to run a transport environment with N agents of type T'''

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
