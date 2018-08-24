from __future__ import absolute_import, division, print_function

import json
import uuid
import argparse

import agent.event as event
from agent.agent import Agent
from agent.outer import Outer
from agent.inner import Inner
from agent.stub import SimulationStub, EnvironmentStub

from environment.nonspatial import EnvironmentNonSpatial
from environment.two_dim_lattice import EnvironmentSpatialLattice

from models.ecoli.sim.simulation import EcoliSimulation

from wholecell.fireworks.firetasks import VariantSimDataTask

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

from wholecell.utils import units

default_kafka_config = {
	'host': '127.0.0.1:9092',
	'simulation_send': 'environment_listen',
	'simulation_receive': 'environment_broadcast',
	'environment_control': 'environment_control',
	'subscribe_topics': []}

class BootOuter(object):

	"""
	Initialize the `EnvironmentStub`, pass it to the `Outer` agent and launch the process.

	This is a demonstration of how to initialize an Outer component. In the place of 
	`EnvironmentStub` you would substitute your own environment class that meets the interface
	defined in `Outer`. 
	"""

	def __init__(self, kafka_config):
		volume = 1
		concentrations = {
			'yellow': 5,
			'green': 11,
			'red': 44,
			'blue': 12}

		self.environment = EnvironmentStub(volume, concentrations)
		self.outer = Outer(kafka_config, self.environment)

class BootInner(object):

	"""
	Initialize the `SimulationStub`, pass it to the `Inner` agent and launch the process.

	This is a demonstration of how to initialize an Inner component. When creating your 
	own simulation you would supply a class that meets the same interface as the `SimulationStub`
	that would be driven by the Inner agent in response to messages from its corresponding 
	Outer agent.
	"""

	def __init__(self, id, kafka_config):
		self.id = id
		self.simulation = SimulationStub()
		self.inner = Inner(
			kafka_config,
			self.id,
			self.simulation)


class BootEnvironmentNonSpatial(object):
	def __init__(self, kafka_config):
		raw_data = KnowledgeBaseEcoli()
		# create a dictionary with all saved environments
		self.environment_dict = {}
		for label in dir(raw_data.condition.environment):
			if label.startswith("__"):
				continue
			self.environment_dict[label] = {}
			# initiate all molecules with 0 concentrations
			for row in raw_data.condition.environment_molecules:
				self.environment_dict[label].update({row["molecule id"]: 0}) #* (units.mmol / units.L)})
			# update non-zero concentrations
			molecule_concentrations = getattr(raw_data.condition.environment, label)
			for row in molecule_concentrations:
				self.environment_dict[label].update({row["molecule id"]: row["concentration"].asNumber()})

		# TODO (Eran) don't hardcode initial environment, get this from timeseries
		concentrations = self.environment_dict['minimal']

		self.environment = EnvironmentNonSpatial(concentrations)
		self.outer = Outer(kafka_config, self.environment)


class BootEnvironmentSpatialLattice(object):
	def __init__(self, kafka_config):
		raw_data = KnowledgeBaseEcoli()
		# create a dictionary with all saved environments
		self.environment_dict = {}
		for label in dir(raw_data.condition.environment):
			if label.startswith("__"):
				continue
			self.environment_dict[label] = {}
			# initiate all molecules with 0 concentrations
			for row in raw_data.condition.environment_molecules:
				self.environment_dict[label].update({row["molecule id"]: 0}) #* (units.mmol / units.L)})
			# update non-zero concentrations
			molecule_concentrations = getattr(raw_data.condition.environment, label)
			for row in molecule_concentrations:
				self.environment_dict[label].update({row["molecule id"]: row["concentration"].asNumber()})

		# TODO (Eran) don't hardcode initial environment, get this from timeseries
		concentrations = self.environment_dict['minimal']

		self.environment = EnvironmentSpatialLattice(concentrations)
		self.outer = Outer(kafka_config, self.environment)


class BootEcoli(object):
	def __init__(self, id, kafka_config):
		self.id = id

		sim_data_fit = '/Users/eranagmon/code/wcEcoli/out/manual/kb/simData_Most_Fit.cPickle'
		sim_data_variant = '/Users/eranagmon/code/wcEcoli/out/manual/wildtype_000000/kb/simData_Modified.cPickle'

		# copy the file simData_Most_Fit.cPickle to simData_Modified.cPickle
		task = VariantSimDataTask(
			variant_function='wildtype',
			variant_index=0,
			input_sim_data=sim_data_fit,
			output_sim_data=sim_data_variant,
			variant_metadata_directory='/Users/eranagmon/code/wcEcoli/out/manual/wildtype_000000/metadata',
		)
		task.run_task({})

		with open(sim_data_variant, "rb") as input_sim_data:
			sim_data = cPickle.load(input_sim_data)

		options = {}
		options["simData"] = sim_data
		options["outputDir"] = '/Users/eranagmon/code/wcEcoli/out/manual/sim_' + self.id + '/simOut'
		options["logToDisk"] = True
		options["overwriteExistingFiles"] = True

		options["seed"] = 0
		options["lengthSec"] = 10800
		options["timeStepSafetyFraction"] = 1.3
		options["maxTimeStep"] = 0.9
		options["updateTimeStepFreq"] = 5
		options["logToShell"] = True
		options["logToDiskEvery"] = 1
		options["massDistribution"] = True
		options["growthRateNoise"] = False
		options["dPeriodDivision"] = False
		options["translationSupply"] = True

		self.simulation = EcoliSimulation(**options)
		self.inner = Inner(
			kafka_config,
			self.id,
			self.simulation)


class EnvironmentControl(Agent):

	"""
	Send messages to the other agents in the system to trigger execution and/or shutdown
	the Outer agent (which sends messages to shutdown all the associated Inner agents) or
	shutdown specific Inner agents directly (which then report back to the Outer agent and
	then terminate).
	"""

	def __init__(self, kafka_config=default_kafka_config):
		id = 'environment_control'
		super(EnvironmentControl, self).__init__(id, kafka_config)

	def trigger_execution(self):
		self.send(self.kafka_config['environment_control'], {
			'event': event.TRIGGER_EXECUTION})

	def shutdown_environment(self):
		self.send(self.kafka_config['environment_control'], {
			'event': event.SHUTDOWN_ENVIRONMENT})

	def shutdown_simulation(self, id):
		self.send(self.kafka_config['simulation_receive'], {
			'event': event.SHUTDOWN_SIMULATION,
			'inner_id': id})

def main():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	parser = argparse.ArgumentParser(
		description='Boot the various agents for the environmental context simulation')

	parser.add_argument(
		'command',
		choices=['inner', 'outer', 'ecoli', 'nonspatial', 'lattice', 'trigger', 'shutdown'],
		help='which command to boot')

	parser.add_argument(
		'--id',
		help='unique identifier for simulation agent')

	parser.add_argument(
		'--kafka-host',
		default='127.0.0.1:9092',
		help='address for Kafka server')

	parser.add_argument(
		'--environment-control',
		default='environment_control',
		help='topic the environment will receive control messages on')

	parser.add_argument(
		'--simulation-receive',
		default='environment_broadcast',
		help='topic the simulations will receive messages on')

	parser.add_argument(
		'--simulation-send',
		default='environment_listen',
		help='topic the simulations will send messages on')

	args = parser.parse_args()
	kafka_config = {
		'host': args.kafka_host,
		'environment_control': args.environment_control,
		'simulation_receive': args.simulation_receive,
		'simulation_send': args.simulation_send,
		'subscribe_topics': []}

	if args.command == 'inner':
		if not args.id:
			raise ValueError('--id must be supplied for inner command')

		inner = BootInner(args.id, kafka_config)

	elif args.command == 'outer':
		outer = BootOuter(kafka_config)

	if args.command == 'ecoli':
		if not args.id:
			raise ValueError('--id must be supplied for inner command')

		inner = BootEcoli(args.id, kafka_config)

	elif args.command == 'nonspatial':
		outer = BootEnvironmentNonSpatial(kafka_config)

	elif args.command == 'lattice':
		outer = BootEnvironmentSpatialLattice(kafka_config)

	elif args.command == 'trigger':
		control = EnvironmentControl(kafka_config)
		control.trigger_execution()
		control.shutdown()

	elif args.command == 'shutdown':
		control = EnvironmentControl(kafka_config)

		if not args.id:
			control.shutdown_environment()
		else:
			control.shutdown_simulation(args.id)
		control.shutdown()

if __name__ == '__main__':
	main()
