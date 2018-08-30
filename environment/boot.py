from __future__ import absolute_import, division, print_function

import os
import argparse

from agent.outer import Outer

from environment.batch_culture_nonspatial import EnvironmentBatchNonSpatial
from environment.two_dim_lattice import EnvironmentSpatialLattice

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli


class BootEnvironmentBatch(object):
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
				self.environment_dict[label].update({row["molecule id"]: row["concentration"].asNumber()}) # TODO (eran) pass units?

		# TODO (Eran) don't hardcode initial environment, get this from timeseries
		concentrations = self.environment_dict['minimal']

		self.environment = EnvironmentBatchNonSpatial(concentrations)
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
		self.outer = Outer(str(self.environment.agent_id), kafka_config, self.environment)


def main():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	parser = argparse.ArgumentParser(
		description='Boot the various agents for the environmental context simulation')

	parser.add_argument(
		'command',
		choices=['batch', 'lattice'],
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

	parser.add_argument(
		'--working-dir',
		default=os.getcwd(),
		help='the directory containing the project files'
	)

	args = parser.parse_args()
	kafka_config = {
		'host': args.kafka_host,
		'environment_control': args.environment_control,
		'simulation_receive': args.simulation_receive,
		'simulation_send': args.simulation_send,
		'subscribe_topics': []}

	if args.command == 'batch':
		outer = BootEnvironmentBatch(kafka_config)

	elif args.command == 'lattice':
		outer = BootEnvironmentSpatialLattice(kafka_config)


if __name__ == '__main__':
	main()
