import json

from environment.agent import Agent
from environment.outer import Outer
from environment.inner import Inner
from environment.stub import SimulationStub
from confluent_kafka import Producer

default_kafka_config = {
	'host': 'localhost:9092',
	'simulation_send': 'environment_listen',
	'simulation_receive': 'environment_broadcast',
	'environment_control': 'environment_control'}

class BootOuter(object):
	def __init__(self):
		self.outer = Outer(
			default_kafka_config,
			['yellow', 'green', 'red', 'blue'],
			1,
			{'yellow': 5,
			 'green': 11,
			 'red': 44,
			 'blue': 12})

class BootInner(object):
	def __init__(self, id):
		self.id = id
		self.simulation = SimulationStub()
		self.inner = Inner(
			self.id,
			self.simulation,
			default_kafka_config)

class EnvironmentControl(Agent):
	def __init__(self):
		id = 'environment_control'
		kafka = default_kafka_config.copy()
		kafka['subscribe_topics'] = []

		super(EnvironmentControl, self).__init__(id, kafka)

	def trigger_execution(self):
		self.send(self.kafka['environment_control'], {
			'event': 'TRIGGER_EXECUTION'})

	def shutdown_environment(self):
		self.send(self.kafka['environment_control'], {
			'event': 'SHUTDOWN_ENVIRONMENT'})
