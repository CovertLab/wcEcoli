import json

from environment.outer import Outer
from environment.inner import Inner
from environment.stub import SimulationStub
from kafka import KafkaProducer

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

class EnvironmentControl(object):
	def __init__(self):
		self.kafka = default_kafka_config
		self.producer = KafkaProducer(
			bootstrap_servers=self.kafka['host'],
			value_serializer=lambda v: json.dumps(v).encode('utf-8'))

	def send(self, message):
		self.producer.send(
			self.kafka['environment_control'],
			message)
