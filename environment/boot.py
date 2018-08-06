import json

from environment.outer import Outer
from environment.inner import Inner
from environment.stub import SimulationStub
from confluent_kafka import Producer

default_kafka_config = {
	'host': 'localhost:9092',
	'simulation_send': 'environment_listen',
	'simulation_receive': 'environment_broadcast',
	'environment_control': 'environment_control'}

def delivery_report(err, msg):
	if err is not None:
		print('message delivery failed: {}'.format(msg))
		print('failed message: {}'.format(err))

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
		self.producer = Producer({
			'bootstrap.servers': self.kafka['host']})

	def send(self, message):
		self.producer.flush()
		self.producer.poll(0)
		self.producer.produce(
			self.kafka['environment_control'],
			json.dumps(message).encode('utf-8'),
			callback=delivery_report)
