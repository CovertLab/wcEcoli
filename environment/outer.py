import json
from kafka import KafkaProducer, KafkaConsumer

class Outer(object):
	def __init__(self, kafka, molecule_ids, run_for, concentrations, id='environment'):
		self.id = id
		self.kafka = kafka
		self.concentrations = concentrations
		self.molecule_ids = concentrations.keys()
		self.run_for = run_for
		self.time = 0

		self.producer = KafkaProducer(
			bootstrap_servers=self.kafka['host'],
			value_serializer=lambda v: json.dumps(v).encode('utf-8'))

		self.consumer = KafkaConsumer(
			bootstrap_servers=self.kafka['host'],
			value_deserializer=json.loads,
			auto_offset_reset='latest',
			group_id=self.id)

		self.consumer.subscribe(
			[self.kafka['simulation_send'],
			 self.kafka['environment_control']])

		self.simulations = {}

		for message in self.consumer:
			print(message)
			self.receive(message.value)

	def send(self, message):
		self.producer.send(
			self.kafka['simulation_receive'],
			message)

	def send_concentrations(self, concentrations, run_for):
		for id, simulation in self.simulations.iteritems():
			simulation['message_id'] += 1
			self.send({
				'id': id,
				'message_id': simulation['message_id'],
				'event': 'ENVIRONMENT_UPDATED',
				'molecule_ids': self.molecule_ids,
				'concentrations': concentrations,
				'run_for': run_for})

	def ready_to_advance(self, time):
		ready = True
		for id, simulation in self.simulations.iteritems():
			if simulation['message_id'] > simulation['last_message_id']:
				ready = False
				break

		return ready

	def receive(self, message):
		if message['event'] == 'SIMULATION_INITIALIZED':
			self.simulations[message['id']] = {
				'time': 0,
				'message_id': -1,
				'last_message_id': -1}

		if message['event'] == 'TRIGGER_EXECUTION':
			self.time += self.run_for
			self.send_concentrations(self.concentrations, self.run_for)

		if message['event'] == 'SIMULATION_ENVIRONMENT':
			self.simulations[message['id']]['changes'] = message['changes']
			self.simulations[message['id']]['time'] = message['time']
			self.simulations[message['id']]['last_message_id'] = message['message_id']

			if self.ready_to_advance(self.time):
				self.time += self.run_for
				self.send_concentrations(self.concentrations, self.run_for)

		if message['event'] == 'TRIGGER_SHUTDOWN':
			for id, simulation in self.simulations.iteritems():
				self.send({
					'id': id,
					'event': 'ENVIRONMENT_SHUTDOWN'})
