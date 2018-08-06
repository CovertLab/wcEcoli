import json
from confluent_kafka import Producer, Consumer, KafkaError

def delivery_report(err, msg):
	if err is not None:
		print('message delivery failed: {}'.format(msg))
		print('failed message: {}'.format(err))

class Outer(object):
	def __init__(self, kafka, molecule_ids, run_for, concentrations, id='environment'):
		self.id = id
		self.kafka = kafka
		self.concentrations = concentrations
		self.molecule_ids = concentrations.keys()
		self.run_for = run_for
		self.time = 0
		self.simulations = {}

		self.producer = Producer({
			'bootstrap.servers': self.kafka['host']})

		self.consumer = Consumer({
			'bootstrap.servers': self.kafka['host'],
			'enable.auto.commit': True,
			'group.id': str(id),
			'default.topic.config': {
				'auto.offset.reset': 'smallest'}})

		self.consumer.subscribe(
			[self.kafka['simulation_send'],
			 self.kafka['environment_control']])

		self.poll()

	def poll(self):
		running = True
		while running:
			raw = self.consumer.poll()
			if raw is None:
				continue
			if raw.error():
				if raw.error().code() == KafkaError._PARTITION_EOF:
					continue
				else:
					print(raw.error())
					running = False

			message = json.loads(raw.value().decode('utf-8'))
			print(message)
			if message['event'] == 'TRIGGER_SHUTDOWN':
				self.finalize()
				running = False
			else:
				self.receive(message)

	def finalize(self):
		for id, simulation in self.simulations.iteritems():
			self.send({
				'id': id,
				'event': 'ENVIRONMENT_SHUTDOWN'})

	def send(self, message):
		self.producer.flush()
		self.producer.poll(0)
		self.producer.produce(
			self.kafka['simulation_receive'],
			json.dumps(message).encode('utf-8'),
			callback=delivery_report)

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
