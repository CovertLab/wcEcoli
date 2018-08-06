import json
from kafka import KafkaProducer, KafkaConsumer

class Inner(object):
	def __init__(self, id, simulation, kafka):
		self.id = id
		self.simulation = simulation
		self.kafka = kafka
		self.simulation.track_environment_change()

		self.producer = KafkaProducer(
			bootstrap_servers=self.kafka['host'],
			value_serializer=lambda v: json.dumps(v).encode('utf-8'))

		self.consumer = KafkaConsumer(
			bootstrap_servers=self.kafka['host'],
			value_deserializer=json.loads,
			auto_offset_reset='latest',
			group_id=('simulation-' + str(id)))

		self.consumer.subscribe(
			[self.kafka['simulation_receive']])

		self.send({
			'event': 'SIMULATION_INITIALIZED',
			'id': id})

		for message in self.consumer:
			if message.value['id'] == self.id:
				print(message)
				self.receive(message.value)
			self.consumer.commit()

	def send(self, message):
		self.producer.send(
			self.kafka['simulation_send'],
			message)

	def receive(self, message):
		if message['event'] == 'ENVIRONMENT_UPDATED':
			self.simulation.set_local_environment(
				message['molecule_ids'],
				message['concentrations'])

			self.simulation.run_incremental(message['run_for'] + self.simulation.time())

			stop = self.simulation.time()
			changes = self.simulation.get_environment_change()

			self.send({
				'event': 'SIMULATION_ENVIRONMENT',
				'id': self.id,
				'message_id': message['message_id'],
				'time': stop,
				'changes': changes})

		if message['event'] == 'ENVIRONMENT_SHUTDOWN':
			self.simulation.finalize()
