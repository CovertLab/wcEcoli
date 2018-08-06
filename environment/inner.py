import json
from confluent_kafka import Producer, Consumer, KafkaError

def delivery_report(err, msg):
	if err is not None:
		print('message delivery failed: {}'.format(msg))
		print('failed message: {}'.format(err))

class Inner(object):
	def __init__(self, id, simulation, kafka):
		self.id = id
		self.simulation = simulation
		self.kafka = kafka
		self.simulation.track_environment_change()

		self.producer = Producer({
			'bootstrap.servers': self.kafka['host']})

		self.consumer = Consumer({
			'bootstrap.servers': self.kafka['host'],
			'enable.auto.commit': True,
			'group.id': 'simulation-' + str(id),
			'default.topic.config': {
				'auto.offset.reset': 'smallest'}})

		self.consumer.subscribe(
			[self.kafka['simulation_receive']])

		self.send({
			'event': 'SIMULATION_INITIALIZED',
			'id': id})

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
			if message['id'] == self.id:
				print(message)
				if message['event'] == 'ENVIRONMENT_SHUTDOWN':
					self.finalize()
					running = False
				else:
					self.receive(message)

	def finalize(self):
		self.simulation.finalize()
		self.consumer.close()

	def send(self, message):
		self.producer.flush()
		self.producer.poll(0)
		self.producer.produce(
			self.kafka['simulation_send'],
			json.dumps(message).encode('utf-8'),
			callback=delivery_report)

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
