import json
from confluent_kafka import Producer, Consumer, KafkaError

def delivery_report(err, msg):
	if err is not None:
		print('message delivery failed: {}'.format(msg))
		print('failed message: {}'.format(err))

class Agent(object):
	def __init__(self, id, kafka):
		self.id = id
		self.kafka = kafka

		self.producer = Producer({
			'bootstrap.servers': self.kafka['host']})

		self.has_consumer = len(self.kafka['subscribe_topics']) > 0
		if self.has_consumer:
			self.consumer = Consumer({
				'bootstrap.servers': self.kafka['host'],
				'enable.auto.commit': True,
				'group.id': 'simulation-' + str(id),
				'default.topic.config': {
					'auto.offset.reset': 'latest'}})

		self.initialize()

		if self.has_consumer:
			self.consumer.subscribe(
				self.kafka['subscribe_topics'])

			self.poll()

	def initialize(self):
		pass

	def poll(self):
		self.running = True
		while self.running:
			raw = self.consumer.poll()
			if raw is None:
				continue
			if raw.error():
				if raw.error().code() == KafkaError._PARTITION_EOF:
					continue
				else:
					print(raw.error())
					self.running = False

			message = json.loads(raw.value().decode('utf-8'))

			if message['event'] == 'GLOBAL_SHUTDOWN':
				self.shutdown()
			else:
				self.receive(raw.topic(), message)

	def send(self, topic, message):
		self.producer.flush()
		self.producer.poll(0)
		self.producer.produce(
			topic,
			json.dumps(message).encode('utf-8'),
			callback=delivery_report)

	def receive(self, topic, message):
		pass

	def shutdown(self):
		self.running = False
		self.finalize()

		self.producer.flush()
		self.producer.poll(0)

		if self.has_consumer:
			self.consumer.commit()
			self.consumer.close()

	def finalize(self):
		pass
