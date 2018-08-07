import json

from environment.agent import Agent

class Inner(Agent):
	def __init__(self, id, simulation, kafka):
		self.simulation = simulation
		self.simulation.initialize_local_environment()
		kafka['subscribe_topics'] = [kafka['simulation_receive']]

		super(Inner, self).__init__(id, kafka)

	def initialize(self):
		self.send(self.kafka['simulation_send'], {
			'event': 'SIMULATION_INITIALIZED',
			'id': self.id})

	def finalize(self):
		self.simulation.finalize()

	def receive(self, topic, message):
		if message['id'] == self.id:
			print(topic + ': ' + str(message))

			if message['event'] == 'ENVIRONMENT_UPDATED':
				self.simulation.set_local_environment(
					message['molecule_ids'],
					message['concentrations'])

				self.simulation.run_incremental(message['run_for'] + self.simulation.time())

				stop = self.simulation.time()
				changes = self.simulation.get_environment_change()

				self.send(self.kafka['simulation_send'], {
					'event': 'SIMULATION_ENVIRONMENT',
					'id': self.id,
					'message_id': message['message_id'],
					'time': stop,
					'changes': changes})

			if message['event'] == 'SHUTDOWN_SIMULATION':
				self.send(self.kafka['simulation_send'], {
					'event': 'SIMULATION_SHUTDOWN',
					'id': self.id})

				self.shutdown()
