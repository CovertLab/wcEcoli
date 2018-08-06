import json

from environment.agent import Agent

class Outer(Agent):
	def __init__(self, kafka, molecule_ids, run_for, concentrations, id='environment'):
		self.concentrations = concentrations
		self.molecule_ids = concentrations.keys()
		self.run_for = run_for
		self.time = 0
		self.simulations = {}

		kafka['subscribe_topics'] = [
			kafka['simulation_send'],
			kafka['environment_control']]

		super(Outer, self).__init__(id, kafka)

	def finalize(self):
		print('environment shutting down')

	def send_concentrations(self, concentrations, run_for):
		for id, simulation in self.simulations.iteritems():
			simulation['message_id'] += 1
			self.send(self.kafka['simulation_receive'], {
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

	def receive(self, topic, message):
		print(topic + ': ' + str(message))

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

		if message['event'] == 'SHUTDOWN_ENVIRONMENT':
			for id, simulation in self.simulations.iteritems():
				self.send(self.kafka['simulation_receive'], {
					'id': id,
					'event': 'SHUTDOWN_SIMULATION'})

		if message['event'] == 'SIMULATION_SHUTDOWN':
			gone = self.simulations.pop(message['id'], {})
			print('simulation shutdown: ' + str(gone))

			if not any(self.simulations):
				self.shutdown()
