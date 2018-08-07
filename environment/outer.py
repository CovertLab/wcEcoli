import json

from environment.agent import Agent

class Outer(Agent):

	"""
	Outer: coordinate the communication between individual simulations and the larger
	environmental context.

	This class represents the larger environmental context for each of the individual simulation
	agents. The general flow is that each simulation will declare its existence to the 
	environment until the environment receives a signal from the control process to begin execution.
	Once this signal is received the environment will broadcast the local environmental
	concentrations to each simulation, at which point they will perform their local calculations
	and report back with their updated local environment. Once the environmental context receives
	a message from each simulation these changes are integrated, and the new updated environmental
	concentrations will be sent back. This loop will continue until the environment receives a 
	message to shutdown, when it will send a message to each simulation to shutdown and wait for
	acknowledgements, at which point it will shutdown itself.
	"""

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
		""" Send updated concentrations to each individual simulation agent. """

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
		"""
		Predicate to determine if the environment has heard back from all known simulations,
		in which case the environment can proceed to the next step.
		"""

		ready = True
		for id, simulation in self.simulations.iteritems():
			if simulation['message_id'] > simulation['last_message_id']:
				ready = False
				break

		return ready

	def receive(self, topic, message):
		"""
		Receive messages from associated simulation agents.

		The environment receives messages from both its associated simulation agents and also
		the control agent.

		Control messages:

		* TRIGGER_EXECUTION: Send messages to all associated simulations to begin execution.
		* SHUTDOWN_ENVIRONMENT: Send messages to simulations notifying them that the environment
		    is shutting, and wait for acknowledgement before exiting.

		Simulation messages:

		* SIMULATION_INITIALIZED: Received when a simulation agent is created. Upon the 
		    TRIGGER_EXECUTION message these are the simulations that will be communicated
		    with to perform the computation of the environmental changes.
		* SIMULATION_ENVIRONMENT: Received from each simulation when it has computed its 
		    environmental changes up to the specified `run_until`. The environment will wait
		    until it has heard from each simulation, integrate their changes and then calculate
		    the new local environment for each simulation and respond with an `ENVIRONMENT_UPDATED`
		    message.
		"""

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
