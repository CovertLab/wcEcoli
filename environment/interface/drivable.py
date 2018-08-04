# This is an interface to enable any simulation to be driven by an external agent.
#
# (as python does not have real interfaces, it will be implemented by simply defining
# methods with the given signatures in the simulation class, so this class is not intended
# to be extended, but rather provide a definition of the methods to implement).
# 
# First, the simulation object, once created, will be initialized with an id and
# whatever state it requires (currently sim_data)
#
#     simulation.initialize(id, sim_data)
# 
# This `sim_data` object in our case is traditionally whatever is output by the fitter.
#
# Once the simulation has been initialized, it will receive a message containing
# its local environment and also the time step to be run until:
#
#     simulation.set_local_environment(molecule_ids, concentrations, run_for)
#
# Once this is invoked, the simulation will set its local environment and trigger
# the execution of the simulation until the given `run_for`. During this execution
# some means must be taken to record the deltas of the `molecule_ids` defined earlier.
#
# Once the simulation has executed for the given time period, the deltas will be read by calling `get_environment_change()`:
#
#     simulation.get_environment_change()
#
# This will return a dictionary of molecule names to integers representing the counts
# of how much each molecule changed with respect to its local environment.

from wholecell.sim.simulation import Simulation

import wholecell.states.environment


class Drivable(object):
	def initialize(self, id, sim_data):
		# state - dictionary of strings to arbitrary structures (whatever is required for initialization)
		# molecules_of_interest - array of strings (molecule names)
		# ----------------------------------------
		#   apply the initial state to the model and prepare it for execution
		#   (this is the state traditionally supplied by the output of the fitter),
		#   and set up the given `molecules_of_interest` to be tracked so that later the deltas
		#   for these molecules can be supplied.
		pass


	def set_local_environment(self, molecule_ids, concentrations, run_for):
		#
		# molecule_ids (str)
		# concentrations (float)
		# run_for (float) - time (sec) to run the WCM
		# ---------------------
		#   accept state of molecular concentrations and apply it to the local environment
		#   then, trigger execution for `run_for` steps (run simulation loop until that time step)

		Environment.setLocalEnvironment(molecule_ids, concentrations)
		Simulation.runIncremental(run_for)


	def get_environment_change(self):
		# return - dictionary of strings (molecule keys) to integers (molecule counts)
		# ----------------------------
		#   return a dictionary containing the deltas for each molecule previously passed
		#   to `initialize` as `molecules_of_interest`
		pass
