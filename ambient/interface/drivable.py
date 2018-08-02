# This is an interface to enable any simulation to be driven by an external agent.
#
# (as python does not have real interfaces, it will be implemented by simply defining
# methods with the given signatures in the simulation class, so this class is not intended
# to be extended, but rather provide a definition of the methods to implement).
# 
# First, the simulation object, once created, will be initialized with whatever state it requires,
# and also an array detailing which molecules the larger simulation will be concerned with:
#
#     simulation.initialize(state, molecules_of_interest)
# 
# This `state` object in our case is traditionally whatever is output by the fitter.
# To keep the interface clean let's have this be the instantiated object and not the pickled file
# (so, we do the work of reading an unpickling the file before it is passed to `initialize`
#
# Once the simulation has been initialized, it will receive a message containing its local environment
# and also the time step to be run until:
#
#     simulation.receive_environment(environment, time_to_run)
#
# Once this is invoked, the simulation will set its local environment and trigger the execution
# of the simulation until the given `time_to_run`. During this execution some means must be taken to record
# the deltas of the `molecules_of_interest` defined earlier.
#
# Once the simulation has executed for the given time period, the deltas will be read by calling `read_deltas()`:
#
#     simulation.read_deltas()
#
# This will return a dictionary of molecule names to integers representing the counts of how much each molecule
# changed with respect to its local environment.

class Drivable(object):
	def initialize(self, state, molecules_of_interest):
		# state - dictionary of strings to arbitrary structures (whatever is required for initialization)
		# molecules_of_interest - array of strings (molecule names)
		# ----------------------------------------
		#   apply the initial state to the model and prepare it for execution
		#   (this is the state traditionally supplied by the output of the fitter),
		#   and set up the given `molecules_of_interest` to be tracked so that later the deltas
		#   for these molecules can be supplied.
		pass

	def receive_environment(self, environment, time_to_run):
		# environment - dictionary of molecules to counts
		# time_to_run - integer
		# ---------------------
		#   accept environmental information and apply it to the local environment
		#   then, trigger execution for `time_to_run` steps (run simulation loop until that time step),
		#   tracking the deltas for the `molecules_of_interest`
		pass

	def read_deltas(self):
		# return - dictionary of strings (molecule keys) to integers (molecule counts)
		# ----------------------------
		#   return a dictionary containing the deltas for each molecule previously passed
		#   to `initialize` as `molecules_of_interest`
		pass
