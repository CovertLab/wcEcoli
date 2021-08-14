"""
Set the concentration of ppGpp to be constant during simulations at different
concentrations.

Modifies:
	sim_data.growth_rate_parameters.get_ppGpp_conc
	sim_data.constants.k_RelA_ppGpp_synthesis
	sim_data.constants.k_SpoT_ppGpp_degradation
	sim_data.constants.k_SpoT_ppGpp_synthesis

Expected variant indices (dependent on FACTORS):
	0-3: lower concentrations of ppGpp
	4: control
	5-8: higher concentrations of ppGpp
"""

from wholecell.utils import units


FACTORS = [0.1, 0.2, 0.5, 0.8, 1, 1.2, 1.5, 2, 5]


class ppGpp():
	def __init__(self, conc):
		self.conc = conc

	def get_ppGpp_conc(self, doubling_time):
		"""
		Function to replace get_ppGpp_conc in growth_rate_parameters.  Needs to
		have the same function signature but we want to fix the output
		conentration to a new value.
		"""
		return self.conc


def ppgpp_conc(sim_data, index):
	control_conc = sim_data.growth_rate_parameters.get_ppGpp_conc(sim_data.doubling_time)
	factor = FACTORS[index]
	new_conc = factor * control_conc

	# Replace function to only return the desired concentration
	sim_data.growth_rate_parameters.get_ppGpp_conc = ppGpp(new_conc).get_ppGpp_conc

	# Set rate constants to 0 to prevent ppGpp concentration changes during sim
	sim_data.constants.k_RelA_ppGpp_synthesis *= 0
	sim_data.constants.k_SpoT_ppGpp_degradation *= 0
	sim_data.constants.k_SpoT_ppGpp_synthesis *= 0

	return dict(
		shortName=f'ppGpp:{factor}x',
		desc=f'ppGpp conc adjusted by {factor}x to {new_conc.asNumber(units.umol / units.L)} uM.'
		), sim_data
