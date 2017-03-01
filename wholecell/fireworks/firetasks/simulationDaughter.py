import time

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.simulation import EcoliDaughterSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

@explicit_serialize
class SimulationDaughterTask(FireTaskBase):

	_fw_name = "SimulationDaughterTask"
	required_params = ["input_sim_data", "output_directory", "inherited_state_path"]
	optional_params = ["seed", "length_sec", "timestep_safety_frac", "timestep_max", "timestep_update_freq", "log_to_shell", "log_to_disk_every", "mass_distribution", "growth_rate_noise", "d_period_division"]

	def run_task(self, fw_spec):
		return