import pickle
import time

from fireworks import FiretaskBase, explicit_serialize
from models.ecoli.sim.simulation import EcoliSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS


@explicit_serialize
class SimulationTask(FiretaskBase):

	_fw_name = "SimulationTask"
	required_params = [
		"input_sim_data",
		"output_directory"]
	optional_params = [
		"seed",
		"timeline",
		"length_sec",
		"timestep_safety_frac",
		"timestep_max",
		"timestep_update_freq",
		"adjust_timestep_for_charging",
		"log_to_shell",
		"log_to_disk_every",
		"jit",
		"mass_distribution",
		"d_period_division",
		'variable_elongation_transcription',
		'variable_elongation_translation',
		"translation_supply",
		"trna_charging",
		"aa_supply_in_charging",
		"ppgpp_regulation",
		"disable_ppgpp_elongation_inhibition",
		"superhelical_density",
		"recycle_stalled_elongation",
		"mechanistic_replisome",
		"mechanistic_translation_supply",
		"mechanistic_aa_transport",
		"trna_attenuation",
		"raise_on_time_limit"]

	def _get_default(self, key, default_key=''):
		return self.get(key, DEFAULT_SIMULATION_KWARGS[default_key or key])

	def run_task(self, fw_spec):
		print("%s: Running simulation" % time.ctime())

		# load the sim_data from the output of the parameter calculator (parca)
		# TODO(spanglry): make the parca output JSON and this load from JSON instead
		with open(self["input_sim_data"], "rb") as input_sim_data:
			sim_data = pickle.load(input_sim_data)

		options = {}

		options["simData"] = sim_data
		options["outputDir"] = self["output_directory"]
		options["logToDisk"] = True

		options["seed"] = int(self._get_default("seed"))
		options["timeline"] = self._get_default("timeline")
		options["lengthSec"] = self._get_default("length_sec", "lengthSec")
		options["timeStepSafetyFraction"] = self._get_default("timestep_safety_frac", "timeStepSafetyFraction")
		options["maxTimeStep"] = self._get_default("timestep_max", "maxTimeStep")
		options["updateTimeStepFreq"] = self._get_default("timestep_update_freq", "updateTimeStepFreq")
		options["adjust_timestep_for_charging"] = self._get_default("adjust_timestep_for_charging")
		options["logToShell"] = self._get_default("log_to_shell", "logToShell")
		options["logToDiskEvery"] = self._get_default("log_to_disk_every", "logToDiskEvery")
		options["jit"] = self._get_default("jit")
		options["massDistribution"] = self._get_default("mass_distribution", "massDistribution")
		options["dPeriodDivision"] = self._get_default("d_period_division", "dPeriodDivision")
		options["translationSupply"] = self._get_default("translation_supply", "translationSupply")
		options["variable_elongation_transcription"] = self._get_default("variable_elongation_transcription")
		options["variable_elongation_translation"] = self._get_default("variable_elongation_translation")
		options["trna_charging"] = self._get_default("trna_charging")
		options["aa_supply_in_charging"] = self._get_default("aa_supply_in_charging")
		options["ppgpp_regulation"] = self._get_default("ppgpp_regulation")
		options["disable_ppgpp_elongation_inhibition"] = self._get_default("disable_ppgpp_elongation_inhibition")
		options["superhelical_density"] = self._get_default("superhelical_density")
		options["recycle_stalled_elongation"] = self._get_default("recycle_stalled_elongation")
		options["mechanistic_replisome"] = self._get_default("mechanistic_replisome")
		options["mechanistic_translation_supply"] = self._get_default("mechanistic_translation_supply")
		options["mechanistic_aa_transport"] = self._get_default("mechanistic_aa_transport")
		options["trna_attenuation"] = self._get_default("trna_attenuation")
		options["raise_on_time_limit"] = self._get_default("raise_on_time_limit")

		sim = EcoliSimulation(**options)

		sim.run()
