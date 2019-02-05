from __future__ import division

import os
import time
import datetime
import numpy as np
from scipy.optimize import minimize

from helper_functions import get_v_translation, get_v_synthetase, get_saturation, get_saturation_CI, get_indices

# Describe paths
OPTIMIZE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(OPTIMIZE_DIR, "out")
PLOT_FILE = os.path.join(OPTIMIZE_DIR, "out", "optimization.pdf")

if not os.path.exists(OUT_DIR):
	os.mkdir(OUT_DIR)

# Constants
DEBUG = False

atp_id = "atp"
rib_id = "rib"

iteration = 0
n_condition = 2
n_aa = 2
f_charged_range = np.linspace(0.05, 0.95, 3)

def optimize_parameters(kb):
	def objective(x_log_space, constants, targets, v_weight, k_weight, parameter_ids):
		x = np.exp(x_log_space)

		# Parameters shared across all aa systems and conditions
		k_m_s_p = x[parameter_ids.index("k_m_s_p")]  # atp kM
		k_c_r = x[parameter_ids.index("k_c_r")]  # ribosome kcat

		# Constants shared across all aa systems and conditions
		fs = constants["f"]

		v_translation_errors = []
		v_steady_state_errors = []
		v_total_trna_errors = []

		# Compute condition dependent errors
		for j in range(n_condition):

			# Values shared across all aa systems
			c_r = constants[j][rib_id]
			c_p = constants[j][atp_id]

			# Compute translation flux
			c_t_all_i = x[get_indices(parameter_ids, "c_t", j = j)]
			c_ta_all_i = x[get_indices(parameter_ids, "c_ta", j = j)]
			k_m_r_t_all_i = x[get_indices(parameter_ids, "k_m_r_t")]
			k_m_r_ta_all_i = x[get_indices(parameter_ids, "k_m_r_ta")]

			v_rib = get_v_translation(k_c_r, c_r, c_t_all_i, c_ta_all_i, k_m_r_t_all_i, k_m_r_ta_all_i, fs, n_aa)

			# Compute error from translation flux target
			v_rib_error = get_flux_error(v_rib, targets[j][rib_id])
			v_translation_errors.append(v_rib_error)

			# Compute aa system specific errors
			for i in range(n_aa):

				# Amino acid system specific values
				k_c_s = x[get_indices(parameter_ids, "k_c_s", i = i)]
				if j == 0:
					c_s = constants[j]["synthetase"][i]
				else:
					c_s = x[get_indices(parameter_ids, "c_s", i = i, j = j)]
				c_t = x[get_indices(parameter_ids, "c_t", i = i, j = j)]
				k_m_s_t = x[get_indices(parameter_ids, "k_m_s_t", i = i)]
				c_a = constants[j]["aa"][i]
				k_m_s_a = x[get_indices(parameter_ids, "k_m_s_a", i = i)]
				c_ta = x[get_indices(parameter_ids, "c_ta", i = i, j = j)]
				f = fs[i]

				# Compute rate of aminoacylation
				v_ta = k_c_s * c_s * get_saturation(c_t, k_m_s_t) * get_saturation(c_a, k_m_s_a) * get_saturation(c_p, k_m_s_p)

				# Compute error from steady state
				v_steady_state_errors.append(get_flux_error(v_ta, v_rib * f))

				# Compute error from total trna
				trna_total_target = constants[j]["trna_total"][i]
				v_total_trna_errors.append(get_rel_error(c_t + c_ta, trna_total_target))

		# Compute v_error and k_error
		v_error = np.sum(v_translation_errors) + np.sum(v_steady_state_errors) + np.sum(v_total_trna_errors)
		k_error = 0

		total_error = (v_weight * v_error) + (k_weight * k_error)

		global iteration
		if DEBUG and iteration % 100 == 0:
			print("\ttotal error: {0:.2e}".format(total_error))
			print("\t\ttranslation  error: {0:.2e}".format(sum(v_translation_errors)))
			print("\t\tsteady state error: {0:.2e}".format(sum(v_steady_state_errors)))
			print("\t\ttrna total   error: {0:.2e}".format(sum(v_total_trna_errors)))
		iteration += 1

		return total_error

	def get_rel_error(x, x_target):
		return np.power(1. - (x / x_target), 2)

	def get_flux_error(v, v_target):
		return np.power(1. - (v / v_target), 2)

	def round_to_sig_dig(x):
		if x == 0:
			return 0
		else:
			return np.round(x, -int(np.floor(np.log10(x))))

	def get_initial_guess(constants, f_charged):
		# Generate initial guesses
		k_c_0 = 1e2 # 1/s (general enzyme kcat)

		parameter_ids = []
		initial_guess = []

		# Add kcat for ribosome
		parameter_ids.append("k_c_r")
		initial_guess.append(k_c_0)

		# Add km for ATP to tRNA synthetase
		parameter_ids.append("k_m_s_p")
		initial_guess.append(round_to_sig_dig(constants[0][atp_id]))

		for i in range(n_aa):
			## Add kinetic parameters (since they are condition independent)
			# kcat, synthetase
			parameter_ids.append("k_c_s__{0}".format(i))
			initial_guess.append(k_c_0)

			# kM, to synthetase, amino acid
			parameter_ids.append("k_m_s_a__{0}".format(i))
			initial_guess.append(constants[0]["aa"][i])

			## Add condition dependent parameters
			for j in range(n_condition):
				if j > 0:
					# Non-basal conditions optimize for synthetase
					parameter_ids.append("c_s__{0}__{1}".format(i, j))
					initial_guess.append(constants[j]["synthetase"][i])

				# tRNA
				trna_total = constants[j]["trna_total"][i]
				trna_free = round_to_sig_dig((1. - f_charged) * trna_total)
				trna_charged = round_to_sig_dig(f_charged * trna_total)

				parameter_ids.append("c_t__{0}__{1}".format(i, j))
				initial_guess.append(trna_free)

				parameter_ids.append("c_ta__{0}__{1}".format(i, j))
				initial_guess.append(trna_charged)

				if j == 0:
					# Use basal condition for kM initial guesses
					parameter_ids.append("k_m_s_t__{0}".format(i))
					initial_guess.append(trna_free)

					parameter_ids.append("k_m_r_t__{0}".format(i))
					initial_guess.append(trna_free)

					parameter_ids.append("k_m_r_ta__{0}".format(i))
					initial_guess.append(trna_charged)

		return parameter_ids, initial_guess

	def evaluate_optimization(optimized_parameters, parameter_ids):
		"""
		kcat: between 1e-2 and 1e6 (catalase is 1e5 1/s)
		kM and molecules: between 1e-9 and 1e-2 mol/L
		"""
		def evaluate(id_tag, lower_bound, upper_bound):
			val = optimized_parameters[get_indices(parameter_ids, id_tag)]
			in_range = (np.all(val > lower_bound) and np.all(val < upper_bound))
			order_of_magnitude_min = np.floor(np.log10(val.min()))
			order_of_magnitude_max = np.floor(np.log10(val.max()))
			return in_range, order_of_magnitude_min, order_of_magnitude_max

		k_c_in_range, k_c_min, k_c_max = evaluate("k_c", 1e-2, 1e6)
		k_m_in_range, k_m_min, k_m_max = evaluate("k_m", 1e-9, 1e2)
		c_t_in_range, c_t_min, c_t_max = evaluate("c_t", 1e-9, 1e2)
		c_ta_in_range, c_ta_min, c_ta_max = evaluate("c_ta", 1e-9, 1e2)
		c_s_in_range, c_s_min, c_s_max = evaluate("c_s", 1e-9, 1e2)

		evaluations_passed = [k_c_in_range, k_m_in_range, c_t_in_range, c_ta_in_range, c_s_in_range]

		for id_tag in ["k_c", "k_m", "c_t", "c_ta", "c_s"]:
			if id_tag == "k_c":
				lower_bound = 1e-2
				upper_bound = 1e6
			else:
				lower_bound = 1e-9
				upper_bound = 1e2
			in_range, order_of_magnitude_min, order_of_magnitude_max = evaluate(id_tag, lower_bound, upper_bound)
			print("\t{0}\t: {1} to {2}".format(id_tag, order_of_magnitude_min, order_of_magnitude_max))

			if not in_range:
				print("\t\t!OUTSIDE RANGE!")

		return np.all(evaluations_passed)

	# Build constants in the optimization problem
	constants = {}
	targets = {}
	for condition_index in range(n_condition):
		constants[condition_index] = {}
		constants[condition_index][rib_id] = kb.ribosome[condition_index]
		constants[condition_index][atp_id] = kb.atp[condition_index]
		constants[condition_index]["synthetase"] = kb.synthetase[condition_index]
		constants[condition_index]["trna_total"] = kb.trna_total[condition_index]
		constants[condition_index]["aa"] = kb.aa[condition_index]

		targets[condition_index] = {}
		if n_aa != 21:
			f_total = kb.f[condition_index][:n_aa].sum()
			targets[condition_index][rib_id] = f_total * kb.target_v_translation[condition_index]
		else:
			targets[condition_index][rib_id] = kb.target_v_translation[condition_index]

	if n_aa != 21:
		f_n_aa = kb.f[condition_index][:n_aa]
		f_scaled = f_n_aa / f_n_aa.sum()
		constants["f"] = f_scaled
	else:
		constants["f"] = kb.f[condition_index]

	# Optimize from several initial conditions
	opt_par_obj_values = []
	opt_par_evaluation = []
	opt_par_fs_charged = []
	opt_par_values = []
	for f_charged in f_charged_range:
		start_sec = time.clock()
		print("\t" + "-" * 76)
		print("\tfraction charged initial: {0}".format(f_charged))

		# Get initial guess
		parameter_ids, initial_guess = get_initial_guess(constants, f_charged)

		# Optimize
		result = minimize(
			objective,
			np.log(initial_guess),
			args = (constants, targets, 1, 0, parameter_ids), # v_weight = 1, k_weight = 0
		)

		# Evaluate objective function value
		opt_par_obj_values.append(result.fun)
		print("\n\tobjective: {0:.2e}".format(result.fun))

		# Evaluate optimized parameter values - order of magnitude reasonable
		opt_par = np.exp(result.x)
		opt_par_values.append(opt_par)
		opt_par_evaluation.append(evaluate_optimization(opt_par, parameter_ids))
		
		# Evaluate optimized parameter values - fraction charged
		c_t_indices = get_indices(parameter_ids, "c_t")
		opt_par_t = opt_par[c_t_indices]
		opt_par_ta = opt_par[get_indices(parameter_ids, "c_ta")]
		opt_par_f_charged = opt_par_ta / (opt_par_t + opt_par_ta)
		opt_par_fs_charged.append(opt_par_f_charged)

		print("\n\tfraction charged")
		for molecule, f_charged in zip(np.array(parameter_ids)[c_t_indices], opt_par_f_charged):
			print("\t\t{0}: {1:.2f}".format(molecule, f_charged))

		end_sec = time.clock()
		elapsed = datetime.timedelta(seconds = (end_sec - start_sec))

		print("\trun in {}h {}m {}s".format(*str(elapsed).split(':')))

	return (parameter_ids, opt_par_values, opt_par_obj_values, opt_par_evaluation, opt_par_fs_charged)

