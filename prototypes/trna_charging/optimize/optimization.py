from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import minimize

from wholecell.utils import units
from knowledge_base import kb

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

def get_v_rib(k_r, r, all_t, all_ta, all_k_r_t, all_k_r_ta, fs):
	data_by_aa_system = [(all_ta[i], all_t[i], all_k_r_ta[i], all_k_r_t[i]) for i in range(n_aa)]
	v_rib = k_r * r * np.prod(np.power([get_saturation_CI(x) for x in data_by_aa_system], fs), axis = 0)
	return v_rib

def get_saturation(substrate, k):
	return substrate / (substrate + k)

def get_saturation_CI(x):
	"""
    Model: competitive inhibition
    Note: "substrate1" is the substrate on the numerator
    """
	substrate1, substrate2, k1, k2 = x
	return (substrate1 / k1) / (1 + (substrate1 / k1) + (substrate2 / k2))


def optimize_parameters(kb):
	def objective(x_log_space, constants, targets, v_weight, k_weight):
		x = np.exp(x_log_space)

		# Parameters shared across all aa systems
		k_p = x[id_to_indices[atp_id]] # atp kM
		k_r = x[id_to_indices[rib_id]] # ribosome kcat

		# Constants shared across all aa systems
		r = constants[rib_id]
		p = constants[atp_id]
		fs = constants["f"]

		# Compute ribosome flux
		all_t = x[id_to_indices["t"]]
		all_ta = x[id_to_indices["ta"]]
		all_k_r_t = x[id_to_indices["k_r_t"]]
		all_k_r_ta = x[id_to_indices["k_r_ta"]]
		v_rib = get_v_rib(k_r, r, all_t, all_ta, all_k_r_t, all_k_r_ta, fs)

		# Compute error from ribosome target
		v_rib_error = get_flux_error(v_rib, targets[rib_id])

		# Compute aa system specific errors
		v_ta_all = []
		v_ta_error = []
		trna_total_error = []
		for i in range(n_aa):
			# Amino acid system specific parameters
			k_s, t, ta, k_t, k_a, k_r_t, k_r_ta = x[id_to_indices[i]]

			# Amino acid system specific constants
			s = constants["synthetase"][i]
			a = constants["aa"][i]
			f = constants["f"][i]

			# Compute rate of aminoacylation
			v_ta = k_s * s * get_saturation(t, k_t) * get_saturation(a, k_a) * get_saturation(p, k_p)
			v_ta_all.append(v_ta)

			# Compute error from steady state
			v_ta_error.append(get_flux_error(v_ta, v_rib * f))

			# Compute error from total trna
			trna_total = constants["trna_total"][i]
			trna_total_error.append(get_rel_error(t + ta, trna_total))
		v_ta_error = np.array(v_ta_error)
		trna_total_error = np.array(trna_total_error)

		# Compute total error
		v_error = v_rib_error + v_ta_error.sum() + trna_total_error.sum()
		k_error = 0
		total_error = (v_weight * v_error) + (k_weight * k_error)

		if DEBUG:
			print "total error: {0:.2e}".format(total_error)
			print "\tv_rib_error: {0:.2e}".format(v_rib_error)
			print "\tv_ta_error sum: {0:.2e}".format(v_ta_error.sum())
			print "\ttrna_total_error sum: {0:.2e}".format(trna_total_error.sum())

		return total_error

	def get_rel_error(x, x_target):
		return np.power(1. - (x / x_target), 2)

	def get_flux_error(v, v_target):
		return np.power(1. - (v / v_target), 2)

	def as_moles_per_liter(data_array):
		return np.array([x.asNumber(units.mol / units.L) for x in data_array])

	def round_to_sig_dig(x):
		if x == 0:
			return 0
		else:
			return np.round(x, -int(np.floor(np.log10(x))))

	def get_initial_guess(aa, trna_total, atp, f_charged):
		aa = as_moles_per_liter(aa)
		trna_total = as_moles_per_liter(trna_total)

		# Generate initial guesses
		k_s_0 = 1e2 # 1/s (tRNA synthetase kcat)
		k_p_0 = round_to_sig_dig(atp.asNumber(units.mol / units.L)) # ATP kM to tRNA synthetase

		x_lin_space = []
		id_to_indices = {}
		for molecule_id in ["t", "ta", "k_r_t", "k_r_ta", "kcat", "km", "km_substrate"]:
			id_to_indices[molecule_id] = []

		x_lin_space.append(k_p_0)
		id_to_indices[atp_id] = 0

		x_lin_space.append(k_s_0) # 1/s (ribosome kcat)
		id_to_indices[rib_id] = 1

		for i in range(n_aa):
			t_0 = round_to_sig_dig((1. - f_charged) * trna_total[i]) # free form
			ta_0 = round_to_sig_dig(f_charged * trna_total[i]) # charged form
			k_t_0 = t_0 #* 1e0
			k_a_0 = round_to_sig_dig(aa[i]) #* 1e1
			k_r_t_0 = t_0
			k_r_ta_0 = ta_0
			x_lin_space += [k_s_0, t_0, ta_0, k_t_0, k_a_0, k_r_t_0, k_r_ta_0]

			# Create mappings for use in objective function
			id_to_indices[i] = range(2 + 7*i, 2 + 7*i + 7)
			id_to_indices["t"].append(id_to_indices[i][1])
			id_to_indices["ta"].append(id_to_indices[i][2])
			id_to_indices["k_r_t"].append(id_to_indices[i][5])
			id_to_indices["k_r_ta"].append(id_to_indices[i][6])

			# Create mappings for use in post-optimization analysis
			id_to_indices["kcat"].append(id_to_indices[i][0])
			id_to_indices["km"] += id_to_indices[i][3:]
			id_to_indices["km_substrate"] += np.array(id_to_indices[i])[[1, 2, 1, 2]].tolist()

		return x_lin_space, id_to_indices

	def evaluate_optimization(result, id_to_indices):
		"""
		expected ranges
			kcat: between 1e-2 and 1e6 (catalase is 1e5 1/s)
			kM and substrate: between 1e-9 and 1e-1 mol/L (glutamate, 96 mM)
		"""
		optimized_parameters = np.exp(result.x)

		# Evaluate kcats
		opt_kcats = optimized_parameters[id_to_indices["kcat"]]
		opt_kcats = np.append(opt_kcats, optimized_parameters[id_to_indices[rib_id]])
		kcats_in_range = (np.all(opt_kcats > 1e-2) and np.all(opt_kcats < 1e6))

		# Evaluate kMs
		opt_kms = optimized_parameters[id_to_indices["km"]]
		opt_kms_order = np.floor(np.log10(opt_kms))
		kms_in_range = np.all(np.all(opt_kms_order >= -9) and np.all(opt_kms_order <= -1))

		# Evaluate substrates
		opt_sub = optimized_parameters[id_to_indices["km_substrate"]]
		opt_sub_order = np.floor(np.log10(opt_sub))
		sub_in_range = np.all(np.all(opt_sub_order >= -9) and np.all(opt_sub_order <= -1))

		if not (kcats_in_range and kms_in_range and sub_in_range):
			for name, res, opt_order in zip(
					["kcats", "kms", "sub"],
					[kcats_in_range, kms_in_range, sub_in_range],
					[opt_kcats, opt_kms_order, opt_sub_order]
			):
				print "\t{0}:\t{1}".format(name, res)
				if not res:
					print "\t\t{0} range: {1} to {2}".format(name, opt_order.min(), opt_order.max())

		return kcats_in_range and kms_in_range and sub_in_range

	# Get constants
	constants = {}
	constants[rib_id] = kb.ribosome.asNumber(units.mol / units.L)
	constants[atp_id] = kb.atp.asNumber(units.mol / units.L)
	constants["synthetase"] = as_moles_per_liter(kb.synthetase)
	constants["trna_total"] = as_moles_per_liter(kb.trna_total)
	constants["aa"] = as_moles_per_liter(kb.aa)
	constants["f"] = kb.f
	if n_aa != kb.aa.shape[0]:
		# Optimizing sub-set of full aa system - rescale f to sum to 1
		f_subset = constants["f"][:n_aa]
		f_scaled = f_subset / f_subset.sum()
		constants["f"] = f_scaled

	# Get targets
	targets = {}
	targets[rib_id] = kb.target_v_ribosome.asNumber(units.mol / units.L / units.s)
	if n_aa != kb.aa.shape[0]:
		# Optimizing sub-set of full aa system - rescale translation rate target
		target_rib_original = targets[rib_id]
		target_rib_adjusted = (constants["f"][:n_aa] * target_rib_original).sum()
		targets[rib_id] = target_rib_adjusted

	# Optimizing from several initial conditions
	obj_fun_values = []
	passed_eval = []
	opt_f_charged = []
	for f_charged in f_charged_range:
		print "\t" + "-"*75
		print "\tinitialized with fraction charged: {0}".format(f_charged)

		# Get initial guess
		initial_guess, id_to_indices = get_initial_guess(kb.aa, kb.trna_total, kb.atp, f_charged)

		# Optimize
		result = minimize(
			objective,
			np.log(initial_guess),
			args = (constants, targets, 1, 0), # v_weight = 1, k_weight = 0
		)

		obj_fun_values.append(result.fun)
		print "\tobjective value: {0:.2e}".format(result.fun)

		# Evaluate result
		passed_eval.append(evaluate_optimization(result, id_to_indices))
		opt_par = np.exp(result.x)
		opt_par_t = opt_par[id_to_indices["t"]]
		opt_par_ta = opt_par[id_to_indices["ta"]]
		opt_par_f_charged = opt_par_ta / (opt_par_t + opt_par_ta)
		opt_f_charged.append(opt_par_f_charged)

	return (obj_fun_values, passed_eval, opt_f_charged)


# Plot optimization solutions
fig, axesList = plt.subplots(2, 1, figsize = (8.5, 11))
ax1, ax2 = axesList

n_aa = 21
print "n_aa: {0}".format(n_aa)

f_charged_range = np.linspace(0.05, 0.95, 10)
obj_fun_values, passed_eval, opt_f_charged = optimize_parameters(kb)
ax1.scatter(f_charged_range, np.log10(obj_fun_values), marker = "o")

colors = [plt.cm.Spectral(i) for i in np.linspace(0, 1, n_aa)]

for i, trial in enumerate(opt_f_charged):
	ax2.scatter(np.array([f_charged_range[i]] * n_aa) + np.linspace(0, 0.01, n_aa), trial, c = colors, s = 8)

patches = [mpatches.Patch(color = colors[i], label = kb.aa_ids[i]) for i in range(len(colors))]
plt.legend(handles = patches, loc = (1.1, 0))
plt.subplots_adjust(right = 0.65)

ax1.set_xticks(f_charged_range)
ax1.set_xlabel("fraction charged initial condition")
ax1.set_ylabel("log10(value of objective)")

ax2.set_facecolor((0.2, 0.2, 0.2))
ax2.set_xticks(f_charged_range)
ax2.set_xlabel("fraction charged initial condition")
ax2.set_ylabel("fraction charged in optimized result")
plt.savefig(PLOT_FILE)
plt.close("all")
