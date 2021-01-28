#! /usr/bin/env python

"""
Compare basal and delta probabilities of expression from two different sim_data
objects. Requires the path to
"""

import argparse
import os
import pickle

import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objs as go

from wholecell.utils.filepath import DEBUG_OUT_DIR


def parse_args():
	parser = argparse.ArgumentParser(
		description='Compare expression probabilities from two Parca outputs.')
	parser.add_argument('sim_dir', nargs=2,
		help='The two out/ sim-dirs to compare.')
	return parser.parse_args()

def load_sim_data(paths):
	sim_data = []
	for path in paths:
		with open(path, 'rb') as f:
			sim_data.append(pickle.load(f))
	return sim_data

def plot_comparison(sim_data1, sim_data2):
	t_reg1 = sim_data1.process.transcription_regulation
	t_reg2 = sim_data2.process.transcription_regulation

	# RNA and TF ids
	rnas1 = sim_data1.process.transcription.rna_data['id']
	rnas2 = sim_data2.process.transcription.rna_data['id']
	tfs1 = t_reg1.tf_ids
	tfs2 = t_reg2.tf_ids

	# Basal probabilities
	basal1 = {
		name: prob
		for name, prob in zip(rnas1, t_reg1.basal_prob)
		}
	basal2 = {
		name: prob
		for name, prob in zip(rnas2, t_reg2.basal_prob)
		}
	basal_ids = sorted(set(rnas1) | set(rnas2))
	x_basal = np.array([basal1.get(id_, 0) for id_ in basal_ids])
	y_basal = np.array([basal2.get(id_, 0) for id_ in basal_ids])
	min_basal = min(x_basal[x_basal > 0].min(), y_basal[y_basal > 0].min())
	x_basal_log = np.log10(x_basal + min_basal / 10)
	y_basal_log = np.log10(y_basal + min_basal / 10)

	# Delta probabilities
	delta1 = {
		(rnas1[i], tfs1[j]): v
		for i, j, v in zip(t_reg1.delta_prob['deltaI'], t_reg1.delta_prob['deltaJ'], t_reg1.delta_prob['deltaV'])
		}
	delta2 = {
		(rnas2[i], tfs2[j]): v
		for i, j, v in zip(t_reg2.delta_prob['deltaI'], t_reg2.delta_prob['deltaJ'], t_reg2.delta_prob['deltaV'])
		}
	delta_ids = sorted(delta1.keys() | delta2.keys())
	x_delta = np.array([delta1.get(id_, 0) for id_ in delta_ids])
	y_delta = np.array([delta2.get(id_, 0) for id_ in delta_ids])
	min_delta_pos = min(x_delta[x_delta > 0].min(), y_delta[y_delta > 0].min())
	x_delta_log_pos = np.log10(x_delta + min_delta_pos / 10)
	y_delta_log_pos = np.log10(y_delta + min_delta_pos / 10)
	min_delta_neg = -max(x_delta[x_delta < 0].max(), y_delta[y_delta < 0].max())
	x_delta_log_neg = np.log10(-x_delta + min_delta_neg / 10)
	y_delta_log_neg = np.log10(-y_delta + min_delta_neg / 10)

	# Plot probabilities
	fig = make_subplots(rows=3, cols=1)

	fig.add_traces([
		go.Scatter(x=x_basal_log, y=y_basal_log, mode='markers', text=basal_ids),
		go.Scatter(x=x_delta_log_pos, y=y_delta_log_pos, mode='markers', text=delta_ids),
		go.Scatter(x=x_delta_log_neg, y=y_delta_log_neg, mode='markers', text=delta_ids)],
		rows=[1, 2, 3], cols=[1, 1, 1],
		)

	os.makedirs(DEBUG_OUT_DIR, exist_ok=True)
	out = os.path.join(DEBUG_OUT_DIR, 'expression_probabilities.html')
	fig.write_html(out)

if __name__ == '__main__':
	args = parse_args()
	sim_data1, sim_data2 = load_sim_data(args.sim_dir)
	plot_comparison(sim_data1, sim_data2)
