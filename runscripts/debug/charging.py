#! /usr/bin/env python
"""
Tools and analysis to debug charging problems.
"""

from __future__ import annotations

import argparse
import os
import pickle
from typing import Dict, Optional, Tuple
import webbrowser

import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import plotly.graph_objs as go
import plotly.subplots

from models.ecoli.processes.polypeptide_elongation import (calculate_trna_charging,
	CONC_UNITS, get_charging_params)
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, scriptBase


# Set of charging parameters that should not be modified by a slider
CONSTANT_PARAMS = {'charging_mask'}

PORT = 8050

# Element IDs
GRAPH_ID = 'graph-id'
LOW_TIMESTEP = 'low-timestep'
HIGH_TIMESTEP = 'high-timestep'


class ChargingDebug(scriptBase.ScriptBase):
	def define_parameters(self, parser):
		super().define_parameters(parser)

		# Path args
		self.define_parameter_sim_dir(parser)
		self.define_parameter_variant_index(parser)
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='The initial simulation number (int). The value will get'
				 ' formatted as a subdirectory name like "000000". Default = 0.')
		parser.add_argument('-g', '--generation', type=int, default=0,
			help='The generation number (int). The value will get formatted'
				 ' as a subdirectory name like "generation_000000". Default = 0.')
		parser.add_argument('-d', '--daughter', type=int, default=0,
			help='The daughter number (int). The value will get formatted as'
				 ' a subdirectory name like "000000". Default = 0.')

		# Sim options
		self.define_parameter_bool(parser, 'variable_elongation_translation',
			default_key='variable_elongation_translation',
			help='set if sims were run with variable_elongation_translation')

		# Debug options
		parser.add_argument('--validation', type=int, default=1,
			help='Number of time steps to run for validation. If < 0, will run all.')
		parser.add_argument('--interactive', action='store_true',
			help='If set, runs interactive analysis plots for debugging.')

	def update_args(self, args):
		super().update_args(args)

		# Extract data from args
		variant_dir_name, _, _ = args.variant_dir
		seed_str = '%06d' % (args.seed,)
		gen_str = 'generation_%06d' % (args.generation,)
		daughter_str = '%06d' % (args.daughter,)
		dirs = os.path.join(seed_str, gen_str, daughter_str)
		input_variant_directory = os.path.join(args.sim_path, variant_dir_name)

		# Set paths from args
		args.sim_data_file = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		args.sim_out_dir = os.path.join(input_variant_directory, dirs, 'simOut')

	def load_data(self, sim_data_file: str, sim_out_dir: str) -> None:
		"""
		Loads sim_data and simulation output data from files and saves it as
		instance variables.

		Args:
			sim_data_file: path to the sim_data file for the simulation
			sim_out_dir: path to the simOut dir for the simulation
		"""

		with open(sim_data_file, 'rb') as f:
			self.sim_data = pickle.load(f)
		self.aa_from_trna = self.sim_data.process.transcription.aa_from_trna

		# Get charging parameters and separate to one that will be adjusted or not
		charging_params = get_charging_params(self.sim_data,
			variable_elongation=self.variable_elongation)
		self.adjustable_charging_params = {
			param: value
			for param, value in charging_params.items()
			if param not in CONSTANT_PARAMS
			}
		self.constant_charging_params = {
			param: value
			for param, value in charging_params.items()
			if param in CONSTANT_PARAMS
			}

		# Listeners used
		growth_reader = TableReader(os.path.join(sim_out_dir, 'GrowthLimits'))
		main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))

		# Load data
		self.synthetase_conc = CONC_UNITS * growth_reader.readColumn('synthetase_conc')[1:, :]
		self.uncharged_trna_conc = CONC_UNITS * growth_reader.readColumn('uncharged_trna_conc')[1:, :]
		self.charged_trna_conc = CONC_UNITS * growth_reader.readColumn('charged_trna_conc')[1:, :]
		self.aa_conc = CONC_UNITS * growth_reader.readColumn('aa_conc')[1:, :]
		self.ribosome_conc = CONC_UNITS * growth_reader.readColumn('ribosome_conc')[1:]
		self.fraction_aa_to_elongate = growth_reader.readColumn('fraction_aa_to_elongate')[1:, :]
		self.charging_fraction = growth_reader.readColumn('fraction_trna_charged')[1:, :]

		self.time_step_sizes = main_reader.readColumn('timeStepSec')[1:]
		self.n_time_steps = len(self.time_step_sizes)

	def solve_timestep(self,
			timestep: int,
			synthetase_adjustments: Optional[np.ndarray] = None,
			trna_adjustments: Optional[np.ndarray] = None,
			aa_adjustments: Optional[np.ndarray] = None,
			param_adjustments: Optional[Dict] = None,
			) -> Tuple[np.ndarray, float]:
		"""
		Calculates charging and elongation rate for a given timestep.

		Args:
			timestep: simulation timestep to select data from
			synthetase_adjustments: adjustments to scale synthetase concentrations
				up or down
			trna_adjustments: adjustments to scale tRNA concentrations
				up or down
			aa_adjustments: adjustments to scale amino acid concentrations
				up or down
			param_adjustments: adjustments to charging parameters

		Returns:
			fraction charged of all tRNAs for each amino acid
			ribosome elongation rate
		"""

		n_aas = len(self.sim_data.molecule_groups.amino_acids)

		if synthetase_adjustments is None:
			synthetase_adjustments = np.ones(n_aas)
		if trna_adjustments is None:
			trna_adjustments = np.ones(n_aas)
		if aa_adjustments is None:
			aa_adjustments = np.ones(n_aas)
		if param_adjustments is None:
			param_adjustments = {}

		charging_params = {}
		for param, value in self.adjustable_charging_params.items():
			charging_params[param] = value * param_adjustments.get(param, 1)
		charging_params.update(self.constant_charging_params)

		return calculate_trna_charging(
			self.synthetase_conc[timestep, :] * synthetase_adjustments,
			self.uncharged_trna_conc[timestep, :] * trna_adjustments,
			self.charged_trna_conc[timestep, :] * trna_adjustments,
			self.aa_conc[timestep, :] * aa_adjustments,
			self.ribosome_conc[timestep],
			self.fraction_aa_to_elongate[timestep, :],
			charging_params,
			time_limit=self.time_step_sizes[timestep],
			)

	def validation(self, n_steps: int) -> None:
		"""
		Performs a validation check to makes sure solving the model from
		loaded data matches the objective from the original solution during
		the simulation.

		Args:
			n_steps: number of timesteps to check
				if 0: does not check
				if <0: runs all timepoints from the simulation
		"""

		if n_steps == 0:
			return
		elif n_steps < 0:
			n_steps = self.n_time_steps
		else:
			n_steps = min(self.n_time_steps, n_steps)

		print('Running validation to check output...')
		for timestep in range(n_steps):
			charging_fraction, _ = self.solve_timestep(timestep)
			if np.any(charging_fraction @ self.aa_from_trna != self.charging_fraction[timestep]):
				raise ValueError(f'Charging fraction does not match for time step {timestep}')
		print('All {} timesteps match the results from the whole-cell model.'.format(n_steps))

	def interactive_debug(self):
		"""
		Run an interactive app in a browser to debug charging.
		"""

		app = self.create_app()
		webbrowser.open_new(f'http://127.0.0.1:{PORT}/')
		app.run_server(port=PORT)

	def create_app(self) -> dash.Dash:
		"""
		Create the Dash app to run in the browser for interactive mode.
		"""

		app = dash.Dash()

		param_ids = sorted(self.adjustable_charging_params)
		aa_ids = self.sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		# Slider elements
		slider_options = dict(value=0, min=-2, max=2, step=0.01, marks={i: {'label': 10**i} for i in range(-2, 3)})
		aa_slider_style = {'display': 'grid', 'grid-template-columns': '15% 25% 25% 25%'}
		aa_slider_headers = [html.Div(style=aa_slider_style, children=[
			html.Plaintext(''),
			html.Plaintext('Synthetases', style={'text-align': 'center'}),
			html.Plaintext('tRNA', style={'text-align': 'center'}),
			html.Plaintext('Amino acids', style={'text-align': 'center'}),
			])]
		aa_sliders = [
			html.Div(style=aa_slider_style, children=[
				html.Plaintext(f'{aa[:-3]}:'),
				dcc.Slider(id=f'{aa}-synthetase', **slider_options),
				dcc.Slider(id=f'{aa}-trna', **slider_options),
				dcc.Slider(id=f'{aa}-aa', **slider_options),
				])
			for aa in aa_ids
			]
		param_slider_style = {'display': 'grid', 'grid-template-columns': '30% 70%'}
		param_headers = [html.Div(style=param_slider_style, children=[
			html.Plaintext(''),
			html.Plaintext('Charging parameters', style={'text-align': 'center'}),
			])]
		param_sliders = [
			html.Div(style=param_slider_style, children=[
				html.Plaintext(f'{param}:'),
				dcc.Slider(id=f'{param}-slider', **slider_options),
				])
			for param in param_ids
			]

		sliders = html.Div(style={'display': 'grid', 'grid-template-columns': '70% 30%'}, children=[
			html.Div(children=aa_slider_headers + aa_sliders),
			html.Div(children=param_headers + param_sliders),
			])

		# Slider inputs
		synthetase_inputs = [
			dash.dependencies.Input(f'{aa}-synthetase', 'value')
			for aa in aa_ids
			]
		trna_inputs = [
			dash.dependencies.Input(f'{aa}-trna', 'value')
			for aa in aa_ids
			]
		aa_inputs = [
			dash.dependencies.Input(f'{aa}-aa', 'value')
			for aa in aa_ids
			]
		param_inputs = [
			dash.dependencies.Input(f'{param}-slider', 'value')
			for param in param_ids
			]

		# Page layout
		app.layout = html.Div(children=[
			html.H2('Charging debugger'),
			html.Plaintext('Timestep limits (lower, upper):'),
			dcc.Input(id=LOW_TIMESTEP, type='number', value=0),
			dcc.Input(id=HIGH_TIMESTEP, type='number', value=10),
			sliders,
			dcc.Graph(id=GRAPH_ID, style={'height': '1200px'})
			])

		# Register callback to update plot when selections change
		# First arg for Output/Input selects the page object
		# Second arg for Output/Input sets or gets a kwarg from the dcc function
		@app.callback(
			dash.dependencies.Output(GRAPH_ID, 'figure'),
			[
				dash.dependencies.Input(LOW_TIMESTEP, 'value'),
				dash.dependencies.Input(HIGH_TIMESTEP, 'value'),
				*synthetase_inputs,
				*trna_inputs,
				*aa_inputs,
				*param_inputs,
			])
		def update_graph(init_t: int, final_t: int, *param_inputs: float) -> Dict:
			"""
			Update the plot based on selection changes.

			Returns:
				plotly figure dict
			"""

			synthetase_adjustments = 10**np.array(param_inputs[:n_aas])
			trna_adjustments = 10**np.array(param_inputs[n_aas:2*n_aas])
			aa_adjustments = 10**np.array(param_inputs[2*n_aas:3*n_aas])
			param_adjustments = {param: 10**value for param, value in zip(param_ids, param_inputs[3*n_aas:])}

			v_rib = []
			f_charged = []
			t = np.arange(init_t, final_t)
			for timestep in t:
				f, v = self.solve_timestep(
					timestep,
					synthetase_adjustments=synthetase_adjustments,
					trna_adjustments=trna_adjustments,
					aa_adjustments=aa_adjustments,
					param_adjustments=param_adjustments
					)
				v_rib.append(v / self.ribosome_conc[timestep].asNumber(CONC_UNITS))
				f_charged.append(f)

			f_charged = np.array(f_charged).T

			fig = plotly.subplots.make_subplots(rows=2, cols=1)
			fig.append_trace(go.Scatter(x=t, y=v_rib, name='Elongation rate'), row=1, col=1)
			for f, aa in zip(f_charged, aa_ids):
				fig.append_trace(go.Scatter(x=t, y=f, name=aa), row=2, col=1)

			fig.update_xaxes(title_text='Timestep', row=2, col=1)
			fig.update_yaxes(title_text='Elongation rate (AA/s)', row=1, col=1)
			fig.update_yaxes(title_text='Fraction charged', row=2, col=1)

			return fig

		return app

	def run(self, args: argparse.Namespace) -> None:
		self.variable_elongation = args.variable_elongation_translation
		self.load_data(args.sim_data_file, args.sim_out_dir)
		self.validation(args.validation)
		if args.interactive:
			self.interactive_debug()


if __name__ == '__main__':
	analysis = ChargingDebug()
	analysis.cli()
