#! /usr/bin/env python

"""
Explore simulation data in an interactive manner with the ability to select
datasets and graph types.

TODO:
	- update names displayed on page
	- default column option? - start with time as x
	- add plotting options - log scale
	- add reducing options - mean across samples
"""

import argparse
import os
import re
import webbrowser

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

from wholecell.io.tablereader import DoesNotExistError, TableReader
from wholecell.utils import scriptBase
import wholecell.utils.filepath as fp


PORT = 8050

# Options from drop down menu
PLOT_OPTIONS = {
	'line': {'function': go.Scatter},
	'scatter': {'function': go.Scatter, 'plot_options': {'mode': 'markers'}},
	'bar': {'function': go.Bar, 'layout_options': {'barmode': 'stack'}},
	}

# Object IDs
GRAPH_ID = 'test-plot'
PLOT_SELECTION = 'plot'
X_DATA_SELECTION_ID = 'x-data'
Y_DATA_SELECTION_ID = 'y-data'
SEPARATOR = '<>'
VALUE_JOIN = f'{{}}{SEPARATOR}{{}}'


def get_vals(d, k):
	if k is None:
		return d

	if type(k) == str:
		k = k.split(SEPARATOR)

	if len(k) > 1:
		return get_vals(d[k[0]], k[1:])
	else:
		return d.get(k[0])

def load_listener(input):
	split = input.split(SEPARATOR)
	path = os.path.join(*split[:-2], 'simOut')
	listener = split[-2]
	column = split[-1]

	reader = TableReader(os.path.join(path, listener))
	data = reader.readColumn2D(column)

	# Read subcolumn attributes
	try:
		subcolumns = reader.readAttribute('subcolumns')
	except DoesNotExistError:
		subcolumns = []
	if column in subcolumns:
		labels = reader.readAttribute(subcolumns[column])
	else:
		labels = list(range(data.shape[1]))

	return data, labels

def create_app(data_structure):
	def data_selection(app, data_structure, id_, multi=False):
		def add_children(children, base_id, parent_id, default, data_structure, count=0):
			vals = get_vals(data_structure, default)
			if vals is not None:
				sub_id = f'{base_id}{count}'
				children.append(dcc.Dropdown(
					id=sub_id,
					multi=multi,
					))

				@app.callback(
					[
						dash.dependencies.Output(sub_id, 'options'),
						dash.dependencies.Output(sub_id, 'value'),
						],
					[dash.dependencies.Input(parent_id, 'value')]
					)
				def update(path):
					vals = sorted(get_vals(data_structure, path))

					options = [{
						'label': v,
						'value': VALUE_JOIN.format(path, v),
						} for v in vals]
					value = VALUE_JOIN.format(path, vals[0])

					return options, value

				default = VALUE_JOIN.format(default, next(iter(vals)))
				return add_children(children, base_id, sub_id, default, data_structure, count=count+1)
			else:
				return children, count

		default_top_level = next(iter(data_structure))

		children = [
			html.H2(id_),
			dcc.Dropdown(
				id=id_,
				options=[{
					'label': os.path.basename(d),
					'value': d,
					} for d in data_structure],
				value=default_top_level,
				)
			]

		children, n_added = add_children(children, id_, id_, default_top_level, data_structure)
		input = dash.dependencies.Input(f'{id_}{n_added-1}', 'value')

		div = html.Div(children=children)

		return div, input

	# Create webpage layout
	app = dash.Dash()
	x_div, x_input = data_selection(app, data_structure, X_DATA_SELECTION_ID)
	y_div, y_input = data_selection(app, data_structure, Y_DATA_SELECTION_ID, multi=False)  # TODO: get multi selection working
	app.layout = html.Div(children=[
		html.H1('Interactive test'),
		html.Div(children=[
			html.H2('Plot selection:'),
			dcc.Dropdown(
				id=PLOT_SELECTION,
				options=[{
					'label': o,  # display name
					'value': o,  # value passed through callback, must be str
					} for o in PLOT_OPTIONS],
				value=next(iter(PLOT_OPTIONS)),
				),
			x_div,
			y_div,
			]),
		dcc.Graph(id=GRAPH_ID),
		])

	# First arg for Output/Input selects the page object
	# Second arg for Output/Input sets or gets a kwarg from the dcc function
	@app.callback(
		dash.dependencies.Output(GRAPH_ID, 'figure'),
		[dash.dependencies.Input(PLOT_SELECTION, 'value'), x_input, y_input])
	def update_graph(plot_id, x_input, y_input):
		x_data, x_labels = load_listener(x_input)
		y_data, y_labels = load_listener(y_input)

		plot = PLOT_OPTIONS[plot_id]
		plot_options = plot.get('plot_options', {})
		layout_options = plot.get('layout_options', {})
		traces = [
			plot['function'](x=x_data[:, 0], y=y_data[:, idx], name=col, **plot_options)
			for idx, col in enumerate(y_labels)
			]

		# Dict used to update 'figure' for dcc.Graph object GRAPH_ID
		return {
			'data': traces,
			'layout': go.Layout(
				title=plot_id,
				xaxis_title=x_input.split(SEPARATOR)[-1],
				yaxis_title=y_input.split(SEPARATOR)[-1],
				**layout_options),
			}

	return app

class AnalysisInteractive(scriptBase.ScriptBase):
	def define_parameters(self, parser):
		super().define_parameters(parser)

		self.define_parameter_sim_dir(parser, default=os.path.join(fp.ROOT_PATH, 'out'))

	def parse_data_structure(self, path: str):
		experiments = {}
		variant_dirs = [v[0] for v in self.list_variant_dirs(path)]
		if len(variant_dirs):
			experiments[path] = {d: {} for d in variant_dirs}
		else:
			for dir in os.listdir(path):
				sim_path = os.path.join(path, dir)
				if not os.path.isdir(sim_path):
					continue

				variant_dirs = [v[0] for v in self.list_variant_dirs(sim_path)]
				if len(variant_dirs):
					experiments[sim_path] = {d: {} for d in variant_dirs}

		seed_regex = re.compile('[0-9]{6}')
		generation_regex = re.compile('generation_[0-9]{6}')
		daughter_regex = re.compile('[0-9]{6}')
		for base, variants in experiments.items():
			for variant, seed_dict in variants.items():
				variant_dir = os.path.join(base, variant)
				for seed in os.listdir(variant_dir):
					if seed_regex.match(seed):
						gen_dict = seed_dict.get(seed, {})
						seed_dir = os.path.join(variant_dir, seed)
						for gen in os.listdir(seed_dir):
							if generation_regex.match(gen):
								daughter_dict = gen_dict.get(gen, {})
								gen_dir = os.path.join(seed_dir, gen)
								for daughter in os.listdir(gen_dir):
									if daughter_regex.match(daughter):
										listener_dict = daughter_dict.get(daughter, {})
										sim_out_dir = os.path.join(gen_dir, daughter, 'simOut')
										for listener in os.listdir(sim_out_dir):
											column_dict = listener_dict.get(listener, {})
											listener_dir = os.path.join(sim_out_dir, listener)
											if os.path.isdir(listener_dir):
												for column in os.listdir(listener_dir):
													if '.' not in column:  # TODO: handle .cPickle columns
														column_dict[column] = None
												listener_dict[listener] = column_dict
										daughter_dict[daughter] = listener_dict
								gen_dict[gen] = daughter_dict
						seed_dict[seed] = gen_dict

		return experiments

	def run(self, args: argparse.Namespace) -> None:
		data_structure = self.parse_data_structure(args.sim_dir)
		app = create_app(data_structure)

		# Serve interactive page (may take a second to load, reload if necessary)
		webbrowser.open_new(f'http://127.0.0.1:{PORT}/')
		app.run_server(port=PORT)


if __name__ == '__main__':
	analysis = AnalysisInteractive()
	analysis.cli()
