"""
Plots proteomics comparison plots for two sets of simulations.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os
from scipy.stats import pearsonr

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
import plotly.express as px
import plotly.graph_objects as go
from wholecell.utils.protein_counts import get_simulated_validation_counts




class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):

	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):

		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		if ap1.n_generation <= 2 or ap2.n_generation <= 2:
			print('Skipping analysis -- not enough sims run.')
			return

		def read_sim_protein_counts(ap, monomer_ids=sim_data1.process.translation.monomer_data["id"]):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))


			monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True).mean(axis=0)

			return monomer_counts

		monomer_counts2 = read_sim_protein_counts(ap2)  # Sim2 (new sim)
		monomer_counts1_raw = read_sim_protein_counts(ap1)  # Sim1 (reference)
		hi = 5

		# Get the monomer IDs:
		monomer_IDs2_raw = sim_data2.process.translation.monomer_data["id"].tolist()
		monomer_IDs1_raw = sim_data1.process.translation.monomer_data["id"].tolist()

		# Strip compartment tags (last 3 characters: [X])
		def strip_compartment(mol_id):
			"""Remove compartment tag [X] from molecule ID"""
			if mol_id.endswith(']') and mol_id[-3] == '[':
				return mol_id[:-3]
			return mol_id

		# Create mapping dicts for looking up by base name
		monomer_base_to_full_sim1 = {strip_compartment(mol_id): mol_id for mol_id in
									 monomer_IDs1_raw}
		monomer_base_to_full_sim2 = {strip_compartment(mol_id): mol_id for mol_id in
									 monomer_IDs2_raw}

		# Use sim2's FULL IDs (with compartments) as reference
		monomer_IDs = monomer_IDs2_raw  # ← This is a list
		monomer_dict = {mol: i for i, mol in enumerate(monomer_IDs)}

		# Align sim1 counts to sim2's monomer order (matching by base name)
		monomer_dict1_raw = {mol_id: monomer_counts1_raw[i] for i, mol_id in
							 enumerate(monomer_IDs1_raw)}

		monomer_counts1 = []
		for mol_id2 in monomer_IDs2_raw:
			base_name = strip_compartment(mol_id2)
			# Find corresponding sim1 molecule (might have different compartment tag)
			mol_id1 = monomer_base_to_full_sim1.get(base_name)
			if mol_id1 and mol_id1 in monomer_dict1_raw:
				monomer_counts1.append(monomer_dict1_raw[mol_id1])
			else:
				monomer_counts1.append(0.0)

		monomer_counts1 = np.array(monomer_counts1)

		# Warn about differences
		missing_base_names = set(monomer_base_to_full_sim2.keys()) - set(
			monomer_base_to_full_sim1.keys())
		if missing_base_names:
			print(f"Note: {len(missing_base_names)} monomers in sim2 not in sim1 (new in model)")
			print(f"  Examples: {list(missing_base_names)[:5]}")

		# Read validation protein counts
		val_monomer_ids = validation_data1.protein.schmidt2015Data["monomerId"]
		val_monomer_counts = validation_data1.protein.schmidt2015Data["glucoseCounts"]

		# Find where validation monomer ids are in sim2 monomer ids
		val_monomer_indexes = [monomer_dict.get(mol_id, None) for mol_id in val_monomer_ids]
		val_monomer_counts_sim2 = np.array(
			[monomer_counts2[idx] if idx is not None else 0.0 for idx in val_monomer_indexes])

		# Helper function for molecule categorization
		def get_monomer_indexes(keys):
			indices = []
			for x in keys:
				if x in monomer_dict:
					indices.append(monomer_dict[x])
				else:
					# Skip molecules not in monomer_IDs (e.g., modified forms, metabolites)
					pass
			return np.array(indices)

		# STEP 1: Get the monomer IDs for each group type:

		# COMPLEXATION COMPLEXES:
		complexation_molecule_ids = sim_data2.process.complexation.molecule_names
		complexation_complex_ids = sim_data2.process.complexation.ids_complexes
		complexation_smm = sim_data2.process.complexation.stoich_matrix_monomers()
		complexation_monomers = {}
		complexes_to_monomers = {}
		for id in complexation_molecule_ids:
			if id in monomer_dict:
				complexation_monomers[id] = monomer_dict[id]
			elif id in complexation_complex_ids:
				# Find constituent monomers of this complex
				complex_idx = complexation_complex_ids.index(id)
				monomer_idxs = np.where(complexation_smm[:, complex_idx] < 0)[0]
				monomers = []
				for idx in monomer_idxs:
					monomer_id = complexation_molecule_ids[idx]
					if monomer_id in monomer_dict:
						complexation_monomers[monomer_id] = monomer_dict[monomer_id]
						monomers.append(monomer_id)
				complexes_to_monomers[id] = monomers

		# EQUILIBRIUM COMPLEXES:
		equilibrium_molecule_ids = sim_data2.process.equilibrium.molecule_names
		equilibrium_complex_ids = sim_data2.process.equilibrium.ids_complexes
		eq_smm = sim_data2.process.equilibrium.stoich_matrix_monomers()
		eq_complexes_to_monomers = {}
		equilibrium_monomers = {}
		for id in equilibrium_molecule_ids:
			if id in monomer_dict:
				equilibrium_monomers[id] = monomer_dict[id]
			elif id in complexation_complex_ids:
				# Use complexation monomers
				if id in complexes_to_monomers:
					for monomer_id in complexes_to_monomers[id]:
						if monomer_id in monomer_dict:
							equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
			elif id in equilibrium_complex_ids:
				# Find constituent molecules of this equilibrium complex
				eq_complex_idx = equilibrium_complex_ids.index(id)
				molecule_idxs = np.where(eq_smm[:, eq_complex_idx] < 0)[0]
				monomers = []
				for idx in molecule_idxs:
					molecule_name = equilibrium_molecule_ids[idx]
					if molecule_name in monomer_dict:
						monomers.append(molecule_name)
						equilibrium_monomers[molecule_name] = monomer_dict[molecule_name]
					elif molecule_name in complexation_complex_ids:
						# Expand complexation complexes to monomers
						if molecule_name in complexes_to_monomers:
							for monomer_id in complexes_to_monomers[molecule_name]:
								if monomer_id in monomer_dict:
									monomers.append(monomer_id)
									equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
				eq_complexes_to_monomers[id] = monomers

		# TCS COMPLEXES:
		two_component_system_molecule_ids = (
			list(sim_data2.process.two_component_system.modified_molecules))
		two_component_system_complex_ids = (
			list(sim_data2.process.two_component_system.complex_to_monomer.keys()))
		tcs_complex_to_monomer_dict = sim_data2.process.two_component_system.complex_to_monomer
		two_component_system_monomers = {}
		tcs_complexes_to_monomers = {}
		for id in two_component_system_molecule_ids:
			monomers = []
			if id in monomer_IDs:
				# This has all the base monomers in it already!
				two_component_system_monomers[id] = get_monomer_indexes([id])[0]
			elif id in complexation_complex_ids:
				# If the molecule is a complexation complex, find its constituent monomers and add them to the list
				monomer_ids = complexes_to_monomers[id]
				tcs_complexes_to_monomers[id] = monomer_ids
			elif id in equilibrium_complex_ids:
				# If the molecule is an equilibrium complex, find its constituent monomers and add them to the list
				monomer_ids = eq_complexes_to_monomers[id]
				tcs_complexes_to_monomers[id] = monomer_ids
			elif id in two_component_system_complex_ids:
				# If the molecule is a two component system complex, find its constituent monomers and add them to the list
				monomer_ids = tcs_complex_to_monomer_dict[id]
				tcs_complexes_to_monomers[id] = monomer_ids.keys()
			for monomer_id in monomers:
				if monomer_id in monomer_IDs:
					two_component_system_monomers[monomer_id] = get_monomer_indexes([monomer_id])[0]


		# RIBOSOMES:
		ribosome_50s_subunits = sim_data2.process.complexation.get_monomers(
			sim_data1.molecule_ids.s50_full_complex)
		ribosome_30s_subunits = sim_data2.process.complexation.get_monomers(
			sim_data1.molecule_ids.s30_full_complex)
		ribosome_subunits = (ribosome_50s_subunits["subunitIds"].tolist() +
								ribosome_30s_subunits["subunitIds"].tolist())
		ribosome_monomers = {}
		for subunit in ribosome_subunits:
			if subunit in monomer_IDs:
				ribosome_monomers[subunit] = get_monomer_indexes([subunit])[0]

		# RNAPs:
		rnap_subunits = sim_data2.process.complexation.get_monomers(
			sim_data1.molecule_ids.full_RNAP)
		rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()
		rnap_monomers = {}
		for subunit in rnap_subunit_ids:
			if subunit in monomer_IDs:
				rnap_monomers[subunit] = get_monomer_indexes([subunit])[0]
			if subunit in complexation_complex_ids:
				monomer_ids = complexes_to_monomers[subunit]
				for id in monomer_ids:
					if id in monomer_IDs:
						rnap_monomers[id] = get_monomer_indexes([id])[0]

		# REPLISOMES:
		replisome_trimer_subunits = sim_data2.molecule_groups.replisome_trimer_subunits
		replisome_monomer_subunits = sim_data2.molecule_groups.replisome_monomer_subunits
		replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits
		replisome_monomers = {}
		for subunit in replisome_subunit_ids:
			if subunit in monomer_IDs:
				replisome_monomers[subunit] = get_monomer_indexes([subunit])[0]
			elif subunit in complexation_complex_ids:
				monomer_ids = complexes_to_monomers[subunit]
				for id in monomer_ids:
					if id in monomer_IDs:
						replisome_monomers[id] = get_monomer_indexes([id])[0]


		# TF MONOMERS
		tfs = sim_data2.process.transcription_regulation.tf_ids
		tf_subunit_ids = [tf_id + f'[{sim_data1.getter.get_compartment(tf_id)[0]}]'
						  for tf_id in tfs]
		tf_monomers = {}
		for subunit in tf_subunit_ids:
			if subunit in monomer_IDs:
				tf_monomers[subunit] = get_monomer_indexes([subunit])[0]
			elif subunit in complexation_complex_ids:
				monomer_ids = complexes_to_monomers[subunit]
				for id in monomer_ids:
					if id in monomer_IDs:
						tf_monomers[id] = get_monomer_indexes([id])[0]

			elif subunit in equilibrium_complex_ids:
				monomer_ids = eq_complexes_to_monomers[subunit]
				for id in monomer_ids:
					if id in monomer_IDs:
						tf_monomers[id] = get_monomer_indexes([id])[0]

			elif subunit in two_component_system_complex_ids:
				monomer_ids = tcs_complexes_to_monomers[subunit]
				for id in monomer_ids:
					if id in monomer_IDs:
						tf_monomers[id] = get_monomer_indexes([id])[0]


		# OK now plot
		def classify_molecule_types_multi(monomer_id, monomer_idx):
			"""Classify molecule into ALL categories it belongs to"""
			categories = []

			if monomer_id in ribosome_monomers and ribosome_monomers[monomer_id] == monomer_idx:
				categories.append('Ribosome')
			if monomer_id in rnap_monomers and rnap_monomers[monomer_id] == monomer_idx:
				categories.append('RNAP')
			if monomer_id in replisome_monomers and replisome_monomers[monomer_id] == monomer_idx:
				categories.append('Replisome')
			if monomer_id in tf_monomers and tf_monomers[monomer_id] == monomer_idx:
				categories.append('transcription factor')
			if monomer_id in complexation_monomers and complexation_monomers[
				monomer_id] == monomer_idx:
				categories.append('Complexation complex')
			if monomer_id in equilibrium_monomers and equilibrium_monomers[
				monomer_id] == monomer_idx:
				categories.append('Equilibrium complex')
			if monomer_id in two_component_system_monomers and two_component_system_monomers[
				monomer_id] == monomer_idx:
				categories.append('two component system')

			if not categories:
				return 'Monomer only'
			elif len(categories) == 1:
				return categories[0] + ' subunit'
			else:
				# Create compound category name
				return ' + '.join(sorted(categories))

		# Create dataframe
		import pandas as pd

		data = []
		for i, monomer_id in enumerate(monomer_IDs):
			mol_type = classify_molecule_types_multi(monomer_id, i)
			data.append({
				'monomer_id': monomer_id,
				'monomer_idx': i,
				'sim1_count': monomer_counts1[i],
				'sim2_count': monomer_counts2[i],
				'sim1_log': np.log10(monomer_counts1[i] + 1),
				'sim2_log': np.log10(monomer_counts2[i] + 1),
				'molecule_type': mol_type
			})

		df = pd.DataFrame(data)

		# Get unique molecule types
		unique_types = df['molecule_type'].unique()

		# Create color palette
		import plotly.express as px
		import plotly.graph_objects as go

		colors = px.colors.qualitative.Plotly + px.colors.qualitative.Set2 + px.colors.qualitative.Pastel

		# Assign colors and styles
		molecule_type_styles = {}
		priority = 1

		# Compound types (multiple categories) - most interesting!
		compound_types = sorted([t for t in unique_types if '+' in t])
		for i, mol_type in enumerate(compound_types):
			molecule_type_styles[mol_type] = {
				'color': colors[i % len(colors)],
				'size': 8,
				'symbol': 'star',
				'priority': priority,
				'opacity': .8,
			}
			priority += 1

		# Single types - with complexation getting lower priority
		single_types = sorted([t for t in unique_types if '+' not in t and t != 'Monomer only'])

		# Separate out complexation
		high_priority_types = [t for t in single_types if t != 'Complexation complex subunit']
		complexation_type = [t for t in single_types if t == 'Complexation complex subunit']

		# Plot high priority single types first
		for i, mol_type in enumerate(high_priority_types):
			molecule_type_styles[mol_type] = {
				'color': colors[(i + len(compound_types)) % len(colors)],
				'size': 5,
				'symbol': 'circle',
				'priority': priority,
				'opacity': .9,
				'line_width': 0.7,
			}
			priority += 1

		# Then complexation (lower priority, plotted earlier/background)
		for i, mol_type in enumerate(complexation_type):
			molecule_type_styles[mol_type] = {
				'color': colors[(len(high_priority_types) + len(compound_types)) % len(colors)],
				'size': 3,  # Slightly smaller too
				'symbol': 'circle',
				'priority': priority,
				'opacity': 0.6
			}
			priority += 1

		# Plain monomers (background)
		if 'Monomer only' in unique_types:
			molecule_type_styles['Monomer only'] = {
				'color': 'lightgray',
				'size': 3,
				'symbol': 'circle',
				'priority': 999,
				'opacity': 0.4
			}

		# Create figure
		fig = go.Figure()

		# Plot each molecule type (reverse priority so important ones on top)
		for mol_type, style in sorted(molecule_type_styles.items(),
									  key=lambda x: x[1]['priority'], reverse=True):
			subset = df[df['molecule_type'] == mol_type]

			if len(subset) == 0:
				continue

			x = subset['sim1_log']
			y = subset['sim2_log']
			hovertext = subset.apply(
				lambda row: f"Monomer ID: {row['monomer_id']}<br>"
							f"Molecule Type: {row['molecule_type']}<br>"
							f"Sim1 Count: {row['sim1_count']:.1f}<br>"
							f"Sim2 Count: {row['sim2_count']:.1f}<br>",
				axis=1
			)

			# Add special label for compound types
			name = f"{mol_type} (n={len(subset)})" if '+' in mol_type else f"{mol_type} (n={len(subset)})"

			fig.add_trace(go.Scatter(
				x=x, y=y,
				hovertext=hovertext,
				mode='markers',
				name=name,
				marker=dict(
					color=style['color'],
					size=style['size'],
					opacity=style['opacity'],
					symbol=style.get('symbol', 'circle'),
					line=dict(width=1, color='black') if '+' in mol_type else dict(width=0)
				)
			))

		# Add y=x line
		max_val = max(df['sim1_log'].max(), df['sim2_log'].max())
		fig.add_trace(go.Scatter(
			x=[0, max_val], y=[0, max_val],
			mode="lines",
			line=dict(color="black", dash="dash"),
			opacity=0.2,
			name="y=x",
			showlegend=True
		))

		sim1 = reference_sim_dir.split('out/')[-1]
		sim2 = input_sim_dir.split('out/')[-1]
		# Update layout
		fig.update_layout(
			title=f"Protein count comparison: Sim1 vs Sim2<br>"
				  f"Sim1 ID: {sim1} | Sim2 ID: {sim2}",
			xaxis_title=f"log10(Sim1 monomer counts + 1)",
			yaxis_title=f"log10(Sim2 monomer counts + 1)",
			autosize=False,
			width=1100,
			height=600,
			showlegend=True,
			plot_bgcolor='white',
			legend_title="Monomer Type (⭐ = multiple categories)"
		)
		fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
		fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)



		# Also save as HTML
		out_name = f"{plotOutFileName}_sim1_{sim1}_sim2_{sim2}.html"
		html_path = os.path.join(plotOutDir, out_name)
		fig.write_html(html_path)

		# Print summary of compound types
		print("\nMolecules belonging to multiple categories:")
		for mol_type in compound_types:
			count = len(df[df['molecule_type'] == mol_type])
			molecules = df[df['molecule_type'] == mol_type]['monomer_id'].tolist()
			print(f"\n{mol_type} ({count} molecules):")
			for mol in molecules[:5]:  # Print first 5
				print(f"  - {mol}")
			if len(molecules) > 5:
				print(f"  ... and {len(molecules) - 5} more")

		# hi
		import matplotlib.pyplot as plt
		import matplotlib.patches as mpatches
		import matplotlib.pyplot as plt
		import matplotlib.patches as mpatches
		import matplotlib.colors as mcolors

		def convert_plotly_color_to_matplotlib(color_str):
			"""Convert Plotly color format to matplotlib-compatible format"""
			if color_str.startswith('rgb('):
				# Extract RGB values from 'rgb(r, g, b)'
				rgb_str = color_str[4:-1]  # Remove 'rgb(' and ')'
				r, g, b = [int(x.strip()) for x in rgb_str.split(',')]
				return (r / 255, g / 255, b / 255)  # Matplotlib uses 0-1 range
			elif color_str.startswith('#'):
				return color_str  # Hex colors work as-is
			else:
				# Try to use as named color
				try:
					return mcolors.to_rgb(color_str)
				except:
					return 'blue'  # Fallback

		# Define the 7 base categories we want to plot
		base_categories = [
			'Ribosome',
			'RNAP',
			'Replisome',
			'transcription factor',
			'Complexation complex',
			'Equilibrium complex',
			'two component system'
		]

		# Function to map molecule types to base categories
		def get_base_category(mol_type):
			"""Extract base category from molecule type string"""
			if mol_type == 'Monomer only':
				return None

			# Handle compound types - assign to first category mentioned
			for base_cat in base_categories:
				if base_cat in mol_type:
					return base_cat

			return None

		# Add validation data to dataframe
		val_dict = {mol_id: val_monomer_counts[i] for i, mol_id in enumerate(val_monomer_ids)}
		df['validation_count'] = df['monomer_id'].map(val_dict).fillna(0.0)
		df['validation_log'] = np.log10(df['validation_count'] + 1)

		# Filter to only molecules with validation data
		df_with_val = df[df['validation_count'] > 0].copy()

		# Add base category to dataframe
		df_with_val['base_category'] = df_with_val['molecule_type'].apply(get_base_category)

		# Filter out None (Monomer only)
		df_with_val = df_with_val[df_with_val['base_category'].notna()].copy()

		print(f"\n=== BASE CATEGORIES ===")
		for base_cat in base_categories:
			count = len(df_with_val[df_with_val['base_category'] == base_cat])
			print(f"{base_cat}: {count} proteins")
		print("=======================\n")

		# Add base category to dataframe
		df_with_val['base_category'] = df_with_val['molecule_type'].apply(get_base_category)

		# Filter out None (Monomer only)
		df_with_val = df_with_val[df_with_val['base_category'].notna()].copy()

		print(f"\n=== BASE CATEGORIES ===")
		for base_cat in base_categories:
			count = len(df_with_val[df_with_val['base_category'] == base_cat])
			print(f"{base_cat}: {count} proteins")
		print("=======================\n")

		# Create 2x4 subplot grid (7 categories + 1 combined)
		fig, axes = plt.subplots(2, 4, figsize=(20, 10))
		axes = axes.flatten()

		# Plot each base category in its own subplot
		for plot_idx, base_cat in enumerate(base_categories):
			ax = axes[plot_idx]

			subset = df_with_val[df_with_val['base_category'] == base_cat]

			if len(subset) == 0:
				ax.text(0.5, 0.5, f'No validation data matches\nfor {base_cat} subunits',
						ha='center', va='center', transform=ax.transAxes)
				ax.set_title(f"{base_cat} (n=0)")
				continue

			# Get ONE color for this entire category
			first_mol_type = subset['molecule_type'].iloc[0]
			style = molecule_type_styles.get(first_mol_type, {})
			category_color = convert_plotly_color_to_matplotlib(style.get('color', 'blue'))

			# Plot all proteins in background (gray)
			ax.scatter(df_with_val['validation_log'], df_with_val['sim2_log'],
					   color='lightgray', s=5, alpha=0.1, zorder=1)

			# Draw arrows from sim1 → sim2 (NOT to validation)
			# Skip arrows for Complexation (too many)
			if base_cat != 'Complexation complex':
				for idx, row in subset.iterrows():
					x_val = row['validation_log']  # Validation x-position
					x_sim1 = row['sim1_log']  # Sim1 x-position (arrow start)
					x_sim2 = row['validation_log']  # Same as validation for x
					y_sim1 = row['sim1_log']  # Sim1 y-position (arrow start)
					y_sim2 = row['sim2_log']  # Sim2 y-position (arrow end)

					# Arrow from (x_val, y_sim1) to (x_val, y_sim2)
					# Shows vertical movement from sim1 to sim2 at the validation x-position
					ax.annotate('',
								xy=(x_val, y_sim2),  # Arrow head at sim2
								xytext=(x_val, y_sim1),  # Arrow tail at sim1
								arrowprops=dict(
									arrowstyle='->',
									color=category_color,  # Same color for all
									lw=1,
									alpha=0.5,
									shrinkA=0, shrinkB=0
								),
								zorder=2)

				# Plot sim1 starting points (small circles)
				ax.scatter(subset['validation_log'], subset['sim1_log'],
						   color=category_color,  # Same color
						   s=10,
						   marker='o',
						   edgecolors='gray',
						   linewidths=0.5,
						   alpha=0.4,
						   zorder=2,
						   label='Sim1')

			# Plot sim2 points (final positions) - all same color
			ax.scatter(subset['validation_log'], subset['sim2_log'],
					   color=category_color,  # Same color for all
					   s=20,
					   marker='o',
					   edgecolors=category_color,
					   linewidths=0.9,
					   alpha=0.5,
					   label='Sim2',
					   zorder=3)

			# Add y=x line
			max_val = max(df_with_val[['validation_log', 'sim2_log', 'sim1_log']].max())
			ax.plot([0, max_val], [0, max_val], 'k--',
					linewidth=0.5, alpha=0.3, zorder=1, label='y=x')

			# Styling
			ax.set_title(f"{base_cat}\n(n={len(subset)})", fontsize=10)
			ax.set_xlabel('log10(Validation counts + 1)', fontsize=8)
			ax.set_ylabel('log10(Simulated counts + 1)', fontsize=8)
			ax.legend(fontsize=7, loc='upper left')
			ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
			ax.set_facecolor('white')

		# 8th subplot: All categories combined (NO arrows, too cluttered)
		ax = axes[7]

		# Plot background
		ax.scatter(df_with_val['validation_log'], df_with_val['sim2_log'],
				   color='lightgray', s=5, alpha=0.2, zorder=1)

		# Plot each base category with its own color
		for base_cat in base_categories:
			subset = df_with_val[df_with_val['base_category'] == base_cat]

			if len(subset) == 0:
				continue

			# Use first molecule type's color for this category
			first_mol_type = subset['molecule_type'].iloc[0]
			style = molecule_type_styles.get(first_mol_type, {})
			color = convert_plotly_color_to_matplotlib(style.get('color', 'blue'))

			ax.scatter(subset['validation_log'], subset['sim2_log'],
					   color=color,
					   s=10,
					   marker='o',
					   alpha=0.7,
					   label=f"{base_cat} (n={len(subset)})",
					   zorder=3)

		# Calculate R²
		from scipy.stats import pearsonr
		r_val, p_val = pearsonr(df_with_val['validation_log'], df_with_val['sim2_log'])
		r2 = r_val ** 2

		# Add text box in bottom right
		ax.text(0.95, 0.05, f'$R^2$ = {r2:.3f}',
				transform=ax.transAxes,
				fontsize=10,
				verticalalignment='bottom',
				horizontalalignment='right',
				bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))

		# Add y=x line
		max_val = max(df_with_val[['validation_log', 'sim2_log']].max())
		ax.plot([0, max_val], [0, max_val], 'k--',
				linewidth=0.5, alpha=0.3, zorder=1)

		# Styling
		ax.set_title('All Subunit Categories Combined', fontsize=10, fontweight='bold')
		ax.set_xlabel('log10(Validation counts + 1)', fontsize=8)
		ax.set_ylabel('log10(Sim2 counts + 1)', fontsize=8)
		ax.legend(fontsize=7, loc='upper left', ncol=1)
		ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
		ax.set_facecolor('white')

		# Overall title
		fig.suptitle(f'Protein Counts: Validation (x-axis) vs Simulated (y-axis)\n'
					 f'Arrows show Sim1 ({sim1}) → Sim2 ({sim2}) change ',
					 fontsize=14, y=0.995)

		plt.tight_layout(rect=[0, 0, 1, 0.99])

		# Save
		validation_plot_filename = plotOutFileName + '_validation_by_subunit_type' + f"_sim1_{sim1}_sim2_{sim2}"
		exportFigure(plt, plotOutDir, validation_plot_filename, metadata)

		print(f"\nSaved validation subplot grid with {len(df_with_val)} proteins")


	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
