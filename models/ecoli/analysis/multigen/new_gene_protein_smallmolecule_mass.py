import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import inf
import pandas as pd

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import determine_new_gene_ids_and_indices
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns, read_bulk_molecule_counts)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def get_monomer_name_with_compartment(self, simDataFile, monomer):
		"""
		Obtain the compartment of the monomer from the simDataFile.
		"""
		sim_data = self.read_pickle_file(simDataFile)
		if not monomer.endswith("]"):
			if sim_data.getter.is_valid_molecule(monomer):
				edited_monomer_name = monomer + "[" + sim_data.getter.get_compartment(monomer)[
					0] + "]"
			else:
				print(f"Monomer {monomer} not found in simDataFile.")
				edited_monomer_name = None
		else:
			# pressumably it already has the compartment tag, so just return it:
			print(f"Monomer {monomer} already has a compartment tag.")
			edited_monomer_name = monomer

		return edited_monomer_name
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile,
				metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		dt = read_stacked_columns(
			cell_paths, 'Main', 'time',
			fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
		num_time_steps = read_stacked_columns(
			cell_paths, 'Main', 'time',
			fun=lambda x: len(x)).squeeze()
		gen_labels = np.repeat(np.arange(len(dt)), num_time_steps)
		unique_gen_labels = np.unique(gen_labels)
		gen_start_index = np.array(
			[gen_labels.tolist().index(i) for i in unique_gen_labels])
		gen_end_index = np.concatenate((
			np.array(gen_start_index[1:] - 1), np.array([len(gen_labels) - 1])))

		small_mol_ids = sim_data.process.metabolism.concentration_updates._all_metabolite_ids

		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_reader.readAttribute('monomerIds')
		monomer_mw = sim_data.getter.get_masses(monomer_ids).asNumber(units.fg / units.count)

		monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', 'monomerCounts',
			remove_first=True, ignore_exception=True)
		dry_masses = read_stacked_columns(
			cell_paths, 'Mass', 'dryMass',
			remove_first=True, ignore_exception=True)

		bulk_reader = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
		bulk_molecule_names = bulk_reader.readAttribute("objectNames")
		for name in small_mol_ids:
			small_molecule_name_now = self.get_monomer_name_with_compartment(simDataFile,
																name)
			name = [small_molecule_name_now]

		small_mol_counts = read_stacked_bulk_molecules(cell_paths, small_mol_ids)
		small_mol_mw =sim_data.getter.get_masses(small_mol_ids).asNumber(units.fg / units.count)


		table_rows = []
		table_rows2 = []


		for g, (start_idx, end_idx) in enumerate(zip(gen_start_index, gen_end_index)):
			final_monomer_counts = monomer_counts[end_idx, :]
			final_monomer_masses = final_monomer_counts * monomer_mw
			final_sm_counts = small_mol_counts[end_idx, :]
			final_sm_masses = final_sm_counts * small_mol_mw
			total_monomer_mass = np.sum(final_monomer_masses)
			total_sm_mass = np.sum(final_sm_masses)

			if g > 0:
				prev_counts = monomer_counts[gen_end_index[g - 1], :]
				prev_masses = prev_counts * monomer_mw
				mass_change = final_monomer_masses - prev_masses
				total_change = np.sum(mass_change)
				prev_counts_sm = small_mol_counts[gen_end_index[g - 1], :]
				prev_masses_sm = prev_counts_sm * small_mol_mw
				mass_change_sm = final_sm_masses - prev_masses_sm
				total_change_sm = np.sum(mass_change_sm)
			else:
				mass_change = np.zeros_like(final_monomer_masses)
				total_change = 0.
				mass_change_sm = np.zeros_like(final_sm_masses)
				total_change_sm = 0.

			for i, monomer_id in enumerate(monomer_ids):
				mass = final_monomer_masses[i]
				change = mass_change[i]
				row = {
					"Generation": g,
					"Monomer": monomer_id,
					"Mass (fg)": mass,
					"Mass Change (fg)": change if g > 0 else np.nan,
					"% of Total Protein Mass": (mass / total_monomer_mass * 100) if total_monomer_mass > 0 else np.nan,
					"% of Protein Mass Change": (
								change / total_change * 100) if g > 0 and total_change != 0 else np.nan
				}
				table_rows.append(row)

			for i, small_mol_id in enumerate(small_mol_ids):
				mass = final_sm_masses[i]
				change = mass_change_sm[i]
				row = {
					"Generation": g,
					"Small Molecule": small_mol_id,
					"Mass (fg)": mass,
					"Mass Change (fg)": change if g > 0 else np.nan,
					"% of Total Small Molecule Mass": (mass / total_sm_mass * 100) if total_sm_mass > 0 else np.nan,
					"% of Small Molecule Mass Change": (
								change / total_change_sm * 100) if g > 0 and total_change_sm != 0 else np.nan
				}
				table_rows2.append(row)


		monomer_mass_table = pd.DataFrame(table_rows)
		sm_mass_table = pd.DataFrame(table_rows2)

		monomer_mass_table_sorted = monomer_mass_table.sort_values(
			by="Mass Change (fg)", ascending=False, key=lambda col: np.abs(col.fillna(0))
		).reset_index(drop=True)
		sm_mass_table_sorted = sm_mass_table.sort_values(
			by="Mass Change (fg)", ascending=False, key=lambda col: np.abs(col.fillna(0))
		).reset_index(drop=True)


		sim_seed = self.ap.get_cell_seed(cell_paths[0])
		monomer_mass_table_sorted.to_csv(os.path.join(plotOutDir, f"monomer_mass_by_gen_for_seed_{sim_seed}.csv"),
				index=False)
		sm_mass_table_sorted.to_csv(
			os.path.join(plotOutDir, f"small_molecule_mass_by_gen_for_seed_{sim_seed}.csv"),
			index=False)






if __name__ == '__main__':
	Plot().cli()


