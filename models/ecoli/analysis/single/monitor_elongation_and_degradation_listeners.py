"""
Temporary script for comparing the counts of one elongation listener method to another.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
from ecoli.library.parquet_emitter import read_stacked_bulk_molecules, read_stacked_columns
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}
		# Extract the protein IDs
		monomerIDs = monomer_counts_reader.readAttribute("monomerIds")

		# Extract the total monomer counts:
		total_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomerCounts")

		# Extract the free monomer counts using the monomer counts listener:
		free_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")

		# Extract how many proteins were removed via degradation over the entire sim length:
		degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomersDegraded")

		# Extract how many counts were added via elongation over the entire sim length:
		elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts',
												"monomersElongated")  # can I just use the monomer counts reader like above?
		# Extract terminated counts:
		terminated_counts = read_stacked_columns(cell_paths, 'MonomerCounts',
												"monomersTerminated")

		#NOTE: I do not expect that for every time step, the free monomer counts will equal exactly
		# the degradaed + the elongated counts becuase complexation does happen and subtracts to the
		# total amount of free monomers available. however, if the listeners are working as expected,
		# then terminated and elongated counts should be equal.

		# Extract the sim time:
		time = read_stacked_bulk_molecules(cell_paths, 'Time', 'time')

		# first comparison:
		def simple_compare(elongated_counts, terminated_counts):
			diff = elongated_counts - terminated_counts
			return diff

		def compare(diff):
			if np.allclose(diff, 0):
				print("Elongated and Terminated counts match exactly at all time points!")
			else:
				print("Elongated and Terminated counts do NOT match exactly at all time points.")
				max_diff = np.max(np.abs(diff))
				print(f"Maximum difference between the two methods at any time point: {max_diff} counts.")

		diff = simple_compare(elongated_counts, terminated_counts)
		compare(diff)

		# Compare the two methods more with the degraded counts too:
		def compare_with_degraded_counts(elongated_counts, terminated_counts, degraded_counts):
			total_via_elongation = elongated_counts - degraded_counts
			total_via_termination = terminated_counts - degraded_counts
			return total_via_elongation, total_via_termination

		total_via_elongation, total_via_termination = compare_with_degraded_counts(
			elongated_counts, terminated_counts, degraded_counts)

		def final_compare(total_via_elongation, total_via_termination):
			if np.allclose(total_via_elongation, total_via_termination):
				print("After accounting for degradation, both methods yield the same total monomer counts at all time points!")
			else:
				print("After accounting for degradation, the two methods do NOT yield the same total monomer counts at all time points.")
				max_diff = np.max(np.abs(total_via_elongation - total_via_termination))
				print(f"Maximum difference between the two methods at any time point after accounting for degradation: {max_diff} counts.")

		final_compare(total_via_elongation, total_via_termination)


		# Finally, just to see what this looks like, see how the free monomer counts change releative to elongation and degradation:
		# this is for my own curiosity, not part of the comparison:
		compare_terminated_with_FMCs = []
		for i in range(len(time)):
			if i == 0:
				# at time 0, just record the free monomer counts
				compare_terminated_with_FMCs.append(free_monomer_counts[i])
			else:
				# for all other time pointts, monitor how fmcs change relative to the value at the previous timestep:
				delta_fmc = free_monomer_counts[i-1] + terminated_counts[i] - degraded_counts[i]
				compare_terminated_with_FMCs.append(delta_fmc)

		# Compare with elongated:
		compare_elongated_with_FMCs = []
		for i in range(len(time)):
			if i == 0:
				# at time 0, just record the free monomer counts
				compare_terminated_with_FMCs.append(free_monomer_counts[i])
			else:
				# for all other time pointts, monitor how fmcs change relative to the value at the previous timestep:
				delta_fmc = free_monomer_counts[i - 1] + elongated_counts[i] - degraded_counts[i]
				compare_elongated_with_FMCs.append(delta_fmc)

		# compare the two:
		compare_both_with_FMCs = np.array(compare_terminated_with_FMCs) - np.array(compare_elongated_with_FMCs)
		if np.allclose(compare_both_with_FMCs, 0):
			print("Both methods yield the same free monomer counts over the simulation when accounting for degradation.")
		else:
			print("Both methods do NOT yield the same free monomer counts over the simulation when accounting for degradation.")
			max_diff = np.max(np.abs(compare_both_with_FMCs))
			print(f"Maximum difference between the two methods' free monomer counts at any time point after accounting for degradation: {max_diff} counts.")




		# also see if you can compare how big the impact of complexes is:
		compare_fmcs_to_previous_timepoint = []
		for i in range(len(time)):
			if i == 0:
				compare_fmcs_to_previous_timepoint.append(np.zeros(len(monomerIDs)))
			else:
				delta_fmc = free_monomer_counts[i] - free_monomer_counts[i-1]
				compare_fmcs_to_previous_timepoint.append(delta_fmc)

		# understand how complexation has a role in changing free monomer counts:
		terminated_minus_degraded = terminated_counts - degraded_counts
		complexation_impact = compare_fmcs_to_previous_timepoint - terminated_minus_degraded
		if np.allclose(complexation_impact, 0):
			print("Complexation has no impact on free monomer counts over the simulation.")
		else:
			print("Complexation impacts free monomer counts over the simulation.")
			max_impact = np.max(np.abs(complexation_impact))
			print(f"Maximum impact of complexation on free monomer counts at any time point:{max_impact} counts.")


	# finally, compare total monomer counts to free monomer counts:
		def compare_total_to_free(total_monomer_counts, free_monomer_counts):
			diff = total_monomer_counts - free_monomer_counts
			return diff
		diff_total_vs_free = compare_total_to_free(total_monomer_counts, free_monomer_counts)
		if np.allclose(diff_total_vs_free, 0):
			print("Total monomer counts equal free monomer counts at all time points!")
		else:
			print("Total monomer counts do NOT equal free monomer counts at all time points.")
			max_diff = np.max(np.abs(diff_total_vs_free))
			print(f"Maximum difference between total and free monomer counts at any time point: {max_diff} counts.")




if __name__ == '__main__':
	Plot().cli()
