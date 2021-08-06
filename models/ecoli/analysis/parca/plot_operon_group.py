"""
grid of all operon
pdf search for an operon
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, input_string, type, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = pickle.load(f)
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		all_valid_tus = [tu for tu in raw_data.transcription_units if sim_data.getter.is_valid_molecule(tu['id'])]
		rna_ids = [cistron_id for cistron_id in sim_data.process.transcription.cistron_data['id'] if
				   sim_data.getter.is_valid_molecule(cistron_id)]
		rna_ids.extend([tu['id'] for tu in all_valid_tus])
		n_rnas = len(rna_ids)
		gene_id_to_rna_id = {gene['id']: gene['rna_id'] for gene in raw_data.genes}
		tu_id_to_cistron_ids = {tu['id']: [gene_id_to_rna_id[gene] for gene in tu['genes']] for tu in all_valid_tus}
		# block ends, tu_id_to_cistron_ids only contains one TU for this test case

		cistron_id = sim_data.process.transcription.cistron_data['id']
		rna_id_to_gene_id = {
			gene['rna_id']: gene['id'] for gene in raw_data.genes}
		gene_id_to_left_end_pos = {
			gene['id']: gene['left_end_pos'] for gene in raw_data.genes
		}
		gene_id_to_right_end_pos = {
			gene['id']: gene['right_end_pos'] for gene in raw_data.genes
		}
		all_cistrons = [
			rna for rna in raw_data.rnas
			if rna['id'] in rna_id_to_gene_id
			   and gene_id_to_left_end_pos[rna_id_to_gene_id[rna['id']]] is not None
			   and gene_id_to_right_end_pos[rna_id_to_gene_id[rna['id']]] is not None
		]
		n_cistrons = len(all_cistrons)

		new_dic = {}
		# initialize dictionary with TU keys
		keylist = tu_id_to_cistron_ids.keys()
		for i in keylist:
			new_dic[i] = []
		overall_extra = 0

		# rearrange new dictionary to desired order of cistrons for each TU
		for TU_units in tu_id_to_cistron_ids.keys():
			cistron_count = []
			for cistron in tu_id_to_cistron_ids[TU_units]:
				cistron_count.append(len(sim_data.process.transcription.cistron_id_to_rna_indexes(
					cistron)))  # get how many TU this cistron is it
			sort = np.array(cistron_count)
			sort_index = np.argsort(sort)
			overall_extra += 1
			for index in sort_index:
				new_dic[TU_units].append(np.array(tu_id_to_cistron_ids[TU_units])[index])

		# keeping track of the TU and cistron (x and y axis) ID
		matrix_TU_id = []
		matrix_cis_id = []

		# initialize matrix
		columns = len(new_dic.keys())  # +len(rna_ids)
		rows = n_cistrons  # +3000 #NEED TO FIX
		new_matrix = [[0 for i in range(columns)] for j in range(rows)]

		# initialize i and j count
		x = 0
		y = 0

		# fill matrix by looping through TU and ordering the cist so cistrons in the same TU are together
		while len(new_dic.keys()) > 0:
			next_TU = [list(new_dic.keys())[0]]
			for TU in next_TU:
				matrix_TU_id.append(TU)
				cistrons = new_dic[TU]
				for cist in cistrons:
					hold = x
					if cist in matrix_cis_id:
						x = matrix_cis_id.index(cist)
						new_matrix[x][y] = 1
						x = hold
					else:
						matrix_cis_id.append(cist)
						new_matrix[x][y] = 1
						x += 1
					overlap = sim_data.process.transcription.cistron_id_to_rna_indexes(cist)
					if len(overlap) != 0:
						next_TU = overlap
				del new_dic[TU]
				y += 1

		#get position of input TU or input cistron in the matrix
		position = -1
		if type == "TU":
			try:
				position = matrix_TU_id.index(input_string)
			except ValueError:
				print("Input transcription unit does not exist in an operon group in the data")
		if type == "cistron":
			try:
				position = matrix_cis_id.index(input_string)
			except ValueError:
				print("Input cistron does not exist in an operon group in the data")

		#With the position of the TU or cistron, continuously search neighbor blocks until boundary of operon group
		range = [0,0,0,0]
		if type == 'TU':
			neigh = True
			y_start = matrix_cis_id.index(new_dic[input_string][0])
			x_start = position
			start = [x_start, y_start]
			top = start
			bottom = start
			left = start
			right = start
			while neigh == True:
				contin = False
				if new_matrix[top[0]][top[1]-1] == 1:
					top = [top[0],top[1]-1]
					contin = True
				if new_matrix[top[0]][top[1]+1] == 1:
					bottom = [bottom[0],bottom[1]+1]
					contin = True
				if new_matrix[top[0]-1][top[1]] == 1:
					left = [left[0]-1, left[1]]
					contin = True
				if new_matrix[top[0]+1][top[1]] == 1:
					right = [right[0]+1,right[1]]
					contin = True
				if contin == False:
					neigh = False #breaks out of search
			corners = [top,bottom,left,right]
			for pos in corners:
				if pos[0] < range[0]:
					range[0] = pos[0]
				if pos[0] > range[1]:
					range[1] = pos[0]
				if pos[1] < range[2]:
					range[2] = pos[1]
				if pos[1] > range[3]:
					range[3] = pos[1]

		if type == 'cistron':


		plt.figure()
		sub_matrix = new_matrix[range[0]:range[1],range[2]:range[3]]
		im.show(sub_matrix)
		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
