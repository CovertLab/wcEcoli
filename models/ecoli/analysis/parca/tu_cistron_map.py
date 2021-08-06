"""
This script generates a plot that shows the transcription units and their respective cistrons,
with the cistrons in the same transcription units adjacent/close to each other, and with TUs with
similar cistrons close to each other. This helps to portray "operons".
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):

	# depth first search
	def dfs(self, visited, graph, node):
		if node not in visited:
			visited.append(node)
			for neighbor in graph[node]:
				self.dfs(visited, graph, neighbor)
		return visited


	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# can omit the following block if tu_id_to_cistron_ids is an attribute in process.transcription
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = pickle.load(f)
		all_valid_tus = [tu for tu in raw_data.transcription_units if sim_data.getter.is_valid_molecule(tu['id'])]
		rna_ids = [cistron_id for cistron_id in sim_data.process.transcription.cistron_data['id'] if
				   sim_data.getter.is_valid_molecule(cistron_id)]
		rna_ids.extend([tu['id'] for tu in all_valid_tus])
		gene_id_to_rna_id = {gene['id']: gene['rna_id'] for gene in raw_data.genes}
		tu_id_to_cistron_ids = {tu['id']: [gene_id_to_rna_id[gene] for gene in tu['genes']]for tu in all_valid_tus}
		# block ends, tu_id_to_cistron_ids only contains one TU for this test case

		rna_id_to_gene_id = {
			gene['rna_id']: gene['id'] for gene in raw_data.genes}
		gene_id_to_left_end_pos = {
			gene['id']: gene['left_end_pos'] for gene in raw_data.genes
		}
		gene_id_to_right_end_pos = {
			gene['id']: gene['right_end_pos'] for gene in raw_data.genes
		}

		new_dic = {}
		keylist = tu_id_to_cistron_ids.keys()
		for i in keylist:
			new_dic[i] = []
		overall_extra = 0

		# rearrange new dictionary to where tu1:[cistron1,cistron2]
		for TU_units in tu_id_to_cistron_ids.keys():
			cistron_count = []
			for cistron in tu_id_to_cistron_ids[TU_units]:
				cistron_count.append(len(sim_data.process.transcription.cistron_id_to_rna_indexes(cistron))) #get how many TU this cistron is it
			sort = np.array(cistron_count)
			sort_index = np.argsort(sort)
			overall_extra += 1
			for index in sort_index:
				new_dic[TU_units].append(np.array(tu_id_to_cistron_ids[TU_units])[index])

		# keeping track of the TU and cistron (x and y axis) ID
		matrix_TU_id = []
		matrix_cis_id = []

		# initialize matrix for output
		columns = len(new_dic.keys())
		rows = 0
		for tu in new_dic.keys():
			for cis in new_dic[tu]:
				rows+=1
		new_matrix = [[0 for i in range(columns)] for j in range(rows)]

		# initialize i and j count
		x = 0
		y = 0

		# generate graph for dfs, the nodes are TU, graph relates one TU to another TU that shares cistron
		graph = {}
		for i in keylist:
			graph[i] = []
		for tu in graph.keys():
			tu_related = []
			cist_list = new_dic[tu]
			for cis in cist_list:
				tu_related_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes(cis)
				for indexes in list(set(tu_related_indexes)):
					if rna_ids[indexes] != tu:
						if not rna_ids[indexes] in tu_related:
							tu_related.append(rna_ids[indexes])
			for each in list(set(tu_related)):
				graph[tu].append(each)



		track = list(graph.keys())
		groups = []
		while len(track) > 0:
			unit = track[0]
			visited = []
			vis = self.dfs(visited,graph,unit) # list
			# get rid of duplicates in list
			temp_list = []
			for i in vis:
				if i not in temp_list:
					temp_list.append(i)
			vis = temp_list
			for group in vis:
				if group in track:
					track.remove(group)
			groups.append(vis)

		for operon in groups:
			for tu in operon:
				matrix_TU_id.append(tu)
				cistrons = new_dic[tu]
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
				y += 1

		# convert from matrix to scatter, block out if want matrix
		b = []
		a = []
		for r in range(0,rows,1):
			for c in range(0,columns,1):
				if new_matrix[r][c] == 1:
					b.append(r)
					a.append(c)

		plt.figure(figsize=(70,150))
		plt.scatter(a, b, c='blue',marker = 's') # block out for matrix in line 159
		plt.title("Transcription units and the corresponding cistrons", fontsize=20)
		plt.xlabel('TU', fontsize=15)
		plt.ylabel('cistron', fontsize=15)
		# plt.imshow(new_matrix) # black and white binary

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
