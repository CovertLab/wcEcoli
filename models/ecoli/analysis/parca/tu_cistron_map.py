"""
Template for parca analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		#can omit the following block if set tu_id_to_cistron_ids as an attribute in process.transcription
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = pickle.load(f)
		all_valid_tus = [tu for tu in raw_data.transcription_units if sim_data.getter.is_valid_molecule(tu['id'])]
		rna_ids = [cistron_id for cistron_id in sim_data.process.transcription.cistron_data['id'] if
				   sim_data.getter.is_valid_molecule(cistron_id)]
		rna_ids.extend([tu['id'] for tu in all_valid_tus])
		n_rnas = len(rna_ids)
		gene_id_to_rna_id = {gene['id']: gene['rna_id'] for gene in raw_data.genes}
		tu_id_to_cistron_ids = {tu['id']: [gene_id_to_rna_id[gene] for gene in tu['genes']]for tu in all_valid_tus}
		#block ends, tu_id_to_cistron_ids only contains one TU for this test case

		#matrix = sim_data.process.transcription.cistron_tu_mapping_matrix()  # this should have the matrix organized by operons blocks
		cistron_id = sim_data.process.transcription.cistron_data['id']

		new_dic = {}
		#initialize dictionary with TU keys
		keylist = tu_id_to_cistron_ids.keys()
		for i in keylist:
			new_dic[i] = []
		overall_extra = 0
		#rearrange new dictionary to desired order of cistrons for each TU
		for TU_units in tu_id_to_cistron_ids.keys():
			cistron_count = []
			for cistron in tu_id_to_cistron_ids[TU_units]:
				cistron_count.append(sim_data.process.transcription.cistron_id_to_tu_indexes(cistron)) #get how many TU this cistron is it
			sort = np.array(cistron_count)
			sort_index = np.argsort(sort)
			overall_extra += 1
			for index in sort_index:
				new_dic[TU_units].append(np.array(tu_id_to_cistron_ids[TU_units])[index])

		#keeping track of the TU and cistron (x and y axis) ID
		matrix_TU_id = []
		matrix_cis_id = []

		#initialize matrix
		columns = 6000 #len(rna_ids)+overall_extra
		rows = 6000 #len(rna_ids)+overall_extra
		new_matrix = [[0 for i in range(columns)] for j in range(rows)]

		#initialize i and j count
		x = 0
		y = 0

		#fill matrix by looping through TU and ordering the cist so cistrons in the same TU are together
		while len(new_dic.keys()) > 0:
			next_TU = [list(new_dic.keys())[0]]
			for TU in next_TU:
				matrix_TU_id.append(TU)
				cistrons = new_dic[TU]
				for cist in cistrons:
					hold = x
					if cist in matrix_cis_id:
						x = matrix_cis_id.index(cist)
						#get the index and set x to that
						new_matrix[x][y] = 1
						#reset x to last x before this
						x = hold
					else:
						matrix_cis_id.append(cist)
						new_matrix[x][y] = 1
						x+=1
					overlap = sim_data.process.transcription.cistron_id_to_tu_indexes(cist)
					if len(overlap) != 0:
						next_TU = overlap
				del new_dic[TU]
				y += 1

		#delete the rna ids that are not mono
		#num = 0
		#while num < len(matrix_cis_id):
			#rna_ids.remove(matrix_cis_id[num][0])
			#num += 1
	    # already not in the list????

		#fill the rest (mono)
		#grab the rest of the TU and their respective cistrons and fill in the rest of the matrix
		for rna_id in rna_ids:
			matrix_TU_id.append(rna_id) #need to fix this
			matrix_cis_id.append(rna_id)
			new_matrix[x][y] = 1
			x += 1
			y += 1

		# desired x and y axis: TU (mRNA) vs gene(individual/cistron)
		# need to organize so that the diagonal shows the operon groups
		# Monocistrons at the end of each operon group
		# ipdb.set_time()
		# bigger plot!
		#for each gene, look at how many TU it is in loop through TU, find the genes in the TU, for each gene, look at how many TU it is in, order the ones not in other TUs first, fill in the new matrix with the genes in TU in order, then continue to the next TU (prioritize the ones that contains the genes in the prior round),
		#sort by rows to make sure the genes that are in the same operon are together
		plt.figure()
		plt.title("Transcription units and the corresponding cistrons")
		plt.xlabel('TU')
		plt.ylabel('cistron')
		plt.imshow(new_matrix) # black and white binary
		plt.show()
		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
