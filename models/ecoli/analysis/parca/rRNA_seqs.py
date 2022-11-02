"""
Plots expression expected with ppGpp versus without.

Useful for troublshooting differences in expression with ppGpp.  Low ribosomal
expression can lead to slower growth than expected (points will appear below
diagonal line).
"""

import os
import pickle

import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objs as go

from models.ecoli.analysis import parcaAnalysisPlot


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		rRNA_ids = sim_data.process.transcription.cistron_data['id'][
            sim_data.process.transcription.cistron_data['is_5S_rRNA']]
		seqs = sim_data.getter.get_sequences(rRNA_ids)
        
		#molecule_groups = sim_data.molecule_groups
		#transcription = sim_data.process
		#seqs = sim_data.getter.get_sequences
        #molecule_groups = sim_data.molecule_groups
        #rRNA_ids = sim_data.process.transcription.cistron_data['id'][sim_data.process.transcription.cistron_data['is_5S_rRNA']]



if __name__ == "__main__":
	Plot().cli()
