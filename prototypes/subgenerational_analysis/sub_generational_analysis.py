"""
This analysis script is intended as a scratch pad to pull out sub-generational gene expression info 
for a set of specified genes, as well as for further data exploration.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/21/2018
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import sys
from datetime import date
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
#from wholecell.containers.bulk_objects_container import BulkObjectsContainer
#from wholecell.analysis.analysis_tools import exportFigure
#from models.ecoli.analysis import multigenAnalysisPlot
#from wholecell.utils.sparkline import whitePadSparklineAxis



def get_todays_date():
	"""
	Want to grab today's date so that data in output folders is not 
	overwritten (just in case its needed later)
	"""
	today = date.today()
	date_str = str(today.month) + str(today.day) + str(today.year)[2:]
	return date_str

def load_simulation_data(sim_data_file):
	f = open(sim_data_file, 'rb')
	sim_data = cPickle.load(f)
	f.close()
	return sim_data

def make_output_dirs(sim_data_path):
	"""
	Store outputs in E. coli output directory, near where simulation data 
	is being kept.

	Make the output directory if it does not already exist
	"""
	output_file_name = get_todays_date()
	output_dir = os.path.join(os.path.dirname(sim_data_path), 
		'subgenerational_analysis', 'outputs', output_file_name)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	return

def get_input_paths(argument): 
	# Check for input file.
	if len(argument) < 2:
		raise Exception, "Oops you didn't provide a path to the data!"
	else:
		sim_data_path = argument[-1]
		sim_data_path_gen = os.path.join(sim_data_path, 'wildtype_000000', '000000')
		#check that parent path is in fact a directory:
		if not os.path.isdir(sim_data_path):
			raise Exception, "input_data_parent_dir does not currently " \
				"exist as a directory"
		sim_data_file = os.path.join(sim_data_path, 'wildtype_000000', 'kb', 
			'simData_Modified.cPickle')
		if not os.path.isfile(sim_data_file):
			raise Exception, "There is something wrong with the data path " \
				"provided. Make sure you are pointing to the main sim output " \
				"directory."
	return sim_data_path, sim_data_path_gen, sim_data_file 

simulation_data_path, simulation_data_path_gen, simulation_data_file = get_input_paths(sys.argv)
sim_data = load_simulation_data(simulation_data_file)
make_output_dirs(simulation_data_path)

#get all cells
ap = AnalysisPaths(simulation_data_path_gen, multi_gen_plot = True)
allDir = ap.get_cells()

rna_ids = sim_data.process.transcription.rnaData["id"]
gene_ids = sim_data.process.transcription.rnaData["geneId"]
protein_ids = sim_data.process.translation.monomerData["id"]

# For each generation
nonzeroSumRnaCounts_allGens = []
nonzeroSumProteinCounts_allGens = []
time = []
time_eachGen = []
for i, simDir in enumerate(allDir):
	simOutDir = os.path.join(simDir, "simOut")
	print 'Got to this point at least'

	time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()
	time_eachGen.append(TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()[0])

	# Read counts of transcripts
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	if i == 0:
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		rnaIndices = np.array([moleculeIds.index(x) for x in rna_ids])
		proteinIndicesAll = np.array([moleculeIds.index(x) for x in proteinIdsall])
	rna_counts = bulkMolecules.readColumn("counts")[:, rnaIndices]
	proteinCountsAll = bulkMolecules.readColumn("counts")[:, proteinIndicesAll]
	bulkMolecules.close()

	# Sum counts over timesteps
	sumRnaCounts = rna_counts.sum(axis = 0)
	sumProteinCounts = proteinCountsAll.sum(axis = 0)

	# Flag where the sum is nonzero (True if nonzero, False if zero)
	nonzeroSumRnaCounts = sumRnaCounts != 0
	nonzeroSumRnaCounts_allGens.append(nonzeroSumRnaCounts)
	nonzeroSumProteinCounts = sumProteinCounts != 0
	nonzeroSumProteinCounts_allGens.append(nonzeroSumProteinCounts)

import pdb; pdb.set_trace()

rna_ids_counts = np.vstack((rna_ids, rna_counts))
protein_ids_counts = np.vstack((proteinIdsall, proteinCountsAll))
with open(os.path.join(plotOutDir, "multi_gen_rnacounts_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, rna_ids_counts, '%s','\t')
with open(os.path.join(plotOutDir, "multi_gen_proteincounts_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, protein_ids_counts, '%s','\t')

# Average (mean) over generations
nonzeroSumRnaCounts_allGens = np.array(nonzeroSumRnaCounts_allGens)
avgRnaCounts = nonzeroSumRnaCounts_allGens.mean(axis = 0)

# Average (mean) monomers over generations
nonzeroSumProteinCounts_allGens = np.array(nonzeroSumProteinCounts_allGens)
avgProteinCounts = nonzeroSumProteinCounts_allGens.mean(axis = 0)

# Identify subgenerationally transcribed genes
subgenRnaIndices = np.where(np.logical_and(avgRnaCounts != 0., avgRnaCounts != 1.))[0]
subgenRnaIds = rnaIds[subgenRnaIndices]
subgenMonomerIndices = [np.where(sim_data.process.translation.monomerData["rnaId"] == x)[0][0] for x in subgenRnaIds]
subgenMonomerIds = sim_data.process.translation.monomerData["id"][subgenMonomerIndices]

# Identify subgenerationally translated genes (not directly tied to subgenerational transcription)
subgenMonomerIndicesAll = np.where(np.logical_and(avgProteinCounts != 0., avgProteinCounts != 1.))[0]
subgenMonomerIdsAll = proteinIdsall[subgenMonomerIndicesAll]

with open(os.path.join(plotOutDir, "subgenerational_rna_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, subgenRnaIds, '%s','\t')
with open(os.path.join(plotOutDir, "subgenerational_protein_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, subgenMonomerIdsAll, '%s','\t')


#load simulation data

