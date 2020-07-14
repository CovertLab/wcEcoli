"""
Plots mRNA expression with operons integration (improved version)

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/7/20
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os
import io

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.io import tsv
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from functools import partial
from reconstruction.spreadsheets import JsonReader, JsonWriter

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        sim_data = cPickle.load(open(simDataFile, "rb"))

        # Listeners used
        main_reader = TableReader(os.path.join(simOutDir, 'Main'))
        bulkMolecules = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
        mRNACounts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))

        # Load data
        mRNA_ids = mRNACounts_reader.readAttribute('mRNA_ids')
        mRNA_counts = mRNACounts_reader.readColumn('mRNA_counts')
        bulkMolecule_ids = bulkMolecules.readAttribute('objectNames')
        bulkMolecule_counts = bulkMolecules.readColumn('counts')
        initial_time = main_reader.readAttribute('initialTime')
        time = main_reader.readColumn('time') - initial_time
        time_total = time.shape[0]

        # Import polycistronic mRNAs from flat file
        FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
        POLY_CISTRON_FILE = os.path.join(FLAT_DIR, 'polycistronic_mrnas_in_model.tsv')
        PC_INFO, pc_fieldnames = self.parse_tsv(POLY_CISTRON_FILE)
        print(PC_INFO)
        print(pc_fieldnames)

    def parse_tsv(self, tsv_file):
        tsv_list = []
        with open(tsv_file) as tsvfile:
            reader = JsonReader(tsvfile)
            fieldnames = reader.fieldnames
            for row in reader:
                tsv_list.append(row)
        return tsv_list, fieldnames


if __name__ == '__main__':
    Plot().cli()
