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
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import itertools
from textwrap import wrap

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.spreadsheets import JsonReader, JsonWriter

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        # sim_data = cPickle.load(open(simDataFile, "rb"))

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

        # Import TUs, mRNAs, proteins, and protein complexes from flat files
        FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
        POLY_CISTRON_FILE = os.path.join(FLAT_DIR, 'polycistronic_mrnas_in_model.tsv')
        RNA_FILE = os.path.join(FLAT_DIR, 'operon_rnas.tsv')
        PROTEIN_FILE = os.path.join(FLAT_DIR, 'proteins.tsv')
        COMPLEXATION_FILE = os.path.join(FLAT_DIR, 'complexationReactions.tsv')
        PC_INFO, pc_fieldnames = self.parse_tsv(POLY_CISTRON_FILE)
        RNA_INFO, rna_fieldnames = self.parse_tsv(RNA_FILE)
        PROTEIN_INFO, protein_fieldnames = self.parse_tsv(PROTEIN_FILE)
        COMPLEX_INFO, complex_fieldnames = self.parse_tsv(COMPLEXATION_FILE)

        # Get all mRNAs transcribed from TUs + protein monomers
        mRNA_protein_dict = {}
        for pc in PC_INFO:
            pc_gene_id = '_'.join(pc['transcription_units'])
            pc_mRNA_ids = [pc_gene_id + "_RNA"]
            non_removed_monomers = list(set(pc['transcription_units']) ^ set(pc['monomers_to_remove']))
            non_removed_mRNA_ids = [x + "_RNA" for x in non_removed_monomers]
            pc_mRNA_ids.extend(non_removed_mRNA_ids) # names of all transcribed mRNAs
            mRNA_protein_dict[pc_gene_id] = {}

            for rna in RNA_INFO:
                if rna['id'] in pc_mRNA_ids:
                    location = rna['location'][0]
                    mRNA_with_location = rna['id'] + "[" + location + "]"
                    pc_mRNA_ids = [mRNA_with_location if x==rna['id'] else x for x in pc_mRNA_ids] # add mRNA location
                    mRNA_protein_dict[pc_gene_id][mRNA_with_location] = []
                    for protein in PROTEIN_INFO:
                        if protein['id'] in rna['monomerSet']:
                            protein_with_location = protein['id'] + "[" + protein['location'][0] + "]"
                            mRNA_protein_dict[pc_gene_id][mRNA_with_location].append(protein_with_location) # add protein monomer with location

        # print(mRNA_protein_dict)

        # Make plot, one for each operon
        with PdfPages(os.path.join(plotOutDir, 'mRNA_protein_expression.pdf')) as pdf:
            for plotIndex, tu in enumerate(mRNA_protein_dict.keys()):
                fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12, 8))
                fig.add_subplot(111, frameon=False)
                plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

                for mRNA in mRNA_protein_dict[tu]:
                    ax1.plot(mRNA_counts[:, mRNA_ids.index(mRNA)], label='\n'.join(wrap(mRNA, 20)))
                    for protein in mRNA_protein_dict[tu][mRNA]:
                        ax2.plot(bulkMolecule_counts[:,bulkMolecule_ids.index(protein)], label='\n'.join(wrap(protein, 20)))
                        # now account for complexes
                        complex_list = self.get_protein_complexes(protein[:-3], COMPLEX_INFO) # truncating protein location
                        if complex_list:
                            for cplx in complex_list:
                                ax2.plot(bulkMolecule_counts[:,bulkMolecule_ids.index(cplx)], label='\n'.join(wrap(cplx, 20)))
                ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
                ax2.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
                ax1.set_title(tu + " mRNAs")
                ax2.set_title(tu + " proteins")
                plt.xlabel("Time (sec)")
                plt.ylabel("count")
                plt.tight_layout()
                pdf.savefig()
                plt.close('all')

    def parse_tsv(self, tsv_file):
        tsv_list = []
        with open(tsv_file) as tsvfile, open('temp.csv', 'w') as tempfile:
            tsvfile_start = tsvfile.tell()
            if not tsvfile.readline().startswith('#'):
                tsvfile.seek(tsvfile_start)
                reader = JsonReader(tsvfile)
                fieldnames = reader.fieldnames
                for row in reader:
                    tsv_list.append(row)
            else:
                for this_line in tsvfile:
                    if not this_line.startswith('#'): tempfile.write(this_line)
                tempfile.close()
                with open('temp.csv') as tempfile:
                    reader = JsonReader(tempfile, dialect="excel-tab")
                    fieldnames = reader.fieldnames
                    for row in reader:
                        tsv_list.append(row)
                if os.path.exists("temp.csv"): os.remove("temp.csv")
        return tsv_list, fieldnames

    def get_protein_complexes(self, protein, complex_info):
        complex_list = []
        for complex_row in complex_info:
            for stoich in complex_row['stoichiometry']:
                if stoich['molecule']==protein:
                    complex_name = complex_row['stoichiometry'][0]['molecule'] + "[" + complex_row['stoichiometry'][0]['location'] + "]"
                    complex_list.append(complex_name)
        return complex_list



if __name__ == '__main__':
    Plot().cli()
