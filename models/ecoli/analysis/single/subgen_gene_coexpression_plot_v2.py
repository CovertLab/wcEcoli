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
from wholecell.io import tsv
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
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
        mRNA_gene_dict = {}
        exp_time_df = pd.DataFrame(columns=['tu_unit', 'gene_id', 'time_expressed', 'percent_time_expressed']) # df to store % time expression for each gene

        for pc in PC_INFO:
            pc_gene_id = '_'.join(pc['transcription_units']) # get TU id
            pc_mRNA_ids = [pc_gene_id + "_RNA"] # TU id -> polycistronic mRNA id

            # not all monomers may be removed. find the ones that are expressed on its own
            non_removed_monomers = list(set(pc['transcription_units']) ^ set(pc['monomers_to_remove']))
            non_removed_mRNA_ids = [x + "_RNA" for x in non_removed_monomers]
            pc_mRNA_ids.extend(non_removed_mRNA_ids) # names of all transcribed single mRNAs

            mRNA_protein_dict[pc_gene_id] = {}
            mRNA_gene_dict[pc_gene_id] = pc['transcription_units']
            polycistron_exp_count = 0
            polycistron_percent_time = 0

            for rna in RNA_INFO:
                if rna['id'] in pc_mRNA_ids:
                    location = rna['location'][0]
                    mRNA_with_location = rna['id'] + "[" + location + "]"
                    pc_mRNA_ids = [mRNA_with_location if x==rna['id'] else x for x in pc_mRNA_ids] # add mRNA location

                    # add gene expression time to dataframe, to plot later
                    mRNA_exp_count = np.count_nonzero(mRNA_counts[:, mRNA_ids.index(mRNA_with_location)])
                    mRNA_percent_time = (mRNA_exp_count / time_total)
                    exp_time_df = exp_time_df.append({'tu_unit': pc_gene_id, 'gene_id': rna['id'][:-4], 'time_expressed': mRNA_exp_count, 'percent_time_expressed': mRNA_percent_time}, ignore_index=True)

                    # store gene expression time specifically for polycistronic gene
                    if rna['id'][:-4] == pc_gene_id: # check that current mRNA is not a non-removed monomer
                        polycistron_exp_count = mRNA_exp_count
                        polycistron_percent_time = mRNA_percent_time

                    mRNA_protein_dict[pc_gene_id][mRNA_with_location] = []
                    for protein in PROTEIN_INFO:
                        if protein['id'] in rna['monomerSet']: # check if protein is expressed by mRNA
                            protein_with_location = protein['id'] + "[" + protein['location'][0] + "]"
                            mRNA_protein_dict[pc_gene_id][mRNA_with_location].append(protein_with_location) # add protein monomer with location

            # add to dataframe: single genes for monomers that were removed
            for gene in pc['transcription_units']:
                if gene not in exp_time_df.values: # gene not already in dataframe
                    exp_time_df = exp_time_df.append({'tu_unit': pc_gene_id, 'gene_id': gene, 'time_expressed': polycistron_exp_count, 'percent_time_expressed': polycistron_percent_time}, ignore_index=True)
                else: # gene in dataframe (i.e. it was a removed monomer and added earlier, see if expression data should be replaced with higher number)
                    gene_row = exp_time_df.loc[exp_time_df['gene_id'] == gene]
                    if polycistron_percent_time > gene_row['percent_time_expressed'].values[0]:
                        exp_time_df.at[gene_row.index.item(), 'time_expressed'] = polycistron_exp_count
                        exp_time_df.at[gene_row.index.item(), 'percent_time_expressed'] = polycistron_percent_time


        # Make multi-page plot, one page per operon
        with PdfPages(os.path.join(plotOutDir, 'mRNA_protein_expression.pdf')) as pdf:
            for plotIndex, tu in enumerate(mRNA_protein_dict.keys()):
                fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12, 8))
                fig.add_subplot(111, frameon=False)
                plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

                # setting up table
                proteins_plotted = []

                for mRNA in mRNA_protein_dict[tu]:
                    mRNA_count = mRNA_counts[:, mRNA_ids.index(mRNA)] # retrieve array of counts over time
                    mRNA_percent_time = (np.count_nonzero(mRNA_count) / time_total) * 100 # get % of time expressed
                    ax1.plot(mRNA_counts[:, mRNA_ids.index(mRNA)], label='\n'.join(wrap(mRNA, 20)))

                    # plotting proteins translated from mRNA
                    for protein in mRNA_protein_dict[tu][mRNA]:
                        if protein not in proteins_plotted:
                            proteins_plotted.append(protein) # add row to table
                            protein_count = bulkMolecule_counts[:,bulkMolecule_ids.index(protein)]
                            protein_percent_time = (np.count_nonzero(protein_count) / time_total) * 100 # get % of time expressed
                            ax2.plot(protein_count, label='\n'.join(wrap(protein + "(% t_total: " + str(protein_percent_time) + "%)", 20)))
                            # now account for complexes
                            complex_list = self.get_protein_complexes(protein[:-3], COMPLEX_INFO) # truncating protein location
                            if complex_list:
                                for cplx in complex_list:
                                    cplx_count = bulkMolecule_counts[:,bulkMolecule_ids.index(cplx)]
                                    cplx_percent_time = (np.count_nonzero(cplx_count) / time_total) * 100 # get % of time expressed
                                    ax2.plot(cplx_count, label='\n'.join(wrap(cplx + "(% t_total: " + str(cplx_percent_time) + "%)", 20)))

                text_y_pos = 0.95
                for gene in mRNA_gene_dict[tu]:
                    gene_row = exp_time_df[exp_time_df['gene_id'] == gene]
                    percent_time_exp = '%.3f'%(gene_row['percent_time_expressed'].values[0] * 100)
                    ax1.text(0.05, text_y_pos, gene + " = " + str(percent_time_exp) + "%", transform=plt.gca().transAxes)
                    text_y_pos = text_y_pos - 0.03

                ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left') # mRNA plot
                ax2.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left') # protein plot

                ax1.set_title(tu + " mRNAs")
                ax2.set_title(tu + " proteins")
                plt.xlabel("Time (sec)")
                plt.ylabel("count")
                plt.tight_layout()
                pdf.savefig()
                plt.close('all')

        self.write_time_tsv(plotOutDir, mRNA_counts, mRNA_ids, mRNA_protein_dict)

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

    def write_time_tsv(self, plotOutDir, mRNA_counts, mRNA_ids, mRNA_protein_dict):
        tsvFile = io.open(os.path.join(plotOutDir, "mRNA_percent_time_expressed.tsv"), "wb")
        output = tsv.writer(tsvFile)
        output.writerow(['tu_id', 'mRNA_id', 'protein_id', 'time_coexp', 'percent_time_coexp'])
        # convert dict to tuple, easier to write to tsv
        info_to_write = [(tu, mRNA, protein) for tu in mRNA_protein_dict for mRNA in mRNA_protein_dict[tu] for protein in mRNA_protein_dict[tu][mRNA]]
        for tu, mRNA, protein in info_to_write:
            mRNA_count = mRNA_counts[:, mRNA_ids.index(mRNA)] # retrieve array of counts over time (to do: make this less redundant with do_plot function)
            mRNA_coexp_time = np.count_nonzero(mRNA_count)
            mRNA_percent_time = (mRNA_coexp_time / len(mRNA_count)) * 100 # get % of time expressed
            output.writerow([tu, mRNA, protein, mRNA_coexp_time, mRNA_percent_time])
        return

if __name__ == '__main__':
    Plot().cli()
