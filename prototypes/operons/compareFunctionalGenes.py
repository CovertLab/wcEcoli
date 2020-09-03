import csv
import os
import numpy as np
from collections import defaultdict
from prototypes.operons.get_operon_rna_list import pc_mc_conversion_dicts
import matplotlib.pyplot as plt
import argparse
from wholecell.io.tablereader import TableReader
from functools import partial
from reconstruction import spreadsheets

from wholecell.analysis.analysis_tools import exportFigure

def inspect_functional_genes(sim1Name, sim1Dir, sim2Name, sim2Dir, pc_rnas_path):
    # check if functional_genes.tsv exists
    if not os.path.isfile('runscripts/reflect/functional_genes.tsv'):
        print("functional_genes.tsv does not exist. create this file using runscripts/reflect/model_inspection")
        return


    with open('runscripts/reflect/functional_genes.tsv') as f:

        tsvreader = csv.reader(f, delimiter='\t')

        process_dict = defaultdict(list)

        lcount = 0
        for line in tsvreader:

            if lcount > 1:
                process_dict[line[3]].append(line[1])

            lcount += 1

    # some of the annotations seem redundant, so let's concatenate them
    process_dict['Transcription'] = process_dict['Transcription'] + process_dict['Transcription Regulation']
    process_dict['Metabolism'] = process_dict['Metabolism'] + process_dict['Metabolism, RNA Decay']
    process_dict['Metabolism'] = process_dict['Metabolism'] + process_dict['Transcription Regulation, Metabolism']
    process_dict['RNA Decay'] = process_dict['RNA Decay'] + process_dict['Metabolism, RNA Decay']
    process_dict['Transcription'] = process_dict['Transcription'] + process_dict['Transcription Regulation, Metabolism']

    # LOAD RNA AND PROTEIN DNA
    mRNA_names = {}
    mRNA_counts = {}
    monomer_counts = {}
    monomer_names = {}

    for simName, simDir in zip([sim1Name, sim2Name], [sim1Dir, sim2Dir]):
        mRNA_counts_reader = TableReader(os.path.join(simDir, 'mRNACounts'))
        mRNA_names[simName] = mRNA_counts_reader.readAttribute('mRNA_ids')
        mRNA_counts[simName] = mRNA_counts_reader.readColumn('mRNA_counts')

        monomer_counts_reader = TableReader(os.path.join(simDir, "MonomerCounts"))
        monomer_counts[simName] = monomer_counts_reader.readColumn("monomerCounts")
        monomer_names[simName] = monomer_counts_reader.readAttribute("monomerIds")


    # ------------------------------------------------------------------------------------------------------------
    # PLOTS
    # plot the RNA and the protein levels for all the functional genes

    keys_to_plot = ['Metabolism', 'RNA Decay', 'Transcription', 'Translation']

    # RNA PLOT
    rna_delta = np.empty(shape=(0,6)) # array % change in rna expression vs sim2, write to tsv later
    rd = np.empty(shape=(0,1))
    rna_r2 = {}
    pc_to_mc_dict, mc_to_pc_dict = pc_mc_conversion_dicts(pc_rnas_path)
    plt.figure(figsize=(15,9))
    fig, axs = plt.subplots(2, 4, dpi=300)
    for p, key in enumerate(keys_to_plot):

        plotdata = np.empty(shape=(0,2)) # save data of each plot to get r2

        for rna in process_dict[key]:

            rna_name = rna + '[c]'

            if rna_name in mc_to_pc_dict:
                rnas_s1 = mc_to_pc_dict[rna_name]
            else:
                rnas_s1 = [rna_name]

            rnas_s2 = rna_name

            for r in rnas_s1:
                x = np.mean(mRNA_counts[sim2Name][:, mRNA_names[sim2Name].index(rnas_s2)])
                y = np.mean(mRNA_counts[sim1Name][:, mRNA_names[sim1Name].index(r)])
                axs[0, p].scatter(x, y, c='w', edgecolor='k', alpha=.7, s=7)

                if x > 0:
                    s1_s2_delta = round(y/x,6)
                else:
                    s1_s2_delta = np.nan

                row = [rnas_s2, r,  s1_s2_delta, round(x,6), round(y,6), key]
                rna_delta = np.vstack((rna_delta, row))
                rd = np.vstack((rd, s1_s2_delta))
                plotdata = np.vstack((plotdata, [x,y]))
        # r2 to write to tsv
        corr_mat = np.corrcoef(plotdata[:,0], plotdata[:,1])
        rna_r2[key] = corr_mat[0,1]



        equalize_axes(axs[0, p])
        axs[0, p].set_aspect('equal')
        axs[0, p].set_xlabel(sim2Name + ' RNA', fontsize=6)
        axs[0, p].set_ylabel(sim1Name + ' RNA', fontsize=6)
        axs[0, p].set_title(key, fontsize=8)
        axs[0,p].tick_params(axis='both', labelsize=6)


    # PROTEIN PLOT
    rna_to_protein_dict = build_rna_to_protein_dict()
    prot_r2 = {}
    prot_delta = np.empty(shape=(0,5))
    for p, key in enumerate(keys_to_plot):
        plotdata = np.empty(shape=(0,2))
        for rna in process_dict[key]:

            x = np.mean(monomer_counts[sim2Name][:, monomer_names[sim2Name].index(rna_to_protein_dict[rna])])
            y = np.mean(monomer_counts[sim1Name][:, monomer_names[sim1Name].index(rna_to_protein_dict[rna])])
            axs[1, p].scatter(np.log10(x+1), np.log10(y+1), c='w', edgecolor='r', alpha=0.7, s=7)

            if x > 0:
                s1_s2_delta = round(y / x, 6)
            else:
                s1_s2_delta = np.nan

            row = [rna.strip('_RNA'), s1_s2_delta, round(x, 6), round(y, 6), key]
            prot_delta = np.vstack((prot_delta, row))
            plotdata = np.vstack((plotdata, [x,y]))

        # r2 to write to tsv
        corr_mat = np.corrcoef(plotdata[:,0], plotdata[:,1])
        prot_r2[key] = corr_mat[0,1]


        equalize_axes(axs[1, p])
        axs[1, p].set_aspect('equal')
        axs[1, p].set_xlabel(sim2Name + ' protein (log)', fontsize=6)
        axs[1, p].set_ylabel(sim1Name + ' protein (log)', fontsize=6)
        axs[1, p].set_title(key, fontsize=8)
        axs[1, p].tick_params(axis='both', labelsize=6)

    plt.suptitle('Average Expression for Functional Genes')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    OutDir = args['sim1Dir'].split('simOut')[0] + 'plotOut/'
    plt.savefig(OutDir + sim1Name + '_' + sim2Name + '_' + 'functional_genes.png')



    # sort rna delta data and save as tsv
    rna_delta = rna_delta[rd[:,0].argsort()]
    with open(OutDir + 'functional_genes_rna_diff.tsv', 'wt') as f:
        tsv_writer = csv.writer(f, delimiter='\t')
        tsv_writer.writerow([sim2Name + ' rnaId', sim1Name + ' rnaId', 'percentage change rna exp', sim2Name + ' sim value', sim1Name + ' sim value','process'])
        for row in rna_delta:
            tsv_writer.writerow(row)
    # write protein data and save to tsv
    prot_delta = prot_delta[np.array(prot_delta[:,1], dtype=float).argsort()]
    with open(OutDir + 'functional_genes_protein_diff.tsv', 'wt') as f:
        tsv_writer = csv.writer(f, delimiter='\t')
        tsv_writer.writerow(['rnaId', 'percentage change protein exp', sim2Name + ' sim value', sim1Name + ' sim value','process'])
        for row in prot_delta:
            tsv_writer.writerow(row)

    # write r2 to tsv
    with open(OutDir + 'functional_genes_r_squared.tsv', 'wt') as f:
        tsv_writer = csv.writer(f, delimiter='\t')
        tsv_writer.writerow(['process', 'rna r2', 'prot r2'])
        for key in keys_to_plot:
            tsv_writer.writerow([key, rna_r2[key], prot_r2[key]])


def equalize_axes(ax):
    lims = ax.get_ylim() + ax.get_xlim()
    max_lim = max(lims)
    min_lim = min(lims)
    ax.set_ylim((min_lim, max_lim))
    ax.set_xlim((min_lim, max_lim))

def build_rna_to_protein_dict():
    protein_tsv = 'reconstruction/ecoli/flat/proteins.tsv'
    DIALECT = "excel-tab"
    JsonReader = partial(spreadsheets.JsonReader, dialect=DIALECT)

    tsv_list = []
    with open(protein_tsv) as tsvfile:
        reader = JsonReader(tsvfile)
        fieldnames = reader.fieldnames
        for row in reader:
            tsv_list.append(row)

    rna_to_protein_dict = {}
    for row in tsv_list:
        rna_to_protein_dict[row['rnaId']] = row['id'] + '[' + row['location'][0] + ']'

    return rna_to_protein_dict



# -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim 1 name', type=str, help='name of first sim')
    parser.add_argument('sim1Dir', type=str, help='directory containing sim (str)')
    parser.add_argument('sim 2 name', type=str, help='name of second sim')
    parser.add_argument('sim2Dir', type=str, help='directory containing master sim (str)')
    parser.add_argument('polycistron file', type=str, help='polycistornic_mrnas_in_model.tsv for sim 1')
    args = vars(parser.parse_args())
    inspect_functional_genes(args['sim 1 name'], args['sim1Dir'], args['sim 2 name'], args['sim2Dir'], args['polycistron file'])
