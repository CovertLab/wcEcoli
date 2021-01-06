from reconstruction.spreadsheets import tsv_reader
import csv
from collections import defaultdict
import copy
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def rna_seq_diff(operon_rnas):
    # READ IN FLAT FILES
    gene_data_file = 'reconstruction/ecoli/flat/genes.tsv'
    with tsv_reader(gene_data_file) as reader:
        gene_data = list(reader)

    # gene id to coordinate, for ordering genes for plotting
    id_to_coord = {x['id']: x['coordinate'] for x in gene_data}
    id_to_dir = {x['id']: x['direction'] for x in gene_data}

    with tsv_reader(operon_rnas) as reader:
        tsv_data = list(reader)

    # make dictionaries to map rnaid <> geneid
    # make list of polycistronic rnas
    geneid_to_rna_names = defaultdict(list)
    rna_name_to_geneids = defaultdict(list)
    polycistron_list = []

    for rna in tsv_data:
        if rna['id'].count('_') > 1:
            polycistron_list.append(rna['id'])
        rna_name_to_geneids[rna['id']] = rna['gene_set']
        for gene in rna['gene_set']:
            geneid_to_rna_names[gene].append(rna['id'])

    # sort by operon length
    polycistron_list= sorted(polycistron_list, key=lambda x: x.count("_"), reverse=True)
    rnas_included = []
    operons = defaultdict(list)
    # gather all rnas that cover a set of genes

    for pc in polycistron_list:
        if pc not in rnas_included:
            for gene in rna_name_to_geneids[pc]:
                for rna in geneid_to_rna_names[gene]:
                    if rna not in operons[pc]:
                        # organize rnas into operons dictionary
                        operons[pc].append(rna)
                        rnas_included.append(rna)
            # import ipdb; ipdb.set_trace()
            operons[pc].sort(key=lambda x: x.count('_'), reverse=True)

    # get all genes contained in an operon
    operon_to_genes = {}
    for operon_rna in operons:
        gene_list = []
        for rna in operons[operon_rna]:
            gene_list = gene_list + rna_name_to_geneids[rna]

        operon_to_genes[operon_rna] = np.unique(gene_list)
        # sort the genes by their position
        try:
            dir = id_to_dir[operon_to_genes[operon_rna][0]]
            coord_list = [id_to_coord[x] for x in operon_to_genes[operon_rna]]
            if dir == '+':
                sorted_index = np.argsort(coord_list)
            elif dir == '-':
                sorted_index = np.argsort(coord_list)[::-1]
            operon_to_genes[operon_rna] = operon_to_genes[operon_rna][sorted_index]
        except:
            import ipdb; ipdb.set_trace()


        # import ipdb; ipdb.set_trace()


    rna_seq_file = 'reconstruction/ecoli/flat/rna_seq_data/rnaseq_rsem_tpm_mean.tsv'
    with tsv_reader(rna_seq_file) as reader:
        tsv_data = list(reader)

    gene_to_count = {}
    for row in tsv_data:
        gene_to_count[row['Gene']] = row['M9 Glucose minus AAs']

    gene_to_count_raw = copy.deepcopy(gene_to_count)

    operon_count_disparity = {}
    operon_count_disparity_norm = {}
    operon_count = {}
    disparity_mean = np.empty(len(operons.keys()), dtype=float)



    for idx, operon_rna in enumerate(operons.keys()):
        operon_count_disparity[operon_rna] = {x : gene_to_count[x] for x in operon_to_genes[operon_rna]}
        operon_count[operon_rna] = {x : gene_to_count[x] for x in operon_to_genes[operon_rna]}
        TUs = sorted(operons[operon_rna], key=lambda x: x.count('_'), reverse=True)

        for rna in TUs:
            genes_in_rna = rna.split('_')[0:-1]
            min_count = min([gene_to_count[x] for x in genes_in_rna])
            if min_count < 0:
                min_count = 0

            for gene in genes_in_rna:
                operon_count_disparity[operon_rna][gene] -= min_count
                gene_to_count[gene] -= min_count

        # TODO: make sure this normalization scheme makes sense
        count_norm = min(list(operon_count[operon_rna].values())) + 1
        # count_norm = mean(list(operon_count[operon_rna].values()))
        operon_count_disparity_norm[operon_rna] = {k: abs(v)/count_norm for k, v in operon_count_disparity[operon_rna].items()}
        disparity_mean[idx] = np.mean(list(operon_count_disparity_norm[operon_rna].values()))


    # DEG RATE PAPER
    # compare this data to mrna deg rate paper
    rna_deg_file = 'prototypes/operons/Dar_Soreck_RNA_deg.csv'
    with open(rna_deg_file, newline='', encoding='utf-8-sig') as csv_file:
        reader = csv.reader(csv_file)
        rna_deg_data = [line for line in reader]


    name_to_id = {x['symbol']: x['id'] for x in gene_data}
    # Alternate names used in this dataset mapped to ecocyc ID
    name_to_id['bepA'] = 'G7311'
    name_to_id['rfc'] = 'EG11982'
    name_to_id['yfbB'] = 'EG12438'
    name_to_id['mraW'] = 'EG11085'
    name_to_id['rfaF'] = 'EG12210'
    name_to_id['rfaC'] = 'EG11189'
    name_to_id['rfaL'] = 'EG11424'
    name_to_id['yjeE'] = 'EG11757'


    rna_deg_dict = {}
    covert_ratio_list = []
    dar_ratio_list = []
    stable_gene_list = []
    for line in rna_deg_data[1:]:
        if line[0] == '':
            break
        name = line[0]
        rna_deg_dict[name] = {}

        # convert operon name to ecocyc ids
        operon_gene_ids = []
        for x in line[0].split('-'):
            operon_gene_ids.append(name_to_id[x])
        # store info as ecocyc ids
        rna_deg_dict[name]['operon ids'] = copy.deepcopy(operon_gene_ids)
        rna_deg_dict[name]['upstream'] = name_to_id[line[1]]
        rna_deg_dict[name]['downstream'] = name_to_id[line[2]]
        rna_deg_dict[name]['stabilized gene'] = rna_deg_dict[name][line[3]]

        rna_counts = []
        operon_gene_ids.remove(rna_deg_dict[name]['stabilized gene'])
        for x in operon_gene_ids:
            rna_counts.append(gene_to_count_raw[x])

        rna_ratio = gene_to_count_raw[rna_deg_dict[name]['stabilized gene']]/np.mean(rna_counts)

        rna_deg_dict[name]['Covert ratio'] = rna_ratio

        rna_deg_dict[name]['Dar Soreck ratio'] = float(line[8])

        rna_deg_dict[name]['ratio diff'] = abs(rna_ratio - float(line[8]))

        # lists for plotting
        covert_ratio_list.append(rna_ratio)
        dar_ratio_list.append(float(line[8]))
        stable_gene_list.append(rna_deg_dict[name]['stabilized gene'])




    # PLOT DAR VS COVERT RNA RATIO
    m, b = np.polyfit(dar_ratio_list, covert_ratio_list, 1)
    y = m*np.array(dar_ratio_list) + b
    plt.figure(dpi=300)
    plt.scatter(dar_ratio_list, covert_ratio_list, color='w', edgecolor='k', alpha=0.6)
    plt.plot(dar_ratio_list, y , color='r', linestyle='-.')
    plt.xlabel('Dar Soreck 2018 mRNA Ratio')
    plt.ylabel('Covert RNAseq mRNA Ratio')
    plt.savefig('prototypes/operons/dar_covert_ratio_corr.png')

    with open('prototypes/operons/dar_soreck_comparison.csv', 'wt') as f:
        tsv_writer = csv.writer(f, delimiter = '\t')
        tsv_writer.writerow(['operon', 'stabilized gene', 'Dar ratio', 'Covert ratio', 'ratio diff'])

        for k in rna_deg_dict.keys():
            row = [rna_deg_dict[k]['operon ids'], rna_deg_dict[k]['stabilized gene'],
                   rna_deg_dict[k]['Dar Soreck ratio'], round(rna_deg_dict[k]['Covert ratio'],2),
                   round(rna_deg_dict[k]['ratio diff'],2)]
            tsv_writer.writerow(row)

    # RNA SEQ DISPARITY
    # -----------------------------------------------------------------------------------------------




    sorted_index = np.argsort(disparity_mean)
    disparity_mean = disparity_mean[sorted_index]
    operon_list = np.array(list(operons.keys()))[sorted_index]

    num_pages = int(np.ceil(len(operon_list)/4))
    plot_total = 0
    with PdfPages('prototypes/operons/operon_rna_seq_graphs.pdf') as pdf:
        for page in range(0,num_pages):

            while plot_total < len(operon_list):
                fig, ax = plt.subplots(2, 2)
                ax = ax.flatten()
                for p in range(0,4):
                    rna = operon_list[plot_total]

                    # add S if rna is stabilized
                    xlabels = list(operon_count[rna].keys())
                    xlables_copy = copy.deepcopy(xlabels)
                    for idx, x in enumerate(xlabels):
                        if x in stable_gene_list:
                            xlabels[idx] = x + ' (S)'

                    # bar graph
                    x = np.arange(len(operon_to_genes[rna]))
                    y = list(operon_count[rna].values())
                    ax[p].bar(x,y, color='darkgrey')
                    ax[p].set_xticks(x)
                    ax[p].set_xticklabels(xlabels, rotation=45, fontsize=8)
                    y2 = [abs(x) for x in operon_count_disparity[rna].values()]
                    ax[p].bar(x, y2, alpha=0.4, color='r')
                    ax[p].set_ylabel('rna seq count')

                    # plot lines to represent operon structure
                    # import ipdb; ipdb.set_trace()

                    max_ylim = ax[p].get_ylim()[1]

                    step = max_ylim * 0.05
                    line_y = max_ylim + step * 5
                    # ax[p].set_prop_cycle('color',[plt.cm.bone(i) for i in np.linspace(0, 1, 6)])
                    for rna_line in operons[rna]:

                        start = xlables_copy.index(rna_line.split('_')[0]) - 0.25
                        end = xlables_copy.index(rna_line.split('_')[-2]) + 0.25

                        ax[p].plot([start, end], [line_y, line_y])

                        line_y += step


                    plot_total += 1

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

    import ipdb; ipdb.set_trace()
    with open('prototypes/operons/' + 'operon_rna_seq_disparity.tsv', 'wt') as f:
        tsv_writer = csv.writer(f, delimiter = '\t')
        tsv_writer.writerow(['operon name', 'mean count disparity', 'genes in operon'])
        for idx, rna in enumerate(operon_list):
            row = [rna, disparity_mean[idx], operon_to_genes[rna]]
            tsv_writer.writerow(row)





# -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    rna_seq_diff('reconstruction/ecoli/flat/operon_rnas.tsv')
    # count_disparity(operons, operon_to_genes)




