import csv
from reconstruction import spreadsheets
from functools import partial


DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

with open('reconstruction/ecoli/flat/proteins.tsv') as f:
    tsvreader = csv.reader(f, delimiter='\t')

    prot_to_gene = {}
    lcount = 0
    for line in tsvreader:

        if lcount > 0:
            prot_to_gene[line[0]] = line[5]

        lcount += 1


with open('runscripts/reflect/functional_genes.tsv') as f:

    tsvreader = csv.reader(f, delimiter='\t')

    functional_gene_list = []

    lcount = 0
    for line in tsvreader:

        if lcount > 1:
            functional_gene_list.append(prot_to_gene[line[2]])

        lcount += 1


with open('prototypes/operons/all_ecocyc_tus/all_polycistrons.tsv') as f:
    reader = JsonReader(f)
    all_pcs = []
    for row in reader:
        all_pcs.append(list(row.values()))


fieldnames = ["transcription_units", "monomers_to_remove"]
with open('prototypes/operons/all_ecocyc_tus/no_functional_genes_polycistron_file.tsv', 'w') as f:
    writer = JsonWriter(f, fieldnames, dialect=csv.excel_tab)
    writer.writeheader()
    for line in all_pcs:
        if not any(item in line[0] for item in functional_gene_list):
            writer.writerow(dict(zip(fieldnames, line)))
        else:
            print(line[0])


