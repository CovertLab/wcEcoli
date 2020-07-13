
import itertools
import os
import csv

currentPath = os.getcwd()

# --- Parse stringList

def parseStringList(string):
    newStrArray = []
    splitStr = string.split('//')
    for st in splitStr:
        newStrArray.append(st.strip())
    return newStrArray

def get_ids(genes):
    object_ids = []
    accession_numbers_1 = []
    accession_numbers_2 = []
    for gene in genes:
        for row in geneData:
            if gene in row[1] or gene in row[2]:
                object_ids.append(row[0])
                accession_numbers_1.append(row[3])
                accession_numbers_2.append(row[4])
        
    return object_ids, accession_numbers_1, accession_numbers_2

# --- Upload TU structure data

tuFile = 'tus_coordinates_directions.txt'
tuPath = currentPath + '/' +  tuFile

tuData = list(csv.reader(open(tuPath, 'r'), delimiter='\t'))[1:]
tuData_header = list(csv.reader(open(tuPath, 'r'), delimiter='\t'))[0]

# --- Upload gene data

geneFile = '070320_gene_ids_syns.txt'
genePath = currentPath + '/' +  geneFile

geneData = list(csv.reader(open(genePath, 'r'), delimiter = '\t'))



#get rid of duplicate rows in tuData

tuData.sort()
tuData = list(tuData for tuData,_ in itertools.groupby(tuData))

tu_gene_ids = []
genes_in_tu = []
all_tu_object_ids = []
all_tu_accession_nums_1 = []
all_tu_accession_nums_2 = []
for row in tuData:
    genes = parseStringList(row[0])
    genes_in_tu.append(genes)
    object_ids, accession_numbers_1, accession_numbers_2 = get_ids(genes)
    all_tu_object_ids.append(object_ids)
    tu_gene_id = '_'.join(object_ids)

tu_obj_ids_filtered = list(filter(None, all_tu_object_ids))

polycistrons_file = []
header = ["transcription_units", "monomers_to_remove"]
polycistrons_file.append(header)
for row in tu_obj_ids_filtered:
    if len(row) > 1:
        polycistrons_file.append([row, row])

with open('polycistrons_file.csv', 'wb') as myfile:
    wr = csv.writer(myfile, delimiter=',', escapechar=' ', quoting=csv.QUOTE_NONE)
    wr.writerows(polycistrons_file)