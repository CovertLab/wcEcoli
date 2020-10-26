from reconstruction.spreadsheets import tsv_reader


def TU_to_operons(polycistron_file):
    with tsv_reader(polycistron_file) as reader:
        tsv_fieldnames = reader.fieldnames
        polycistrons = list(reader)




