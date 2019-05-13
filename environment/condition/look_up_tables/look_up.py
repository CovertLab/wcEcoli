from __future__ import absolute_import, division, print_function

import os

from reconstruction.spreadsheets import JsonReader

LOOKUP_DIR = os.path.join('environment', 'condition', 'look_up_tables')

LIST_CONC_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "avg_concentrations", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "avg_concentrations", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "avg_concentrations", "minimal_plus_amino_acids.tsv"),
)
LIST_FLUX_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "avg_flux", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "avg_flux", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "avg_flux", "minimal_plus_amino_acids.tsv"),
)



class LookUp(object):

	def __init__(self):

		# Load saved wcEcoli average concentrations from all media conditions
		self.avg_conc_lookup = {}
		for file_name in LIST_CONC_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.avg_conc_lookup[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					molecule_id = row["molecule id"]
					conc = row["average concentration mmol/L"]
					self.avg_conc_lookup[media][molecule_id] = conc

		# Load saved wcEcoli average fluxes from all media conditions
		self.avg_flux_lookup = {}
		for file_name in LIST_FLUX_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.avg_flux_lookup[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					reaction_id = row["reaction id"]
					flux = row["average flux mmol/L"]
					self.avg_flux_lookup[media][reaction_id] = flux