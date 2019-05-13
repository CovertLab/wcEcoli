from __future__ import absolute_import, division, print_function

import os
import csv

from reconstruction.spreadsheets import JsonReader

LOOKUP_DIR = os.path.join('environment', 'condition', 'look_up_tables')

AVG_CONC_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "avg_concentrations", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "avg_concentrations", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "avg_concentrations", "minimal_plus_amino_acids.tsv"),
)
AVG_FLUX_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "avg_flux", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "avg_flux", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "avg_flux", "minimal_plus_amino_acids.tsv"),
)
FLUX_DISTRIBUTION_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "flux_distribution", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "flux_distribution", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "flux_distribution", "minimal_plus_amino_acids.tsv"),
)

CSV_DIALECT = csv.excel_tab

class LookUp(object):

	def __init__(self):

		# Load saved wcEcoli average concentrations from all media conditions
		self.avg_concentration = {}
		for file_name in AVG_CONC_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.avg_concentration[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					molecule_id = row["molecule id"]
					conc = row["average concentration mmol/L"]
					self.avg_concentration[media][molecule_id] = conc

		# Load saved wcEcoli average fluxes from all media conditions
		self.avg_flux = {}
		for file_name in AVG_FLUX_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.avg_flux[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					reaction_id = row["reaction id"]
					flux = row["average flux mmol/L"]
					self.avg_flux[media][reaction_id] = flux

		# Load saved wcEcoli average fluxes from all media conditions
		self.flux_distribution = {}
		for file_name in FLUX_DISTRIBUTION_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.flux_distribution[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					reaction_id = row["reaction id"]
					flux = row["average flux mmol/L"]
					self.flux_distribution[media][reaction_id] = flux