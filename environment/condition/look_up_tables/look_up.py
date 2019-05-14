from __future__ import absolute_import, division, print_function

import os
import csv

from reconstruction.spreadsheets import JsonReader

LOOKUP_DIR = os.path.join('environment', 'condition', 'look_up_tables')

CONC_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "transport_concentrations", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "transport_concentrations", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "transport_concentrations", "minimal_plus_amino_acids.tsv"),
)
FLUX_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "transport_fluxes", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "transport_fluxes", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "transport_fluxes", "minimal_plus_amino_acids.tsv"),
)
KCAT_LOOKUP_FILES = (
	os.path.join(LOOKUP_DIR, "transport_kcats", "minimal.tsv"),
	os.path.join(LOOKUP_DIR, "transport_kcats", "minimal_minus_oxygen.tsv"),
	os.path.join(LOOKUP_DIR, "transport_kcats", "minimal_plus_amino_acids.tsv"),
)

CSV_DIALECT = csv.excel_tab

class LookUp(object):

	def __init__(self):

		# Load saved wcEcoli  concentrations from all media conditions
		self.concentration_avg = {}
		self.concentration_dist = {}
		for file_name in CONC_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.concentration_avg[media] = {}
			self.concentration_dist[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					molecule_id = row.get("enzyme id")
					conc_avg = row.get("concentration avg mmol/L")
					conc_dist = row.get("concentration distribution mmol/L")
					self.concentration_avg[media][molecule_id] = conc_avg
					self.concentration_dist[media][molecule_id] = conc_dist

		# Load saved wcEcoli fluxes from all media conditions
		self.flux_avg = {}
		self.flux_dist = {}
		for file_name in FLUX_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.flux_avg[media] = {}
			self.flux_dist[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					reaction_id = row.get("reaction id")
					flux_avg = row.get("flux avg mmol/L/s")
					flux_dist = row.get("flux distribution mmol/L/s")
					self.flux_avg[media][reaction_id] = flux_avg
					self.flux_dist[media][reaction_id] = flux_dist

		# Load estimated wcEcoli k_cats from all media conditions
		self.kcat_avg = {}
		self.kcat_dist = {}
		for file_name in FLUX_LOOKUP_FILES:
			media = file_name.split(os.path.sep)[-1].split(".")[0]
			self.kcat_avg[media] = {}
			self.kcat_dist[media] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					reaction_id = row.get("reaction id")
					# enzyme_id = row.get("enzyme id")
					kcat_avg = row.get("k_cat avg")
					kcat_dist = row.get("k_cat distribution")
					self.kcat_avg[media][reaction_id] = kcat_avg
					self.kcat_dist[media][reaction_id] = kcat_dist
