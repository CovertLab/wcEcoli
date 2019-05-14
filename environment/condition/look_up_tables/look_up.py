from __future__ import absolute_import, division, print_function

import os
import csv
import random

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

					# convert to list of floats
					conc_dist = conc_dist.replace('[', '').replace(']', '').split(', ')
					conc_dist = [float(conc) for conc in conc_dist]

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

					# convert to list of floats
					flux_dist = flux_dist.replace('[', '').replace(']', '').split(', ')
					flux_dist = [float(flux) for flux in flux_dist]

					self.flux_avg[media][reaction_id] = flux_avg
					self.flux_dist[media][reaction_id] = flux_dist


	def get_fluxes(self, lookup_type, media, reaction_ids):
		''' Get a flux for each reaction in reaction_ids'''

		fluxes = {}
		if lookup_type == 'average':
			fluxes = {
				reaction_id: self.flux_avg[media][reaction_id]
				for reaction_id in reaction_ids}
		if lookup_type == 'distribution':
			fluxes = {
				reaction_id: random.choice(self.flux_dist[media][reaction_id])
				for reaction_id in reaction_ids}

		return fluxes

	def get_concs(self, lookup_type, media, molecule_ids):
		''' Get a flux for each reaction in reaction_ids'''

		concentrations = {}
		if lookup_type == 'average':
			concentrations = {
				molecule_id: self.concentration_avg[media][molecule_id]
				for molecule_id in molecule_ids}
		if lookup_type == 'distribution':
			concentrations = {
				molecule_id: random.choice(self.concentration_dist[media][molecule_id])
				for molecule_id in molecule_ids}

		return concentrations
