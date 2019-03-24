'''
Functions for making media
'''

from __future__ import absolute_import, division, print_function

from wholecell.utils import units

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS


class MakeMedia(object):

	def __init__(self):

		raw_data = KnowledgeBaseEcoli()

		# get dicts from knowledge base
		self.environment_molecules_fw = self.get_environment_molecules_fw(raw_data)
		self.stock_media = self.get_stock_media(raw_data)

	def get_environment_molecules_fw(self, raw_data):
		'''get formula weight (units.g / units.mol) for all environmental molecules'''

		environment_molecules_fw = {}
		for row in raw_data.condition.environment_molecules:
			mol = row["molecule id"]
			fw = row["formula weight"]
			if fw == 'None':
				environment_molecules_fw[mol] = None
			else:
				environment_molecules_fw[mol] = float(fw) * (units.g / units.mol)

		return environment_molecules_fw

	def get_stock_media(self, raw_data):
		stock_media = {}
		for label in vars(raw_data.condition.media):
			# initiate all molecules with 0 concentrations
			stock_media[label] = {
				row["molecule id"]: 0.0 * CONC_UNITS
				for row in raw_data.condition.environment_molecules}

			# get non-zero concentrations (assuming units.mmol / units.L)
			molecule_concentrations = getattr(raw_data.condition.media, label)

			environment_non_zero_dict = {
				row["molecule id"]: row["concentration"]
				for row in molecule_concentrations}

			# update environment_dict with non zero concentrations
			stock_media[label].update(environment_non_zero_dict)

		return stock_media

	def combine_media(self, media_1, media_1_volume, media_2, media_2_volume):
		'''
		combines two medias and returns a new media

		Args:
			media_1, media_2 (dict): a dict with {molecule_id: concentration}
			media_1_volume, media_2_volume: the volume of media combined into new_media

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		# initialize new media
		new_media = {mol_id: 0 for mol_id, conc in media_1.iteritems()}
		new_volume = media_1_volume + media_2_volume

		# update media concentrations
		for mol_id, conc in media_1.iteritems():
			if conc.asNumber() == float("inf"):
				new_media[mol_id] = float("inf") * CONC_UNITS
			else:
				counts = conc * media_1_volume
				dillute_conc = counts / new_volume
				new_media[mol_id] = dillute_conc

		# update media concentrations
		for mol_id, conc in media_2.iteritems():
			old_conc = new_media[mol_id]

			if conc.asNumber() == float("inf") or old_conc.asNumber() == float("inf"):
				new_media[mol_id] = float("inf") * CONC_UNITS
			else:
				old_counts = old_conc * new_volume
				added_counts = conc * media_2_volume
				new_counts = old_counts + added_counts
				new_conc = new_counts / new_volume

				# update media
				new_media[mol_id] = new_conc

		return new_media

	def add_ingredients(self, media_1, media_1_volume, ingredients):
		'''

		Args:
			media_1 (dict): {molecule_id: concentrations}
			media_1_volume:
			ingredients (tuple): (mol_id, weight, volume)

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		# initialize new media
		new_media = {mol_id: 0 for mol_id, conc in media_1.iteritems()}

		# get new media volume
		ingredients_volume = 0 * VOLUME_UNITS
		for (mol_id, weight, volume) in ingredients:
			ingredients_volume += volume
		new_volume = media_1_volume + ingredients_volume

		# update media concentrations
		for mol_id, conc in media_1.iteritems():
			if conc.asNumber() == float("inf"):
				new_media[mol_id] = float("inf") * CONC_UNITS
			else:
				counts = conc * media_1_volume
				new_conc = counts / new_volume
				new_media[mol_id] = new_conc

		# add ingredients
		for (mol_id, weight, volume) in ingredients:
			old_conc = new_media[mol_id]
			old_counts = old_conc * new_volume

			# get counts of added ingredient, add to media
			fw = self.environment_molecules_fw[mol_id]
			added_counts = weight / fw
			new_counts = old_counts + added_counts
			new_conc = new_counts / new_volume

			new_media[mol_id] = new_conc

		return new_media

# examples
media_obj = MakeMedia()
base_media = media_obj.stock_media['M9_GLC']
ingredients = [('L-ALPHA-ALANINE', 1.78 * units.g, 0.025 * units.L)]
new_media = media_obj.add_ingredients(base_media, 0.8 * units.L, ingredients)

media1 = media_obj.stock_media['M9_GLC']
media2 = media_obj.stock_media['5X_supplement_EZ']
new_media2 = media_obj.combine_media(media1, 0.8 * units.L, media2, 0.2 * units.L)
