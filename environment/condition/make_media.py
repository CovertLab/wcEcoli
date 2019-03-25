'''
Functions for making media

# example use
media_obj = MakeMedia()

# make ingredients into a new media
base_media = media_obj.stock_media['M9_GLC']
base_media2 = media_obj.stock_media['5X_supplement_EZ']
ingredients = [
	('L-ALPHA-ALANINE', 1.78 * units.g, 0.025 * units.L),
	('ARG', 8.44 * units.g, 0.1 * units.L),
]

new_media1, new_volume = media_obj.ingredients_to_media(ingredients)

# add previously mixed ingredients into a different media
new_media2 = media_obj.combine_media(base_media, 0.8 * units.L, new_media1, new_volume)

# add ingredients directly into an existing media
new_media3 = media_obj.add_ingredients(base_media, 0.8 * units.L, ingredients)

# combine two medias
new_media4 = media_obj.combine_media(base_media, 0.8 * units.L, base_media2, 0.2 * units.L)
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
	''' MakeMedia '''

	def __init__(self):

		raw_data = KnowledgeBaseEcoli()

		# get dicts from knowledge base
		self.environment_molecules_fw = self._get_environment_molecules_fw(raw_data)
		self.stock_media = self._get_stock_media(raw_data)

	def _get_environment_molecules_fw(self, raw_data):
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

	def _get_stock_media(self, raw_data):
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

	def ingredients_to_media(self, ingredients):
		'''
		Args:
			ingredients (tuple): (mol_id, weight, volume)

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		new_media = {mol_id: 0 * CONC_UNITS for mol_id, fws in self.environment_molecules_fw.iteritems()}

		# get total volume
		total_volume = 0 * VOLUME_UNITS
		for (mol_id, weight, volume) in ingredients:
			total_volume += volume

		if total_volume.asNumber() == 0:
			no_volume = True

		# TODO -- handle infinite weight, 0 volume,

		# add ingredients
		for (mol_id, weight, volume) in ingredients:
			if mol_id in self.environment_molecules_fw:
				fw = self.environment_molecules_fw[mol_id]
				added_counts = weight / fw
				new_conc = added_counts / total_volume
				new_media[mol_id] = new_conc

			# else:
			# 	# TODO -- print exception 'fw not defined'
			# 	continue

		return new_media, total_volume

	def add_ingredients(self, media_1, media_1_volume, ingredients):
		'''

		Args:
			media_1 (dict): {molecule_id: concentrations}
			media_1_volume:
			ingredients (tuple): (mol_id, weight, volume)

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		ingredients_media, volume = self.ingredients_to_media(ingredients)
		new_media = self.combine_media(media_1, media_1_volume, ingredients_media, volume)
		return new_media
