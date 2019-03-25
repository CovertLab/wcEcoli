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


class AddIngredientsError(Exception):
	pass

class Media(object):
	''' Media '''

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

	def _dilute_media(self, media_1, media_1_volume, new_volume):
		new_media = {mol_id: 0 for mol_id, conc in media_1.iteritems()}

		# update media concentrations
		for mol_id, conc in media_1.iteritems():
			if conc.asNumber() == float("inf"):
				new_media[mol_id] = float("inf") * CONC_UNITS
			else:
				counts = conc * media_1_volume
				dilute_conc = counts / new_volume
				new_media[mol_id] = dilute_conc

		return new_media

	def combine_media(self, media_1, media_1_volume, media_2, media_2_volume):
		'''
		Combines two medias and returns a new media

		Args:
			media_1, media_2 (dict): dicts with {molecule_id: concentration}
			media_1_volume, media_2_volume: the volume of media_1 and media_2

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		# get new_media volume
		new_volume = media_1_volume + media_2_volume

		# initialize new_media by diluting media_1 to new_volume
		new_media = self._dilute_media(media_1, media_1_volume, new_volume)

		# add media_2
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
		Combines ingredients to existing media. Ingredients are specified as a list,
		with a tuple for each added ingredient specifying (mol_id, weight, volume).
		If weight is Infinity, the concentration is set to inf. If the weight is -Infinity,
		the concentration is set to 0.

		Args:
			media_1 (dict): {molecule_id: concentrations}
			media_1_volume:
			ingredients (list of tuples): [(mol_id, weight, volume)]

		Returns:
			new_media (dict): {molecule_id: concentrations}
		'''

		# get new_media volume
		ingredients_volume = 0 * VOLUME_UNITS
		for (mol_id, weight, volume) in ingredients:
			ingredients_volume += volume
		new_volume = media_1_volume + ingredients_volume

		if new_volume.asNumber() < 0:
			raise AddIngredientsError(
				"Adding negative volume"
			)

		# initialize new_media by diluting media_1 to new_volume
		new_media = self._dilute_media(media_1, media_1_volume, new_volume)

		# add ingredients
		for (mol_id, weight, volume) in ingredients:
			# add infinite concentration of ingredient if weight is Infinity
			if weight.asNumber() == float("inf"):
				new_media[mol_id] = float("inf") * CONC_UNITS

			# remove ingredient from media if weight is -Infinity
			elif weight.asNumber() == float("-inf"):
				new_media[mol_id] = 0.0 * CONC_UNITS

			# if a weight is specified, it needs a formula weight listed in environment_molecules.tsv
			elif weight.asNumber() > 0:
				if self.environment_molecules_fw[mol_id] is not None:
					fw = self.environment_molecules_fw[mol_id]
					added_counts = weight / fw
					new_conc = added_counts / new_volume
					new_media[mol_id] = new_conc
				else:
					raise AddIngredientsError(
						"No fw defined for {} in environment_molecules.tsv".format(mol_id)
					)

			elif weight.asNumber() < 0:
				raise AddIngredientsError(
					"Negative weight for {}".format(mol_id)
				)

			else:
				raise AddIngredientsError(
					"No added weight for {}".format(mol_id)
				)

		return new_media

# # example use
# media_obj = Media()
#
# # make ingredients into a new media
# base_media = media_obj.stock_media['M9_GLC']
# base_media2 = media_obj.stock_media['5X_supplement_EZ']
# ingredients = [
# 	('L-ALPHA-ALANINE', 1.78 * units.g, 0.025 * units.L),
# 	('ARG', 8.44 * units.g, 0.1 * units.L),
# 	('LEU', float("inf") * units.g, 0 * units.L),
# 	('OXYGEN-MOLECULE', float("-inf") * units.g, 0 * units.L),
# 	]
#
# # add ingredients directly into an existing media
# new_media1 = media_obj.add_ingredients(base_media, 0.8 * units.L, ingredients)
#
# # combine two medias
# new_media2 = media_obj.combine_media(base_media, 0.8 * units.L, base_media2, 0.2 * units.L)
