This directory contains files with data on environmental conditions. There is a corresponding directory ```reconstruction/ecoli/flat/condition/``` that includes files with data on how e. coli responds to some of these conditions.
 
In particular, the media ```minimal```, ```minimal_plus_amino_acids```, ```with_aa```, correspond to data in ```reconstruction/ecoli/flat/condition/```. That corresponding naming convention needs to be maintained for proper model functioning.

# example use of make_media

create media object:
> media_obj = Media()

retrieve stock media:
> base_media = media_obj.stock_media['M9_GLC']
> base_media2 = media_obj.stock_media['5X_supplement_EZ']

define a list of ingredients. Each ingredient is a tuple with (mol_id, weight (units.g), volume (units.L)). If weight is Infinity, it sets the final concentration to infinity. If weight is -Infinity, it sets the final concentration to 0.
> ingredients = [
	('L-ALPHA-ALANINE', 1.78 * units.g, 0.025 * units.L),
	('ARG', 8.44 * units.g, 0.1 * units.L),
	('LEU', float("inf") * units.g, 0 * units.L),
	('OXYGEN-MOLECULE', float("-inf") * units.g, 0 * units.L),
	]

add ingredients directly into an existing media:
> new_media1 = media_obj.add_ingredients(base_media, 0.8 * units.L, ingredients)

combine two medias:
> new_media2 = media_obj.combine_media(base_media, 0.8 * units.L, base_media2, 0.2 * units.L)