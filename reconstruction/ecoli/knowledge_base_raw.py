"""
KnowledgeBase for Ecoli
Whole-cell knowledge base for Ecoli. Contains all raw, un-fit data processed
directly from CSV flat files.

"""
import io
import os
import json
from typing import List, Dict
import warnings

from reconstruction.spreadsheets import read_tsv
from wholecell.io import tsv
from wholecell.utils import units  # used by eval()

FLAT_DIR = os.path.join(os.path.dirname(__file__), "flat")
LIST_OF_DICT_FILENAMES = [
	"amino_acid_export_kms.tsv",
	"amino_acid_export_kms_removed.tsv",
	"amino_acid_pathways.tsv",
	"amino_acid_uptake_rates.tsv",
	"amino_acid_uptake_rates_removed.tsv",
	"biomass.tsv",
	"compartments.tsv",
	"complexation_reactions.tsv",
	"complexation_reactions_added.tsv",
	"complexation_reactions_modified.tsv",
	"complexation_reactions_removed.tsv",
	"disabled_kinetic_reactions.tsv",
	"dna_sites.tsv",
	"dry_mass_composition.tsv",
	"endoRNases.tsv",
	"equilibrium_reaction_rates.tsv",
	"equilibrium_reactions.tsv",
	"equilibrium_reactions_added.tsv",
	"equilibrium_reactions_removed.tsv",
	"fold_changes.tsv",
	"fold_changes_nca.tsv",
	"fold_changes_removed.tsv",
	"footprint_sizes.tsv",
	"genes.tsv",
	"growth_rate_dependent_parameters.tsv",
	"linked_metabolites.tsv",
	"metabolic_reactions.tsv",
	"metabolic_reactions_added.tsv",
	"metabolic_reactions_modified.tsv",
	"metabolic_reactions_removed.tsv",
	"metabolism_kinetics.tsv",
	"metabolite_concentrations.tsv",
	"metabolite_concentrations_removed.tsv",
	"metabolites.tsv",
	"metabolites_added.tsv",
	"modified_proteins.tsv",
	"molecular_weight_keys.tsv",
	"ppgpp_fc.tsv",
	"ppgpp_regulation.tsv",
	"ppgpp_regulation_added.tsv",
	"ppgpp_regulation_removed.tsv",
	"protein_half_lives_measured.tsv",
	"protein_half_lives_n_end_rule.tsv",
	"protein_half_lives_pulsed_silac.tsv",
	"proteins.tsv",
	"relative_metabolite_concentrations.tsv",
	"rna_half_lives.tsv",
	"rna_maturation_enzymes.tsv",
	"rnas.tsv",
	"secretions.tsv",
	"sequence_motifs.tsv",
	"transcription_factors.tsv",
	# "transcription_units.tsv",  # special cased in the constructor
	"transcription_units_added.tsv",
	"transcription_units_removed.tsv",
	"transcriptional_attenuation.tsv",
	"transcriptional_attenuation_removed.tsv",
	"tf_one_component_bound.tsv",
	"translation_efficiency.tsv",
	"trna_charging_reactions.tsv",
	"trna_charging_reactions_added.tsv",
	"trna_charging_reactions_removed.tsv",
	"two_component_systems.tsv",
	"two_component_system_templates.tsv",
	os.path.join("mass_fractions", "glycogen_fractions.tsv"),
	os.path.join("mass_fractions", "ion_fractions.tsv"),
	os.path.join("mass_fractions", "LPS_fractions.tsv"),
	os.path.join("mass_fractions", "lipid_fractions.tsv"),
	os.path.join("mass_fractions", "murein_fractions.tsv"),
	os.path.join("mass_fractions", "soluble_fractions.tsv"),
	os.path.join("trna_data", "trna_ratio_to_16SrRNA_0p4.tsv"),
	os.path.join("trna_data", "trna_ratio_to_16SrRNA_0p7.tsv"),
	os.path.join("trna_data", "trna_ratio_to_16SrRNA_1p6.tsv"),
	os.path.join("trna_data", "trna_ratio_to_16SrRNA_1p07.tsv"),
	os.path.join("trna_data", "trna_ratio_to_16SrRNA_2p5.tsv"),
	os.path.join("trna_data", "trna_growth_rates.tsv"),
	os.path.join("rna_seq_data", "rnaseq_rsem_tpm_mean.tsv"),
	os.path.join("rna_seq_data", "rnaseq_rsem_tpm_std.tsv"),
	os.path.join("rna_seq_data", "rnaseq_seal_rpkm_mean.tsv"),
	os.path.join("rna_seq_data", "rnaseq_seal_rpkm_std.tsv"),
	os.path.join("rrna_options", "remove_rrff", "genes_removed.tsv"),
	os.path.join("rrna_options", "remove_rrff", "rnas_removed.tsv"),
	os.path.join("rrna_options", "remove_rrff", "transcription_units_modified.tsv"),
	os.path.join("rrna_options", "remove_rrna_operons", "transcription_units_removed.tsv"),
	os.path.join("condition", "tf_condition.tsv"),
	os.path.join("condition", "condition_defs.tsv"),
	os.path.join("condition", "environment_molecules.tsv"),
	os.path.join("condition", "timelines_def.tsv"),
	os.path.join("condition", "media_recipes.tsv"),
	os.path.join("condition", "media", "5X_supplement_EZ.tsv"),
	os.path.join("condition", "media", "MIX0-55.tsv"),
	os.path.join("condition", "media", "MIX0-57.tsv"),
	os.path.join("condition", "media", "MIX0-58.tsv"),
	os.path.join("condition", "media", "MIX0-844.tsv"),
	os.path.join("base_codes", "amino_acids.tsv"),
	os.path.join("base_codes", "dntp.tsv"),
	os.path.join("base_codes", "nmp.tsv"),
	os.path.join("base_codes", "ntp.tsv"),
	os.path.join("adjustments", "amino_acid_pathways.tsv"),
	os.path.join("adjustments", "balanced_translation_efficiencies.tsv"),
	os.path.join("adjustments", "translation_efficiencies_adjustments.tsv"),
	os.path.join("adjustments", "rna_expression_adjustments.tsv"),
	os.path.join("adjustments", "rna_deg_rates_adjustments.tsv"),
	os.path.join("adjustments", "protein_deg_rates_adjustments.tsv"),
	os.path.join("adjustments", "relative_metabolite_concentrations_changes.tsv")
	]
SEQUENCE_FILE = 'sequence.fasta'
LIST_OF_PARAMETER_FILENAMES = (
	"dna_supercoiling.tsv",
	"parameters.tsv",
	"mass_parameters.tsv",
	)

REMOVED_DATA = {
	'amino_acid_export_kms': 'amino_acid_export_kms_removed',
	'amino_acid_uptake_rates': 'amino_acid_uptake_rates_removed',
	'complexation_reactions': 'complexation_reactions_removed',
	'equilibrium_reactions': 'equilibrium_reactions_removed',
	'fold_changes': 'fold_changes_removed',
	'fold_changes_nca': 'fold_changes_removed',
	'metabolic_reactions': 'metabolic_reactions_removed',
	'metabolite_concentrations': 'metabolite_concentrations_removed',
	'ppgpp_regulation': 'ppgpp_regulation_removed',
	'transcriptional_attenuation': 'transcriptional_attenuation_removed',
	'trna_charging_reactions': 'trna_charging_reactions_removed',
	}
MODIFIED_DATA = {
	'complexation_reactions': 'complexation_reactions_modified',
	'metabolic_reactions': 'metabolic_reactions_modified',
	}

ADDED_DATA = {
	'complexation_reactions': 'complexation_reactions_added',
	'equilibrium_reactions': 'equilibrium_reactions_added',
	'metabolic_reactions': 'metabolic_reactions_added',
	'metabolites': 'metabolites_added',
	'ppgpp_regulation': 'ppgpp_regulation_added',
	'trna_charging_reactions': 'trna_charging_reactions_added',
	}


class DataStore(object):
	def __init__(self):
		pass

class KnowledgeBaseEcoli(object):
	""" KnowledgeBaseEcoli """

	def __init__(self, operons_on: bool, remove_rrna_operons: bool, remove_rrff: bool, new_genes_option: str="off"):
		self.operons_on = operons_on
		self.new_genes_option = new_genes_option

		if not operons_on and remove_rrna_operons:
			warnings.warn("Setting the 'remove_rrna_operons' option to 'True'"
			              " has no effect on the simulations when the 'operon'"
			              " option is set to 'off'.")

		self.compartments: List[dict] = []  # mypy can't track setattr(self, attr_name, rows)
		self.transcription_units: List[dict] = []

		# Make copies to prevent issues with sticky global variables when
		# running multiple operon workflows through Fireworks
		self.list_of_dict_filenames: List[str] = LIST_OF_DICT_FILENAMES.copy()
		self.removed_data: Dict[str, str] = REMOVED_DATA.copy()
		self.modified_data: Dict[str, str] = MODIFIED_DATA.copy()
		self.added_data: Dict[str, str] = ADDED_DATA.copy()

		self.new_gene_added_data: Dict[str,str] = {}

		if self.operons_on:
			self.list_of_dict_filenames.append('transcription_units.tsv')
			if remove_rrna_operons:
				# Use alternative file with all rRNA transcription units if
				# remove_rrna_operons option was used
				self.removed_data.update({
					'transcription_units': 'rrna_options.remove_rrna_operons.transcription_units_removed',
					})
			else:
				self.removed_data.update({
					'transcription_units': 'transcription_units_removed',
				})
			self.added_data.update({
				'transcription_units': 'transcription_units_added',
				})

		if remove_rrff:
			self.removed_data.update({
				'genes': 'rrna_options.remove_rrff.genes_removed',
				'rnas': 'rrna_options.remove_rrff.rnas_removed',
				})
			if self.operons_on:
				self.modified_data.update({
					'transcription_units': 'rrna_options.remove_rrff.transcription_units_modified',
					})

		if self.new_genes_option != 'off':

			new_gene_subdir = new_genes_option
			new_gene_path = os.path.join('new_gene_data',new_gene_subdir)
			assert os.path.isdir(os.path.join(FLAT_DIR,new_gene_path)), \
				"This new_genes_data subdirectory is invalid."
			nested_attr = 'new_gene_data.' + new_gene_subdir + "."

			# These files do not need to be joined to existing files
			self.list_of_dict_filenames.append(os.path.join(new_gene_path, 'insertion_location.tsv'))
			self.list_of_dict_filenames.append(os.path.join(new_gene_path, 'gene_sequences.tsv'))

			# These files need to be joined to existing files
			new_gene_shared_files = ['genes', 'rnas', 'proteins',
									 'rna_half_lives',
									 'protein_half_lives_measured']
			for f in new_gene_shared_files:
				file_path = os.path.join(new_gene_path, f + '.tsv')
				"""
				if these files dont exist, fill in with default values at a
				later point
				"""
				if os.path.isfile(os.path.join(FLAT_DIR, file_path)):
					self.list_of_dict_filenames.append(file_path)
					self.new_gene_added_data.update({f: nested_attr + f})

			rnaseq_path = os.path.join(new_gene_path, 'rnaseq_rsem_tpm_mean.tsv')
			if os.path.isfile(os.path.join(FLAT_DIR,rnaseq_path)):
				self.list_of_dict_filenames.append(rnaseq_path)
				self.new_gene_added_data.update({'rna_seq_data.rnaseq_rsem_tpm_mean':
													 nested_attr + 'rnaseq_rsem_tpm_mean'})

		# Load raw data from TSV files
		for filename in self.list_of_dict_filenames:
			self._load_tsv(FLAT_DIR, os.path.join(FLAT_DIR, filename))

		for filename in LIST_OF_PARAMETER_FILENAMES:
			self._load_parameters(os.path.join(FLAT_DIR, filename))

		self.genome_sequence = self._load_sequence(os.path.join(FLAT_DIR, SEQUENCE_FILE))

		self._prune_data()

		self._join_data()
		self._modify_data()

		if self.new_genes_option != 'off':
			self._check_new_gene_ids(nested_attr)

			insert_pos = self._update_gene_insertion_location(nested_attr)

			insertion_sequence = self._get_new_gene_sequence(nested_attr)

			insert_end = self._update_gene_locations(nested_attr, insert_pos)
			self.new_gene_added_data.update({'genes':
													  nested_attr+'genes'})

			self.genome_sequence = self.genome_sequence[:insert_pos] + \
								   insertion_sequence + \
								   self.genome_sequence[insert_pos:]
			assert self.genome_sequence[insert_pos:(insert_end + 1)] == \
				   insertion_sequence

			self.added_data = self.new_gene_added_data
			self._join_data()

	def _load_tsv(self, dir_name, file_name):
		path = self
		for sub_path in file_name[len(dir_name) + 1 : ].split(os.path.sep)[:-1]:
			if not hasattr(path, sub_path):
				setattr(path, sub_path, DataStore())
			path = getattr(path, sub_path)
		attr_name = file_name.split(os.path.sep)[-1].split(".")[0]
		setattr(path, attr_name, [])

		rows = read_tsv(file_name)
		setattr(path, attr_name, rows)

	def _load_sequence(self, file_path):
		from Bio import SeqIO

		with open(file_path, "r") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				return record.seq

	def _load_parameters(self, file_path):
		attr_name = file_path.split(os.path.sep)[-1].split(".")[0]
		param_dict = {}

		with io.open(file_path, "rb") as csvfile:
			reader = tsv.dict_reader(csvfile)

			for row in reader:
				value = json.loads(row['value'])
				if row['units'] != '':
					# `eval()` the units [risky!] then strip it to just a unit
					# since `a_list * a_float` (like `1.0 [1/s]`) fails, and
					# `a_list * an_int` repeats the list, which is also broken.
					unit = eval(row['units'])   # risky!
					unit = units.getUnit(unit)  # strip
					value = value * unit
				param_dict[row['name']] = value

		setattr(self, attr_name, param_dict)

	def _prune_data(self):
		"""
		Remove rows that are specified to be removed. Data will only be removed
		if all data in a row in the file specifying rows to be removed matches
		the same columns in the raw data file.
		"""

		# Check each pair of files to be removed
		for data_attr, attr_to_remove in self.removed_data.items():
			# Build the set of data to identify rows to be removed
			data_to_remove = getattr(self, attr_to_remove.split('.')[0])
			for attr in attr_to_remove.split('.')[1:]:
				data_to_remove = getattr(data_to_remove, attr)
			removed_cols = list(data_to_remove[0].keys())
			ids_to_remove = set()
			for row in data_to_remove:
				ids_to_remove.add(tuple([row[col] for col in removed_cols]))

			# Remove any matching rows
			data = getattr(self, data_attr)
			n_entries = len(data)
			removed_ids = set()
			for i, row in enumerate(data[::-1]):
				checked_id = tuple([row[col] for col in removed_cols])
				if checked_id in ids_to_remove:
					data.pop(n_entries - i - 1)
					removed_ids.add(checked_id)

			# Print warnings for entries that were marked to be removed that
			# does not exist in the original data file. Fold changes are
			# excluded since the original entries are split between two files.
			if not data_attr.startswith('fold_changes'):
				for unremoved_id in (ids_to_remove - removed_ids):
					print(f'Warning: Could not remove row {unremoved_id} '
						  f'in flat file {data_attr} because the row does not '
						  f'exist.')

	def _join_data(self):
		"""
		Add rows that are specified in additional files. Data will only be added
		if all the loaded columns from both datasets match.
		"""

		# Join data for each file with data to be added
		for data_attr, attr_to_add in self.added_data.items():
			# Get datasets to join
			data = getattr(self, data_attr.split('.')[0])
			for attr in data_attr.split('.')[1:]:
				data = getattr(data, attr)

			added_data = getattr(self, attr_to_add.split('.')[0])
			for attr in attr_to_add.split('.')[1:]:
				added_data = getattr(added_data, attr)

			# Check columns are the same for each dataset
			col_diff = set(data[0].keys()).symmetric_difference(added_data[0].keys())
			if col_diff:
				raise ValueError(f'Could not join datasets {data_attr} and {attr_to_add} '
					f'because columns do not match (different columns: {col_diff}).')

			# Join datasets
			for row in added_data:
				data.append(row)

	def _modify_data(self):
		"""
		Modify entires in rows that are specified to be modified. Rows must be
		identified by their entries in the first column (usually the ID column).
		"""
		# Check each pair of files to be modified
		for data_attr, modify_attr in self.modified_data.items():
			# Build the set of data to identify rows to be modified
			data_to_modify = getattr(self, modify_attr.split('.')[0])
			for attr in modify_attr.split('.')[1:]:
				data_to_modify = getattr(data_to_modify, attr)
			id_col_name = list(data_to_modify[0].keys())[0]

			id_to_modified_cols = {}
			for row in data_to_modify:
				id_to_modified_cols[row[id_col_name]] = row

			# Modify any matching rows with identical IDs
			data = getattr(self, data_attr)

			if list(data[0].keys())[0] != id_col_name:
				raise ValueError(f'Could not modify data {data_attr} with '
					f'{modify_attr} because the names of the first columns '
					f'do not match.')

			modified_entry_ids = set()
			for i, row in enumerate(data):
				if row[id_col_name] in id_to_modified_cols:
					modified_cols = id_to_modified_cols[row[id_col_name]]
					for col_name in data[i]:
						if col_name in modified_cols:
							data[i][col_name] = modified_cols[col_name]
					modified_entry_ids.add(row[id_col_name])

			# Check for entries in modification data that do not exist in
			# original data
			id_diff = set(id_to_modified_cols.keys()).symmetric_difference(
				modified_entry_ids)
			if id_diff:
				raise ValueError(f'Could not modify data {data_attr} with '
					f'{modify_attr} because of one or more entries in '
					f'{modify_attr} that do not exist in {data_attr} '
					f'(nonexistent entries: {id_diff}).')

	def _check_new_gene_ids(self, nested_attr):
		"""
		Check to ensure each new gene, RNA, and protein id starts with NG.
		"""
		nested_data = getattr(self, nested_attr[:-1].split('.')[0])
		for attr in nested_attr[:-1].split('.')[1:]:
			nested_data = getattr(nested_data, attr)

		new_genes_data = getattr(nested_data, 'genes')
		new_RNA_data = getattr(nested_data,'rnas')
		new_protein_data = getattr(nested_data,'proteins')

		for row in new_genes_data:
			assert row['id'].startswith("NG"), "ids of new genes must start " \
											   "with NG"
		for row in new_RNA_data:
			assert row['id'].startswith("NG"), "ids of new gene rnas must " \
											   "start with NG"
		for row in new_protein_data:
			assert row['id'].startswith("NG"), "ids of new gene proteins " \
											   "must start with NG"
		return


	def _update_gene_insertion_location(self, nested_attr):
		"""
		Update insertion location of new genes to prevent conflicts.
		"""

		genes_data = getattr(self, 'genes')
		tu_data = getattr(self, 'transcription_units')

		nested_data = getattr(self, nested_attr[:-1].split('.')[0])
		for attr in nested_attr[:-1].split('.')[1:]:
			nested_data = getattr(nested_data, attr)

		insert_loc_data = getattr(nested_data, 'insertion_location')

		assert len(insert_loc_data) == 1, 'each noncontiguous insertion should' \
										  ' be in its own directory'
		insert_pos = insert_loc_data[0]['insertion_pos']

		if not tu_data:
			# Check if specified insertion location is in another gene
			data_to_check = genes_data
		else:
			# Check if specified insertion location is in a transcription unit
			data_to_check = tu_data

		conflicts = [row for row in data_to_check if
					 ((row['left_end_pos'] is not None) and
					  (row['left_end_pos'] != '')) and
					 ((row['right_end_pos'] is not None) and
					 (row['left_end_pos'] != '')) and
					 (row['left_end_pos'] < insert_pos) and
					 (row['right_end_pos'] >= insert_pos)]
		# Change insertion location to after conflicts
		if conflicts:
			shift = max([ sub['right_end_pos'] for sub in conflicts ]) - \
					insert_pos + 1
			insert_pos = insert_pos + shift

		return insert_pos


	def _update_gene_locations(self, nested_attr, insert_pos):
		"""
		Modify positions of original genes based upon the insertion location
		of new genes. Returns end position of the gene insertion.
		"""

		genes_data = getattr(self, 'genes')
		tu_data = getattr(self, 'transcription_units')

		nested_data = getattr(self, nested_attr[:-1].split('.')[0])
		for attr in nested_attr[:-1].split('.')[1:]:
			nested_data = getattr(nested_data, attr)

		new_genes_data = getattr(nested_data,'genes')
		new_genes_data = sorted(new_genes_data, key=lambda d: d['left_end_pos'])

		for i in range(len(new_genes_data) - 1):
			assert new_genes_data[i+1]['left_end_pos'] == new_genes_data[i][
				'right_end_pos'] + 1, \
				"gaps in new gene insertions are not supported at this time"

		insert_end = new_genes_data[-1]['right_end_pos'] + insert_pos

		# Update global positions of original genes
		insert_len = insert_end - insert_pos + 1
		for row in genes_data:
			left = row['left_end_pos']
			right = row['right_end_pos']
			if (left is not None) and (right is not None) and left >= \
					insert_pos:
					row.update({'left_end_pos': left+insert_len})
					row.update({'right_end_pos': right+insert_len})

		# Update global positions of transcription units
		if tu_data:
			for row in tu_data:
				left = row['left_end_pos']
				right = row['right_end_pos']
				if (left != '') and (right != '') and left >= insert_pos:
					row.update({'left_end_pos': left + insert_len})
					row.update({'right_end_pos': right + insert_len})

		# Change relative insertion positions to global in reference genome
		for row in new_genes_data:
			left = row['left_end_pos']
			right = row['right_end_pos']
			row.update({'left_end_pos': left + insert_pos})
			row.update({'right_end_pos': right + insert_pos})

		return insert_end

	def _get_new_gene_sequence(self, nested_attr):
		"""
		Determine genome sequnce for insertion using the sequences and
		relative locations of the new genes.
		"""
		from Bio import Seq

		nested_data = getattr(self, nested_attr[:-1].split('.')[0])
		for attr in nested_attr[:-1].split('.')[1:]:
			nested_data = getattr(nested_data, attr)

		new_genes_data = getattr(nested_data,'genes')
		seq_data = getattr(nested_data,'gene_sequences')

		insertion_seq = Seq.Seq('')
		new_genes_data = sorted(new_genes_data, key=lambda d: d['left_end_pos'])
		assert new_genes_data[0]['left_end_pos'] == 0, \
			'first gene in new sequence must start at relative coordinate 0'

		for gene in new_genes_data:
			if gene['direction'] == "+":
				seq_row = next((row for row in seq_data
								if row['id'] == gene['id']), None)
				seq_string = seq_row['gene_seq']
				seq_addition = Seq.Seq(seq_string)
				insertion_seq += seq_addition
			else:
				seq_row = next((row for row in seq_data
								if row['id'] == gene['id']), None)
				seq_string = seq_row['gene_seq']
				seq_addition = Seq.Seq(seq_string)
				insertion_seq += seq_addition.reverse_complement()

			assert len(seq_addition) == (gene['right_end_pos'] -
				   gene['left_end_pos'] + 1),\
				"left and right end positions must agree with actual " \
				"sequence length for " + gene['id']

		return insertion_seq
