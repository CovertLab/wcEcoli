"""
SimulationData for transcription process

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015

TODO: add mapping of tRNA to charged tRNA if allowing more than one modified form of tRNA and separate mappings for tRNA and charged tRNA to AA
TODO: handle ppGpp and DksA-ppGpp regulation separately
"""

from __future__ import division

import numpy as np
from scipy import interpolate

from wholecell.utils import units
from wholecell.utils.fitting import normalize
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize
from wholecell.utils.random import make_elongation_rates


RNA_SEQ_ANALYSIS = "rsem_tpm"
KCAT_ENDO_RNASE = 0.001
ESTIMATE_ENDO_RNASES = 5000
MAX_TIMESTEP_LEN = 2  # Determines length of padding values to add to transcript sequence matrix
PPGPP_CONC_UNITS = units.umol / units.L


class Transcription(object):
	"""
	SimulationData for the transcription process
	"""

	def __init__(self, raw_data, sim_data):
		self._build_ppgpp_regulation(raw_data, sim_data)
		self._build_rna_data(raw_data, sim_data)
		self._build_transcription(raw_data, sim_data)
		self._build_charged_trna(raw_data, sim_data)
		self._build_elongation_rates(raw_data, sim_data)

	def _build_ppgpp_regulation(self, raw_data, sim_data):
		"""
		Determine which genes are regulated by ppGpp and store the fold
		change in expression associated with each RNA.
		"""

		def read_value(d, k):
			"""Handle empty values from raw_data as 0"""
			val = d[k]
			return 0 if val == '' else val

		# Read regulation data from raw_data
		# Treats ppGpp and DksA-ppGpp regulation the same
		gene_to_rna = {g['symbol']: g['rnaId'] for g in raw_data.genes}
		regulation = {}
		for reg in raw_data.ppgpp_regulation:
			# Convert to regulated RNA
			gene = reg['Gene']
			rna = gene_to_rna.get(gene, None)
			if rna is None:
				continue

			# Add additional gene symbols for matching FC data
			curated_gene = reg['Curated Gene']
			if gene != curated_gene:
				gene_to_rna[curated_gene] = rna

			# Update value (some genes are repeated in raw_data))
			direction = read_value(reg, 'ppGpp') + read_value(reg, 'DksA-ppGpp')
			regulation[rna] = regulation.get(rna, 0) + direction

		# Read fold change data from raw_data
		valid_categories = {'A', 'B', 'C', 'D'}  # FC data is significant
		sample_time = 5  # Could also be 10
		sample_id = '1+2+ {} min'.format(sample_time)
		rna_fold_changes = {}
		for fc in raw_data.ppgpp_fc:
			# Convert to regulated RNA
			gene = fc['Gene']
			rna = gene_to_rna.get(gene, None)
			if rna is None:
				continue

			category = fc['{} Category'.format(sample_id)]
			if category not in valid_categories:
				continue
			rna_fold_changes[rna] = fc[sample_id]

		# Store arrays of regulation
		regulated_genes = []
		regulation_direction = []
		fold_changes = []
		for rna in sorted(regulation):
			reg_dir = regulation[rna]
			fc_dir = rna_fold_changes.get(rna, 0)

			# Ignore inconsistent regulatory directions
			if reg_dir == 0 or reg_dir * fc_dir < 0:
				continue

			regulated_genes.append(rna)
			regulation_direction.append(np.sign(reg_dir))
			fold_changes.append(fc_dir)
		self.ppgpp_regulated_genes = np.array(regulated_genes)
		regulation_direction = np.array(regulation_direction)

		# Replace fold changes without data with the average
		fold_changes = np.array(fold_changes)
		average_negative_fc = fold_changes[fold_changes < 0].mean()
		average_positive_fc = fold_changes[fold_changes > 0].mean()
		fold_changes[(fold_changes == 0) & (regulation_direction < 0)] = -0.9316  # from solve_km.py
		fold_changes[(fold_changes == 0) & (regulation_direction > 0)] = average_positive_fc
		self.ppgpp_fold_changes = fold_changes

		self.ppgpp_km = 0.024  # TODO: use solve_km script here to explicitly solve or use max fold change

		# Predict growth rate from ppGpp level
		per_dry_mass_to_per_volume = sim_data.constants.cellDensity * (1. - raw_data.mass_parameters['cellWaterMassFraction'])
		growth_rates = np.log(2) / np.array([d['doublingTime'].asNumber(units.s) for d in raw_data.growthRateDependentParameters])
		ppgpp = np.array([d['ppGpp_conc'].asNumber(units.umol / units.g) for d in raw_data.growthRateDependentParameters]) * per_dry_mass_to_per_volume.asNumber(units.g / units.L)
		self.ppgpp_growth_parameters = interpolate.splrep(ppgpp[::-1], growth_rates[::-1], k=1)

	def _build_rna_data(self, raw_data, sim_data):
		"""
		Build RNA-associated simulation data from raw data.
		"""
		assert all([len(rna['location']) == 1 for rna in raw_data.rnas])

		# Loads RNA IDs, degradation rates, lengths, and nucleotide compositions
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location'][0])
            for rna in raw_data.rnas]
		rnaDegRates = np.log(2) / np.array([rna['halfLife'] for rna in raw_data.rnas]) # TODO: units
		rnaLens = np.array([len(rna['seq']) for rna in raw_data.rnas])

		ntCounts = np.array([
			(rna['seq'].count('A'), rna['seq'].count('C'),
			rna['seq'].count('G'), rna['seq'].count('U'))
			for rna in raw_data.rnas
			])

		# Load RNA expression from RNA-seq data
		expression = []

		for rna in raw_data.rnas:
			arb_exp = [x[sim_data.basal_expression_condition]
                for x in eval("raw_data.rna_seq_data.rnaseq_{}_mean".format(RNA_SEQ_ANALYSIS))
                if x['Gene'] == rna['geneId']]

			# If sequencing data is not found for rRNA or tRNA, initialize
            # expression to zero. For other RNA types, raise exception.
			if len(arb_exp) > 0:
				expression.append(arb_exp[0])
			elif rna['type'] == 'mRNA' or rna['type'] == 'miscRNA':
				raise Exception('No RNA-seq data found for {}'.format(rna['id']))
			elif rna['type'] == 'rRNA' or rna['type'] == 'tRNA':
				expression.append(0.)
			else:
				raise Exception('Unknown RNA {}'.format(rna['id']))

		expression = np.array(expression)

		# Calculate synthesis probabilities from expression and normalize
		synthProb = expression*(
			np.log(2) / sim_data.doubling_time.asNumber(units.s)
			+ rnaDegRates
			)
		synthProb /= synthProb.sum()

		# Calculate EndoRNase Km values
		Km = (KCAT_ENDO_RNASE*ESTIMATE_ENDO_RNASES/rnaDegRates) - expression

		# Load molecular weights and gene IDs
		mws = np.array([rna['mw'] for rna in raw_data.rnas]).sum(axis = 1)
		geneIds = np.array([rna['geneId'] for rna in raw_data.rnas])

		# Construct boolean arrays and index arrays for each rRNA type
		n_rnas = len(rnaIds)
		is_23S = np.zeros(n_rnas, dtype = np.bool)
		is_16S = np.zeros(n_rnas, dtype = np.bool)
		is_5S = np.zeros(n_rnas, dtype = np.bool)
		idx_23S = []
		idx_16S = []
		idx_5S = []

		for rnaIndex, rna in enumerate(raw_data.rnas):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is_23S[rnaIndex] = True
				idx_23S.append(rnaIndex)

			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is_16S[rnaIndex] = True
				idx_16S.append(rnaIndex)

			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is_5S[rnaIndex] = True
				idx_5S.append(rnaIndex)

		idx_23S = np.array(idx_23S)
		idx_16S = np.array(idx_16S)
		idx_5S = np.array(idx_5S)

		# Load sequence data
		sequences = [rna['seq'] for rna in raw_data.rnas]
		maxSequenceLength = max(len(sequence) for sequence in sequences)

		# Load IDs of protein monomers
		monomerIds = [rna['monomerId'] for rna in raw_data.rnas]

		# Get index of gene corresponding to each RNA
		gene_index = {gene["rnaId"]: i
			for i, gene in enumerate(raw_data.genes)}

		# Get list of coordinates and directions for each gene
		coordinate_list = [gene["coordinate"] for gene in raw_data.genes]
		direction_list = [gene["direction"] for gene in raw_data.genes]

		oric_coordinate = raw_data.parameters['oriCCenter'].asNumber()
		terc_coordinate = raw_data.parameters['terCCenter'].asNumber()
		genome_length = len(raw_data.genome_sequence)

		def get_relative_coordinates(coordinates):
			relative_coordinates = ((coordinates - terc_coordinate)
				% genome_length + terc_coordinate - oric_coordinate
				)

			if relative_coordinates < 0:
				relative_coordinates += 1

			return relative_coordinates

		# Location of transcription initiation relative to origin
		replicationCoordinate = [
			get_relative_coordinates(coordinate_list[gene_index[rna["id"]]])
			for rna in raw_data.rnas]

		# Direction of transcription
		direction = [
			(direction_list[gene_index[rna["id"]]] == "+")
			for rna in raw_data.rnas]

		# Set the lengths, nucleotide counts, molecular weights, and sequences
		# of each type of rRNAs to be identical to those of the first rRNA
		# operon. Later in the sim, transcription of all rRNA genes are set to
		# produce the rRNAs of the first operon. This is done to simplify the
		# complexation reactions that form ribosomes. In reality, all of these
		# genes produce rRNA molecules with slightly different sequences and
		# molecular weights.
		rnaLens[idx_23S] = rnaLens[idx_23S[0]]
		rnaLens[idx_16S] = rnaLens[idx_16S[0]]
		rnaLens[idx_5S] = rnaLens[idx_5S[0]]

		ntCounts[idx_23S, :] = ntCounts[idx_23S[0], :]
		ntCounts[idx_16S, :] = ntCounts[idx_16S[0], :]
		ntCounts[idx_5S, :] = ntCounts[idx_5S[0], :]

		mws[idx_23S] = mws[idx_23S[0]]
		mws[idx_16S] = mws[idx_16S[0]]
		mws[idx_5S] = mws[idx_5S[0]]

		for idx in idx_23S[1:]:
			sequences[idx] = sequences[idx_23S[0]]

		for idx in idx_16S[1:]:
			sequences[idx] = sequences[idx_16S[0]]

		for idx in idx_5S[1:]:
			sequences[idx] = sequences[idx_5S[0]]

		rnaData = np.zeros(
			n_rnas,
			dtype = [
				('id', 'a50'),
				('degRate', 'f8'),
				('length', 'i8'),
				('countsACGU', '4i8'),
				('mw', 'f8'),
				('isMRna', 'bool'),
				('isMiscRna', 'bool'),
				('isRRna', 'bool'),
				('isTRna', 'bool'),
				('isRRna23S', 'bool'),
				('isRRna16S', 'bool'),
				('isRRna5S', 'bool'),
				('isRProtein', 'bool'),
				('isRnap',	'bool'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				('geneId', 'a50'),
				('KmEndoRNase', 'f8'),
				('replicationCoordinate', 'int64'),
				('direction', 'bool'),
				]
			)

		rnaData['id'] = rnaIds
		rnaData['degRate'] = rnaDegRates
		rnaData['length'] = rnaLens
		rnaData['countsACGU'] = ntCounts
		rnaData['mw'] = mws
		rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in raw_data.rnas]
		rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in raw_data.rnas]
		rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in raw_data.rnas]
		rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in raw_data.rnas]
		rnaData['isRProtein'] = [
            "{}[c]".format(x) in sim_data.moleculeGroups.rProteins
            for x in monomerIds]
		rnaData['isRnap'] = [
            "{}[c]".format(x) in sim_data.moleculeGroups.rnapIds
            for x in monomerIds]
		rnaData['isRRna23S'] = is_23S
		rnaData['isRRna16S'] = is_16S
		rnaData['isRRna5S'] = is_5S
		rnaData['sequence'] = sequences
		rnaData['geneId'] = geneIds
		rnaData['KmEndoRNase'] = Km
		rnaData['replicationCoordinate'] = replicationCoordinate
		rnaData['direction'] = direction

		field_units = {
			'id': None,
			'degRate': 1 / units.s,
			'length': units.nt,
			'countsACGU': units.nt,
			'mw': units.g / units.mol,
			'isMRna': None,
			'isMiscRna': None,
			'isRRna': None,
			'isTRna': None,
			'isRRna23S': None,
			'isRRna16S': None,
			'isRRna5S':	None,
			'isRProtein': None,
			'isRnap': None,
			'sequence': None,
			'geneId': None,
			'KmEndoRNase': units.mol / units.L,
			'replicationCoordinate': None,
			'direction': None,
			}

		self.rnaExpression = {}
		self.rnaSynthProb = {}

		# Set basal expression and synthesis probabilities - conditional values
        # are set in the parca.
		self.rnaExpression["basal"] = expression / expression.sum()
		self.rnaSynthProb["basal"] = synthProb / synthProb.sum()

		self.rnaData = UnitStructArray(rnaData, field_units)


	def _build_transcription(self, raw_data, sim_data):
		"""
		Build transcription-associated simulation data from raw data.
		"""
		sequences = self.rnaData["sequence"] # TODO: consider removing sequences

		# Construct transcription sequence matrix
		maxLen = np.int64(
			self.rnaData["length"].asNumber().max()
			+ MAX_TIMESTEP_LEN * sim_data.growthRateParameters.rnaPolymeraseElongationRate.asNumber(units.nt/units.s)
			)

		self.transcriptionSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.transcriptionSequences.fill(polymerize.PAD_VALUE)

		ntMapping = {ntpId: i for i, ntpId in enumerate(["A", "C", "G", "U"])}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.transcriptionSequences[i, j] = ntMapping[letter]

		# Calculate weights of transcript nucleotide monomers
		self.transcriptionMonomerWeights = (
			(
				sim_data.getter.getMass(sim_data.moleculeGroups.ntpIds)
				- sim_data.getter.getMass(["PPI[c]"])
				)
			/ sim_data.constants.nAvogadro
			).asNumber(units.fg)

		self.transcriptionEndWeight = ((sim_data.getter.getMass(["PPI[c]"])
            / sim_data.constants.nAvogadro).asNumber(units.fg))

	def _build_charged_trna(self, raw_data, sim_data):
		'''
		Loads information and creates data structures necessary for charging of tRNA

		Note:
			Requires self.rnaData so can't be built in translation even if some
			data structures would be more appropriate there.
		'''

		# Create list of charged tRNAs
		trna_names = self.rnaData['id'][self.rnaData['isTRna']]
		charged_trnas = [x['modifiedForms'] for x in raw_data.rnas if x['id'] + '[c]' in trna_names]
		filtered_charged_trna = []
		for charged_list in charged_trnas:
			for trna in charged_list:
				# Skip modified forms so only one charged tRNA per uncharged tRNA
				if 'FMET' in trna or 'modified' in trna:
					continue

				assert('c' in sim_data.getter.getLocation([trna])[0])
				filtered_charged_trna += [trna + '[c]']

		self.charged_trna_names = filtered_charged_trna
		assert(len(self.charged_trna_names) == len(trna_names))

		# Create mapping of each tRNA/charged tRNA to associated AA
		trna_dict = {
			'RNA0-300[c]': 'VAL',
			'RNA0-301[c]': 'LYS',
			'RNA0-302[c]': 'LYS',
			'RNA0-303[c]': 'LYS',
			'RNA0-304[c]': 'ASN',
			'RNA0-305[c]': 'ILE',
			'RNA0-306[c]': 'MET',
			}
		aa_names = sim_data.moleculeGroups.aaIDs
		aa_indices = {aa: i for i, aa in enumerate(aa_names)}
		trna_indices = {trna: i for i, trna in enumerate(trna_names)}
		self.aa_from_trna = np.zeros((len(aa_names), len(trna_names)))
		for trna in trna_names:
			aa = trna[:3].upper()
			if aa == 'ALA':
				aa = 'L-ALPHA-ALANINE'
			elif aa == 'ASP':
				aa = 'L-ASPARTATE'
			elif aa == 'SEL':
				aa = 'L-SELENOCYSTEINE'
			elif aa == 'RNA':
				aa = trna_dict[trna]

			assert('c' in sim_data.getter.getLocation([aa])[0])
			aa += '[c]'
			if aa in aa_names:
				aa_idx = aa_indices[aa]
				trna_idx = trna_indices[trna]
				self.aa_from_trna[aa_idx, trna_idx] = 1

		# Arrays for stoichiometry and synthetase mapping matrices
		molecules = []

		# Sparse matrix representation - i, j are row/column indices and v is value
		stoich_matrix_i = []
		stoich_matrix_j = []
		stoich_matrix_v = []

		synthetase_names = []
		synthetase_mapping_aa = []
		synthetase_mapping_syn = []

		# Create stoichiometry matrix for charging reactions
		for reaction in raw_data.modificationReactions:
			# Skip reactions from modificationReactions that don't have both an uncharged and charged tRNA
			no_charged_trna_in_reaction = True
			no_trna_in_reaction = True
			for mol in [molecule['molecule'] + '[' + molecule['location'] + ']' for molecule in reaction['stoichiometry']]:
				if mol in self.charged_trna_names:
					no_charged_trna_in_reaction = False

				if mol in trna_names:
					no_trna_in_reaction = False

			if no_charged_trna_in_reaction or no_trna_in_reaction:
				continue

			assert reaction['process'] == 'rna'
			assert reaction['dir'] == 1

			# Get uncharged tRNA name for the given reaction
			trna = None
			for mol in [molecule['molecule'] + '[' + molecule['location'] + ']' for molecule in reaction['stoichiometry']]:
				if mol in trna_names:
					trna = mol
					break

			if trna is None:
				continue
			trna_index = trna_indices[trna]

			# Get molecule information
			aa_idx = None
			for molecule in reaction['stoichiometry']:
				molecule_prefix = molecule['molecule']
				if molecule['type'] == 'metabolite':
					molecule_prefix = molecule_prefix.upper()

				molecule_name = '{}[{}]'.format(
					molecule_prefix,
					molecule['location']
					)

				if molecule_name not in molecules:
					molecules.append(molecule_name)
					molecule_index = len(molecules) - 1
				else:
					molecule_index = molecules.index(molecule_name)

				aa_idx = aa_indices.get(molecule_name, aa_idx)

				coefficient = molecule['coeff']

				assert coefficient % 1 == 0

				stoich_matrix_i.append(molecule_index)
				stoich_matrix_j.append(trna_index)
				stoich_matrix_v.append(coefficient)

			assert aa_idx is not None

			# Create mapping for synthetases catalyzing charging
			for synthetase in reaction['catBy']:
				synthetase = '{}[{}]'.format(synthetase, molecule['location'])

				if synthetase not in synthetase_names:
					synthetase_names.append(synthetase)

				synthetase_mapping_aa.append(aa_idx)
				synthetase_mapping_syn.append(synthetase_names.index(synthetase))

		# Save matrices and related lists of names
		self._stoich_matrix_i = np.array(stoich_matrix_i)
		self._stoich_matrix_j = np.array(stoich_matrix_j)
		self._stoich_matrix_v = np.array(stoich_matrix_v)

		self.aa_from_synthetase = np.zeros((len(aa_names), len(synthetase_names)))
		self.aa_from_synthetase[synthetase_mapping_aa, synthetase_mapping_syn] = 1

		self.synthetase_names = synthetase_names
		self.charging_molecules = molecules

	def charging_stoich_matrix(self):
		'''
		Creates stoich matrix from i, j, v arrays

		Returns 2D array with rows of metabolites for each tRNA charging reaction on the column
		'''

		shape = (self._stoich_matrix_i.max() + 1, self._stoich_matrix_j.max() + 1)

		out = np.zeros(shape, np.float64)
		out[self._stoich_matrix_i, self._stoich_matrix_j] = self._stoich_matrix_v

		return out
	def _build_elongation_rates(self, raw_data, sim_data):
		self.max_elongation_rate = sim_data.constants.dnaPolymeraseElongationRateMax
		self.RRNA_indexes = np.where(self.rnaData['isRRna'])[0]

	def make_elongation_rates(self, random, base, time_step, variable_elongation=False):
		return make_elongation_rates(
			random,
			self.transcriptionSequences.shape[0],
			base,
			self.RRNA_indexes,
			self.max_elongation_rate,
			time_step,
			variable_elongation)

	def set_ppgpp_expression(self, sim_data):
		"""
		Called during the parca to determine expression of each gene for ppGpp
		bound and free RNAP.

		Attributes set:
			exp_ppgpp (ndarray[float]): expression for each gene when RNAP
				is bound to ppGpp
			exp_free (ndarray[float]): expression for each gene when RNAP
				is not bound to ppGpp
			ppgpp_km (float with units of mol / vol): KM for ppGpp binding to RNAP
			ppgpp_km_squared (float): squared and unitless version of KM for
				faster computation in other functions
		"""

		# TODO: calculate based on growth data
		km = 0.079**2
		ppgpp_aa = 0.058**2
		ppgpp_basal = 0.13**2
		# TODO: create function to determine
		f_ppgpp_aa = ppgpp_aa / (km + ppgpp_aa)
		f_ppgpp_basal = ppgpp_basal / (km + ppgpp_basal)

		rna_idx = {r[:-3]: i for i, r in enumerate(self.rnaData['id'])}
		fcs = np.zeros(len(self.rnaData))
		for rna, fc in zip(self.ppgpp_regulated_genes, self.ppgpp_fold_changes):
			fcs[rna_idx[rna]] = fc
		exp = self.rnaExpression['basal']
		self.exp_ppgpp = ((2**fcs * exp * (1 - f_ppgpp_aa) / (1 - f_ppgpp_basal))
			/ (1 - 2**fcs * (f_ppgpp_aa - f_ppgpp_basal * (1 - f_ppgpp_aa) / (1 - f_ppgpp_basal))))
		self.exp_free = (exp - self.exp_ppgpp*f_ppgpp_basal) / (1 - f_ppgpp_basal)
		self.exp_free[self.exp_free < 0] = 0  # fold change is limited by KM, can't have very high positive fold changes
		self.ppgpp_km = (units.umol / units.g * np.sqrt(km)
			* sim_data.constants.cellDensity * sim_data.mass.cellDryMassFraction)  # umol / L
		self._ppgpp_km_squared = self.ppgpp_km.asNumber(PPGPP_CONC_UNITS)**2  # save computation

	def adjust_polymerizing_ppgpp_expression(self):
		"""
		Adjust ppGpp expression based on fit for ribosome and RNAP physiological constraints
		using least squares fit for 3 conditions with different growth rates/ppGpp.

		Modifies attributes:
			exp_ppgpp (ndarray[float]): expression for each gene when RNAP
				is bound to ppGpp, adjusted for necessary RNAP and ribosome
				expression, normalized to 1
			exp_free (ndarray[float]): expression for each gene when RNAP
				is not bound to ppGpp, adjusted for necessary RNAP and ribosome
				expression, normalized to 1
		"""

		rna_data = self.rnaData
		exp_free = self.exp_free
		exp_ppgpp = self.exp_ppgpp

		# TODO: don't have hardcoded
		km = 0.079**2
		ppgpp_aa = 0.058**2
		ppgpp_basal = 0.13**2
		ppgpp_an = 0.316**2
		f_ppgpp_aa = ppgpp_aa / (km + ppgpp_aa)
		f_ppgpp_basal = ppgpp_basal / (km + ppgpp_basal)
		f_ppgpp_an = ppgpp_an / (km + ppgpp_an)

		adjusted_mask = rna_data['isRnap'] | rna_data['isRProtein'] | rna_data['isRRna']
		F = np.array([[1- f_ppgpp_aa, f_ppgpp_aa], [1 - f_ppgpp_basal, f_ppgpp_basal], [1 - f_ppgpp_an, f_ppgpp_an]])
		Flst = np.linalg.inv(F.T.dot(F)).dot(F.T)
		exp = np.array([
			self.rnaExpression['with_aa'],
			self.rnaExpression['basal'],
			self.rnaExpression['no_oxygen']])
		adjusted_free, adjusted_ppgpp = Flst.dot(exp)
		exp_free[adjusted_mask] = adjusted_free[adjusted_mask]
		exp_ppgpp[adjusted_mask] = adjusted_ppgpp[adjusted_mask]

		ppgpp_regulated = np.array([g[:-3] in self.ppgpp_regulated_genes for g in rna_data['id']])
		scale_free_by = (1 - exp_free[ppgpp_regulated].sum()) / exp_free[~ppgpp_regulated].sum()
		exp_free[~ppgpp_regulated] *= scale_free_by
		assert(scale_free_by > 0)
		scale_ppgpp_by = (1 - exp_ppgpp[ppgpp_regulated].sum()) / exp_ppgpp[~ppgpp_regulated].sum()
		exp_ppgpp[~ppgpp_regulated] *= scale_ppgpp_by
		assert(scale_ppgpp_by > 0)

	def fraction_rnap_bound_ppgpp(self, ppgpp):
		"""
		Calculates the fraction of RNAP expected to be bound to ppGpp
		at a given concentration of ppGpp.

		Args:
			ppgpp (float with or without mol / volume units): concentration of ppGpp,
				if unitless, should represent the concentration of PPGPP_CONC_UNITS

		Returns:
			float: fraction of RNAP that will be bound to ppGpp
		"""

		if units.hasUnit(ppgpp):
			ppgpp = ppgpp.asNumber(PPGPP_CONC_UNITS)

		return ppgpp**2 / (self._ppgpp_km_squared + ppgpp**2)



	def expression_from_ppgpp(self, ppgpp):
		"""
		Calculates the expression of each gene at a given concentration of ppGpp.

		Args:
			ppgpp (float with or without mol / volume units): concentration of ppGpp,
				if unitless, should represent the concentration of PPGPP_CONC_UNITS

		Returns:
			ndarray[float]: normalized expression for each gene
		"""

		f_ppgpp = self.fraction_rnap_bound_ppgpp(ppgpp)
		return normalize(self.exp_free * (1 - f_ppgpp) + self.exp_ppgpp * f_ppgpp)


	def synth_prob_from_ppgpp(self, ppgpp, copy_number):
		"""
		Calculates the synthesis probability of each gene at a given concentration
		of ppGpp.

		Args:
			ppgpp (float with mol / volume units): concentration of ppGpp
			copy_number (Callable[float, int]): function that gives the expected copy
				number given a doubling time and gene replication coordinate

		Returns
			ndarray[float]: normalized synthesis probability for each gene

		Note:
			copy_number should be sim_data.process.replication.get_average_copy_number
			but saving the functino handle as a clasa attribute prevents pickling of sim_data
			without additional handling
		"""

		ppgpp = ppgpp.asNumber(PPGPP_CONC_UNITS)
		f_ppgpp = self.fraction_rnap_bound_ppgpp(ppgpp)

		growth = max(interpolate.splev(ppgpp, self.ppgpp_growth_parameters), 0)
		tau = np.log(2) / growth / 60
		loss = growth + self.rnaData['degRate'].asNumber(1 / units.s)

		n_avg_copy = copy_number(tau, self.rnaData['replicationCoordinate'])

		return normalize((self.exp_free * (1 - f_ppgpp) + self.exp_ppgpp * f_ppgpp) * loss / n_avg_copy)
