"""
Mass

Mass listener. Represents the total cellular mass.
"""

# TODO: generalize this logic for use with a generic simulation

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener
from wholecell.utils import units
import six
from wholecell.utils.migration.write_json import write_json

class Mass(wholecell.listeners.listener.Listener):
	""" Mass """

	_name = 'Mass'

	# Constructor
	def __init__(self, *args, **kwargs):
		# References to other internal states
		self.internal_states = None

		# NOTE: molecule weight is converted to femtograms/molecule from
		# grams/mol in BulkMolecules
		self.massUnits = 'fg'

		super(Mass, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Mass, self).initialize(sim, sim_data)

		self.internal_states = sim.internal_states

		self.processNames = list(sim.processes.keys())

		self.cellCycleLen = sim_data.condition_to_doubling_time[sim_data.condition].asNumber(units.s)

		self.rnaIndexes = np.array([
			sim_data.submass_name_to_index[name]
			for name in ["rRNA", "tRNA", "mRNA", "miscRNA", "nonspecific_RNA"]
			])

		self.rRnaIndex = sim_data.submass_name_to_index["rRNA"]
		self.smallMoleculeIndex = sim_data.submass_name_to_index["metabolite"]
		self.tRnaIndex = sim_data.submass_name_to_index["tRNA"]
		self.mRnaIndex = sim_data.submass_name_to_index["mRNA"]
		self.dnaIndex = sim_data.submass_name_to_index["DNA"]
		self.proteinIndex = sim_data.submass_name_to_index["protein"]
		self.waterIndex = sim_data.submass_name_to_index["water"]

		self.projection_index = sim_data.compartment_id_to_index["CCO-CELL-PROJECTION"]
		self.cytosol_index = sim_data.compartment_id_to_index["CCO-CYTOSOL"]
		self.extracellular_index = sim_data.compartment_id_to_index["CCO-EXTRACELLULAR"]
		self.flagellum_index = sim_data.compartment_id_to_index["CCO-FLAGELLUM"]
		self.membrane_index = sim_data.compartment_id_to_index["CCO-MEMBRANE"]
		self.outer_membrane_index = sim_data.compartment_id_to_index["CCO-OUTER-MEM"]
		self.periplasm_index = sim_data.compartment_id_to_index["CCO-PERI-BAC"]
		self.pilus_index = sim_data.compartment_id_to_index["CCO-PILUS"]
		self.inner_membrane_index = sim_data.compartment_id_to_index["CCO-PM-BAC-NEG"]


		self.cellDensity = sim_data.constants.cell_density.asNumber(units.g / units.L)

		# Set initial values
		self.setInitial = False

		self.dryMass = 0.0
		# TODO: set initial masses based on some calculations of the expected
		# mother cell (divided by two) in the last time step

		# Register logged quantities
		self.registerLoggedQuantity(
			"Cell mass\n(fg)",
			"cellMass",
			".2f"
			)

		self.registerLoggedQuantity(
			"Dry mass\n(fg)",
			"dryMass",
			".2f"
			)

		self.registerLoggedQuantity(
			"Dry mass\nfold change",
			"dryMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Expected\nfold change",
			"expectedMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Growth\n(fg/s)",
			"growth",
			".4f"
			)

		self.registerLoggedQuantity(
			"Protein\nfraction",
			"proteinMassFraction",
			".3f"
			)

		self.registerLoggedQuantity(
			"Protein\nfold change",
			"proteinMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"RNA\nfraction",
			"rnaMassFraction",
			".3f"
			)

		self.registerLoggedQuantity(
			"RNA\nfold change",
			"rnaMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Small mol\nfold change",
			"smallMoleculeFoldChange",
			".3f"
			)

		self.saved = False
		self.update_to_save = {}


	def update(self):
		oldDryMass = self.dryMass
		self.update_to_save['old_dry_mass'] = oldDryMass

		all_submasses = sum(
			state.mass() for state in six.viewvalues(self.internal_states))

		compartment_submasses = sum(
			state.compartment_mass() for state in six.viewvalues(
			self.internal_states))

		self.cellMass = all_submasses.sum()  # sum over all submasses
		self.update_to_save['cell_mass'] = self.cellMass

		self.waterMass = all_submasses[self.waterIndex]
		self.dryMass = self.cellMass - self.waterMass
		self.rnaMass = all_submasses[self.rnaIndexes].sum()
		self.rRnaMass = all_submasses[self.rRnaIndex]
		self.tRnaMass = all_submasses[self.tRnaIndex]
		self.mRnaMass = all_submasses[self.mRnaIndex]
		self.dnaMass = all_submasses[self.dnaIndex]
		self.proteinMass = all_submasses[self.proteinIndex]
		self.smallMoleculeMass = all_submasses[self.smallMoleculeIndex]

		self.update_to_save['water_mass'] = self.waterMass
		self.update_to_save['dry_mass'] = self.dryMass
		self.update_to_save['rnaMass'] = self.rnaMass
		self.update_to_save['rRnaMass'] = self.rRnaMass
		self.update_to_save['tRnaMass'] = self.tRnaMass
		self.update_to_save['mRnaMass'] = self.mRnaMass
		self.update_to_save['dnaMass'] = self.dnaMass
		self.update_to_save['proteinMass'] = self.proteinMass
		self.update_to_save['smallMoleculeMass'] = self.smallMoleculeMass

		self.projection_mass = compartment_submasses[self.projection_index, :].sum()
		self.cytosol_mass = compartment_submasses[self.cytosol_index, :].sum()
		self.extracellular_mass = compartment_submasses[self.extracellular_index, :].sum()
		self.flagellum_mass = compartment_submasses[self.flagellum_index, :].sum()
		self.membrane_mass = compartment_submasses[self.membrane_index, :].sum()
		self.outer_membrane_mass = compartment_submasses[self.outer_membrane_index, :].sum()
		self.periplasm_mass = compartment_submasses[self.periplasm_index, :].sum()
		self.pilus_mass = compartment_submasses[self.pilus_index, :].sum()
		self.inner_membrane_mass = compartment_submasses[self.inner_membrane_index, :].sum()

		self.update_to_save['projection_mass'] = self.projection_mass
		self.update_to_save['cytosol_mass'] = self.cytosol_mass
		self.update_to_save['extracellular_mass'] = self.extracellular_mass
		self.update_to_save['flagellum_mass'] = self.flagellum_mass
		self.update_to_save['membrane_mass'] = self.membrane_mass
		self.update_to_save['outer_membrane_mass'] = self.outer_membrane_mass
		self.update_to_save['periplasm_mass'] = self.periplasm_mass
		self.update_to_save['pilus_mass'] = self.pilus_mass
		self.update_to_save['inner_membrane_mass'] = self.inner_membrane_mass

		# TODO (Eran) use this volume everywhere in the codebase that is currently calculating volume
		self.volume = self.cellMass / self.cellDensity
		self.update_to_save['volume'] = self.volume

		self.processMassDifferences = sum(
			state.process_mass_diffs() for state in six.viewvalues(self.internal_states)
			).sum(axis=1)

		self.update_to_save['process_mass_diffs'] = self.processMassDifferences

		if self.simulationStep() > 0:
			self.growth = self.dryMass - oldDryMass

		else:
			self.growth = np.nan

		self.update_to_save['growth'] = self.growth

		self.instantaniousGrowthRate = self.growth / self.timeStepSec() / self.dryMass
		self.update_to_save['instantanious_growth_rate'] = self.instantaniousGrowthRate

		self.proteinMassFraction = self.proteinMass / self.dryMass
		self.rnaMassFraction = self.rnaMass / self.dryMass

		self.update_to_save['protein_mass_fraction'] = self.proteinMassFraction
		self.update_to_save['rna_mass_fraction'] = self.rnaMassFraction

		if not self.setInitial:
			self.setInitial = True

			self.timeInitial = self.time()

			self.dryMassInitial = self.dryMass
			self.proteinMassInitial = self.proteinMass
			self.rnaMassInitial = self.rnaMass
			self.smallMoleculeMassInitial = self.smallMoleculeMass

			self.update_to_save['dry_mass_initial'] = self.dryMassInitial
			self.update_to_save['protein_mass_initial'] = self.proteinMassInitial
			self.update_to_save['rna_mass_initial'] = self.rnaMassInitial
			self.update_to_save['small_molecule_mass_initial'] = self.smallMoleculeMassInitial


		self.dryMassFoldChange = self.dryMass / self.dryMassInitial
		self.proteinMassFoldChange = self.proteinMass / self.proteinMassInitial
		self.rnaMassFoldChange = self.rnaMass / self.rnaMassInitial
		self.smallMoleculeFoldChange = self.smallMoleculeMass / self.smallMoleculeMassInitial

		self.update_to_save['dry_mass_fold_change'] = self.dryMassFoldChange
		self.update_to_save['protein_mass_fold_change'] = self.proteinMassFoldChange
		self.update_to_save['rna_mass_fold_change'] = self.rnaMassFoldChange
		self.update_to_save['small_molecule_fold_change'] = self.smallMoleculeFoldChange

		self.expectedMassFoldChange = np.exp(np.log(2) * (self.time() - self.timeInitial) / self.cellCycleLen)
		self.update_to_save['expected_mass_fold_change'] = self.expectedMassFoldChange

		if not self.saved:
			write_json(f'out/migration/mass_listener_update_t{int(self._sim.time())}.json',
					   self.update_to_save)
					   
			self.saved = True

			import ipdb; ipdb.set_trace()


	def tableCreate(self, tableWriter):
		# Store units as metadata
		tableWriter.writeAttributes(
			cell_units = self.massUnits,
			cellDry_units = self.massUnits,
			metabolite_units = self.massUnits,
			growth_units = self.massUnits,
			rna_units = self.massUnits,
			protein_units = self.massUnits,
			water_units = self.massUnits,
			nucleoid_units = self.massUnits,
			processNames = self.processNames,
			smallMoleculeMass = self.smallMoleculeMass,
			subcolumns = {})


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			cellMass = self.cellMass,
			growth = self.growth,
			dryMass = self.dryMass,
			rnaMass = self.rnaMass,
			rRnaMass = self.rRnaMass,
			tRnaMass = self.tRnaMass,
			mRnaMass = self.mRnaMass,
			dnaMass = self.dnaMass,
			proteinMass = self.proteinMass,
			waterMass = self.waterMass,
			projection_mass = self.projection_mass,
			cytosol_mass = self.cytosol_mass,
			extracellular_mass = self.extracellular_mass,
			flagellum = self.flagellum_mass,
			membrane_mass = self.membrane_mass,
			outer_membrane_mass = self.outer_membrane_mass,
			periplasm_mass = self.periplasm_mass,
			pilus_mass = self.pilus_mass,
			inner_membrane_mass = self.inner_membrane_mass,
			processMassDifferences = self.processMassDifferences.astype(np.float64),
			smallMoleculeMass = self.smallMoleculeMass,
			instantaniousGrowthRate = self.instantaniousGrowthRate,
			cellVolume = self.volume
			)
