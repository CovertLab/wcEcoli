import os
import cPickle
import pandas as pd
import numpy as np
from wholecell.io.tablereader import TableReader

class Reader:
	_sim_data_path = None
	_sim_out_path = None
	_sim_data = None
	_loaded = None

	def __init__(self, sim_data, sim_out, load=None):
		self._sim_data_path = sim_data
		self._sim_out_path = sim_out
		self._load_sim_data()

		self.load_main()
		self._load_data(load)

	def _load_sim_data(self):
		with open(self._sim_data_path, "rb") as f:
		    self._sim_data = cPickle.load(f)

	def _load_data(self, load=None):
		toLoad = self._loaders.values()

		if load is not None:
			toLoad = [v for k, v in self._loaders.iteritems() if k in load]

		for l in toLoad:
			l(self)

		self._loaded = toLoad

	def reload(self, load=None):
		self.load_main()

		if load is not None:
			self._load_data(load)
		else:
			if self._loaded is None:
				self._load_data()
			else:
				for l in self._loaded:
					l(self)

	def load_main(self):
		tr = TableReader(os.path.join(self._sim_out_path, "Main"))
		time = np.cumsum(tr.readColumn('timeStepSec'))
		time = np.insert(time[:-1], 0, 0.)
		time = pd.to_datetime(time, unit="s")

		self.time = time

	def load_bulk_molecules(self):
		tr = TableReader(os.path.join(self._sim_out_path, "BulkMolecules"))

		bm = ReaderData()
		bm.counts = pd.DataFrame(tr.readColumn('counts'), index=self.time, columns=tr.readAttribute('objectNames'))
		self.bulkMolecules = bm

	def load_enzyme_kinetics(self):
		tr = TableReader(os.path.join(self._sim_out_path, "EnzymeKinetics"))

		ek = ReaderData()
		reactionIDs = tr.readAttribute('reactionIDs')
		constraintToReaction = tr.readAttribute('constraintToReactionDict')
		constraintIDs = tr.readAttribute('constraintIDs')

		# XXX: This behaviour depends on the sim using this to get its list of metabolite IDs
		# a metabolitePoolIDs attribute appears magically in simData for completed sims, but
		# is not accessible for a running sim :/
		metaboliteIDs = sorted(self._sim_data.process.metabolism.concDict) 

		ek.constraintsLimits = pd.DataFrame(tr.readColumn('allConstraintsLimits'), index=self.time, columns=constraintIDs)
		ek.countsToMolar = pd.DataFrame(tr.readColumn('countsToMolar'), index=self.time)
		ek.metaboliteConcentrations = pd.DataFrame(tr.readColumn('metaboliteConcentrations'), index=self.time, columns=metaboliteIDs)
		ek.metaboliteCountsInitial = pd.DataFrame(tr.readColumn('metaboliteCountsInit'), index=self.time, columns=metaboliteIDs)
		ek.metaboliteCountsFinal = pd.DataFrame(tr.readColumn('metaboliteCountsFinal'), index=self.time, columns=metaboliteIDs)
		ek.reactionRates = pd.DataFrame(tr.readColumn('reactionRates'), index=self.time, columns=reactionIDs)

		self.enzymeKinetics = ek

	def load_fba(self):
		tr = TableReader(os.path.join(self._sim_out_path, "FBAResults"))

		fba = ReaderData()
		reactionIDs = tr.readAttribute('reactionIDs')
		outputMoleculeIDs = tr.readAttribute('outputMoleculeIDs')
		externalMoleculeIDs = tr.readAttribute('externalMoleculeIDs')

		fba.outputFluxes = pd.DataFrame(tr.readColumn('outputFluxes'), index=self.time, columns=outputMoleculeIDs)
		fba.objectiveComponents = pd.DataFrame(tr.readColumn('objectiveComponents'), index=self.time, columns=outputMoleculeIDs)
		fba.externalExchangeFluxes = pd.DataFrame(tr.readColumn('externalExchangeFluxes'), index=self.time, columns=externalMoleculeIDs)
		fba.reactionFluxes = pd.DataFrame(tr.readColumn('reactionFluxes'), index=self.time, columns=reactionIDs)

		self.fbaResults = fba

	def load_growth_limits(self):
		tr = TableReader(os.path.join(self._sim_out_path, "GrowthLimits"))

		gl = ReaderData()
		for col in tr.columnNames():
			if col.startswith('aa'):
				setattr(gl, col, pd.DataFrame(tr.readColumn(col), index=self.time, columns=self._sim_data.moleculeGroups.aaIDs))
			elif col.startswith('ntp'):
				setattr(gl, col, pd.DataFrame(tr.readColumn(col), index=self.time, columns=self._sim_data.moleculeGroups.ntpIds))
			elif col.startswith('gtp'):
				setattr(gl, col, pd.DataFrame(tr.readColumn(col), index=self.time, columns=['GTP[c]']))

		self.growthLimits = gl

	def load_mass(self):
		tr = TableReader(os.path.join(self._sim_out_path, "Mass"))

		mass = ReaderData()
		processNames = tr.readAttribute('processNames')

		for col in tr.columnNames():
			if col.endswith("rocessMassDifferences"):
				setattr(mass, col, pd.DataFrame(tr.readColumn(col), index=self.time, columns=processNames))
			else:
				setattr(mass, col, pd.DataFrame(tr.readColumn(col), index=self.time))

		self.mass = mass

	def load_replication_data(self):
		tr = TableReader(os.path.join(self._sim_out_path, "ReplicationData"))

		rd = ReaderData()
		for col in tr.columnNames():
			setattr(rd, col, pd.DataFrame(tr.readColumn(col), index=self.time))

		self.replicationData = rd

	def load_ribosome_data(self):
		tr = TableReader(os.path.join(self._sim_out_path, "RibosomeData"))

		rd = ReaderData()
		for col in tr.columnNames():
			setattr(rd, col, pd.DataFrame(tr.readColumn(col), index=self.time))

			if col.startswith('aa'):
				getattr(rd, col).columns = self._sim_data.moleculeGroups.aaIDs

		self.ribosomeData = rd

	def load_rna_degradation_data(self):
		tr = TableReader(os.path.join(self._sim_out_path, "RnaDegradationListener"))

		rd = ReaderData()
		for col in tr.columnNames():
			setattr(rd, col, pd.DataFrame(tr.readColumn(col), index=self.time))	

		rd.countRnaDegraded.columns = self._sim_data.process.transcription.rnaData["id"]

		self.rnaDegradation = rd

	def load_rnap_data(self):
		tr = TableReader(os.path.join(self._sim_out_path, "RnapData"))

		rd = ReaderData()
		for col in tr.columnNames():
			if col.startswith('ntp'):
				continue # These are broken

			setattr(rd, col, pd.DataFrame(tr.readColumn(col), index=self.time))

		self.rnapData = rd

	def load_unique_molecule_counts(self):
		tr = TableReader(os.path.join(self._sim_out_path, "UniqueMoleculeCounts"))

		uniqueMoleculeCounts = pd.DataFrame(
			tr.readColumn('uniqueMoleculeCounts'), 
			index=self.time, 
			columns=tr.readAttribute('uniqueMoleculeIds'))

		self.uniqueMoleculeCounts = uniqueMoleculeCounts

	_loaders = {
		'BulkMolecules': load_bulk_molecules,
		'EnzymeKinetics': load_enzyme_kinetics,
		'FBA': load_fba,
		'GrowthLimits': load_growth_limits,
		'Mass': load_mass,
		'ReplicationData': load_replication_data,
		'RibosomeData': load_ribosome_data,
		'RNADegradationData': load_rna_degradation_data,
		'RNAPData': load_rnap_data,
		'UniqueMoleculeCounts': load_unique_molecule_counts
	}


class ReaderData:
	pass