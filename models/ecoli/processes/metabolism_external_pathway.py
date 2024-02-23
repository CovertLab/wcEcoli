"""
Metabolism sub-model for new(external) pathway.
This is a separate model to E.coli's metabolism and it is only run if the user has added an external
pathway to E. coli's metabolism.
The corresponding reactions of this model are catalysed by the enzymes produced by the exogenous genes
that we are adding.
The reactions are defined in /wcEcoli/reconstruction/ecoli/flat/new_gene_data.
"""

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

TIME_UNITS = units.s
COUNTS_UNITS = units.mmol
MASS_UNITS = units.g
TIME_UNITS = units.s
VOLUME_UNITS = units.L
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS
CONVERSION_UNITS = MASS_UNITS * TIME_UNITS / VOLUME_UNITS
GDCW_BASIS = units.mmol / units.g / units.h
MODEL_FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

class MetabolismExternalPathway(wholecell.processes.process.Process):
    """ Metabolism External Pathway """

    _name = "MetabolismExternalPathway"

    # Constructor
    def __init__(self):
        super(MetabolismExternalPathway, self).__init__()

    # Construct object graph
    def initialize(self, sim, sim_data):
        super(MetabolismExternalPathway, self).initialize(sim, sim_data)
        # Simulation options
        self.jit = sim._jit

        # Get constants
        self.n_avogadro = sim_data.constants.n_avogadro.asNumber(1 / units.mmol)
        self.cell_density = sim_data.constants.cell_density.asNumber(units.g / units.L)

        self.external_pathway_exists = sim_data.process.metabolism_external_pathway.has_external_pathway

        # check if there is a new metabolic pathway to be added
        if self.external_pathway_exists:
            # Create method
            self.molecules_to_next_time_step = sim_data.process.metabolism_external_pathway.molecules_to_next_time_step

            # Build views
            self.molecule_names = sim_data.process.metabolism_external_pathway.molecule_names
            self.molecules = self.bulkMoleculesView(self.molecule_names)

            enzyme_names = sum(sim_data.process.metabolism_external_pathway.enzymes.values(), [])
            self.clean_enzyme = sorted([*set(list(filter(lambda item: item is not None, enzyme_names)))])
            self.enzymes = self.bulkMoleculesView(self.clean_enzyme)

            self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)


    def calculateRequest(self):
        if self.external_pathway_exists:
            # Get molecule counts
            molecule_counts = self.molecules.total_counts()
            self.molecule_dict = dict(zip(self.molecule_names, molecule_counts))
            enzyme_counts = self.enzymes.total_counts()
            self.enzymes_dict = dict(zip(self.clean_enzyme, enzyme_counts))

            # Get cell mass and volume
            self.dry_mass = (self.readFromListener("Mass", "dryMass") * units.fg).asNumber(units.g)
            self.cell_mass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
            self.cell_volume = self.cell_mass / self.cell_density
            self.dcw_to_volume = self.cell_density * (self.dry_mass / self.cell_mass)
            self.molecules_required, self.all_molecule_changes, self.flux = self.molecules_to_next_time_step(
                self.molecule_dict, self.enzymes_dict,
                self.cell_volume, self.n_avogadro, self.timeStepSec(), method="LSODA")

            # Request counts of molecules needed
            self.molecules.requestIs(self.molecules_required)

    def evolveState(self):
        if self.external_pathway_exists:
            self.writeToListener('FBAResults', 'externalPathwayFluxes', self.flux)
            # Increment changes in molecule counts
            self.molecules.countsInc(self.all_molecule_changes)
