"""
Metabolism for new(external) pathway
Metabolism for new pathway sub-model

"""
from __future__ import absolute_import, division, print_function

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM_NEW


class MetabolismNewPathway(wholecell.processes.process.Process):
    """ Metabolism New Pathway """

    _name = "MetabolismNewPathway"

    # Constructor
    def __init__(self):
        super(MetabolismNewPathway, self).__init__()

    # Construct object graph
    def initialize(self, sim, sim_data):
        super(MetabolismNewPathway, self).initialize(sim, sim_data)

        # Simulation options
        self.jit = sim._jit

        # Get constants
        self.nAvogadro = sim_data.constants.n_avogadro.asNumber(1 / units.mmol)
        self.cellDensity = sim_data.constants.cell_density.asNumber(units.g / units.L)

        #check if there is a new metabolic pathway to be added
        if hasattr(sim_data.process.metabolism_new_pathway, 'molecule_names'):
            self.new_pathway = True
            # Create method
            self.moleculesToNextTimeStep = sim_data.process.metabolism_new_pathway.molecules_to_next_time_step

            # Build views
            self.moleculeNames = sim_data.process.metabolism_new_pathway.molecule_names
            self.molecules = self.bulkMoleculesView(self.moleculeNames)

            self.enzymeNames = sim_data.process.metabolism_new_pathway.enzymes
            cleanEnzyme = [*set(list(filter(lambda item: item is not None, self.enzymeNames)))]
            self.cleanEnzyme = [x+'[c]' for x in cleanEnzyme]
            self.enzymes = self.bulkMoleculesView(self.cleanEnzyme)

            # Set priority to a lower value (but greater priority than metabolism)
            self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM_NEW)
        else:
            self.new_pathway = False


    def calculateRequest(self):
        if self.new_pathway:
            # Get molecule counts
            moleculeCounts = self.molecules.total_counts()
            enzymeCounts = self.enzymes.total_counts()
            EnzymesDict = dict(zip(self.cleanEnzyme, enzymeCounts))

            # Get cell mass and volume
            cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
            self.cellVolume = cellMass / self.cellDensity


            # Solve ODEs to next time step using the BDF solver through solve_ivp.
            # Note: the BDF solver has been empirically tested to be the fastest
            # solver for this setting among the list of solvers that can be used
            # by the scipy ODE suite.
            self.molecules_required, self.all_molecule_changes = self.moleculesToNextTimeStep(
            	moleculeCounts, self.moleculeNames, self.enzymeNames, EnzymesDict, self.cellVolume, self.nAvogadro,
            	self.timeStepSec(), method="LSODA", jit=self.jit)

            # Request counts of molecules needed
            self.molecules.requestIs(self.molecules_required)


    def evolveState(self):
        if self.new_pathway:
            # Get counts of molecules allocated to this process
            moleculeCounts = self.molecules.counts()
            enzymeCounts = self.enzymes.total_counts()
            EnzymesDict = dict(zip(self.cleanEnzyme, enzymeCounts))


            # Check if any molecules were allocated fewer counts than requested
            if (self.molecules_required > moleculeCounts).any():
                # Solve ODEs to a large time step using the the counts of molecules
                # allocated to this process using the BDF solver for stable integration.
                # The number of reactions has already been determined in calculateRequest
                # and rates will be much lower with a fraction of total counts allocated
                # so a much longer time scale is needed.

                _, self.all_molecule_changes = self.moleculesToNextTimeStep(
                    moleculeCounts, self.moleculeNames, self.enzymeNames, EnzymesDict, self.cellVolume, self.nAvogadro,
                    10000, method="BDF", jit=self.jit
                    )

            # Increment changes in molecule counts
            self.molecules.countsInc(self.all_molecule_changes)

