from wholecell.utils import units
import numpy as np
from wholecell.utils import build_ode
import scipy
import scipy.integrate
import sympy as sp
from wholecell.utils import data

#TODO: reactions catalyzed by several enzymes; reactions that have several substrates; spontaenous reactions have to be handled differently

# Alternative methods to try (in order of priority) when solving ODEs to the next time step
IVP_METHODS = ['LSODA']

class MetabolismExternalPathway(object):
    def __init__(self, raw_data, sim_data):
        if raw_data.new_genes_option != "off":
            # Build the abstractions needed for the metabolism corresponding to the external pathway
            molecules = []  # list of all molecules involved in the pathway
            kcats = []
            kms = []
            substrates = []
            rxnIds = []  # ID tied to each rxn equation
            enzymes = {}  # dictionary of all enzymes for each reaction

            stoichMatrixI = []  # Molecule indices
            stoichMatrixJ = []  # Reaction indices
            stoichMatrixV = []  # Stoichometric coefficients

            stoichMatrixMass = []  # molecular mass of molecules in stoichMatrixI
            independentMolecules = []  # list of all specific independent molecule names
            independent_molecule_indexes = []  # index of each of the independent molecules

            new_genes_folder = getattr(raw_data.new_gene_data, raw_data.new_genes_option)
            for reactionIndex, reaction in enumerate(new_genes_folder.metabolic_reactions_external):
                reactionName = reaction["id"]
                if reactionName not in rxnIds:
                    rxnIds.append(reactionName)
                    kcats.append(reaction["kcat"])
                    kms.append(reaction["km"])
                    substrates.append(reaction["substrate"])

                    enzymes[reactionName] = reaction["catalyzed_by"]

                #Build stoichiometry matrix
                for mol_id, coeff in reaction["stoichiometry"].items():
                    mol_id_with_compartment = mol_id

                    if mol_id_with_compartment not in molecules:
                        molecules.append(mol_id_with_compartment)
                        molecule_index = len(molecules) - 1
                    else:
                        molecule_index = molecules.index(mol_id_with_compartment)

                    assert (coeff % 1) == 0

                    stoichMatrixI.append(molecule_index)
                    stoichMatrixJ.append(reactionIndex)
                    stoichMatrixV.append(coeff)

                    # Find molecular mass of the molecule and add to mass matrix
                    molecularMass = sim_data.getter.get_mass(mol_id_with_compartment[:-3]).asNumber(units.g / units.mol)
                    stoichMatrixMass.append(molecularMass)

                self.molecule_names = np.array(molecules, dtype='U')
                self.rxn_ids = rxnIds
                self.flux = np.array(np.zeros(len(rxnIds)))
                self.enzymes = enzymes
                self.kcats = np.array(kcats)
                self.kms = np.array(kms)
                self.substrates = np.array(substrates)

                self.independent_molecules = np.array(independentMolecules, dtype='U')
                self.independent_molecule_indexes = np.array(independent_molecule_indexes)
                self._stoichMatrixI = np.array(stoichMatrixI)
                self._stoichMatrixJ = np.array(stoichMatrixJ)
                self._stoichMatrixV = np.array(stoichMatrixV)

                self._stoich_matrix_mass = np.array(stoichMatrixMass)
                self.balance_matrix = self.stoich_matrix() * self.mass_matrix()

                # Find the mass balance of each equation in the balanceMatrix
                massBalanceArray = self.mass_balance()

                # The stoichometric matrix should balance out to numerical zero. This gives an error with smaller number.
                assert np.max(np.absolute(massBalanceArray)) < 1e-9

                # Build matrices
                self._populate_derivative_and_jacobian()

            if rxnIds==[]:
                print('There are external genes, but no external pathway')
                self.has_external_pathway = False
            else:
                print('There is an external pathway and external genes.')
                self.has_external_pathway = True
        else:
            print('There is no external pathway and no external genes.')
            self.has_external_pathway = False

    def __getstate__(self):
        """Return the state to pickle, omitting derived attributes that
        __setstate__() will recompute, esp. the ode_derivatives
        that don't pickle.
        """
        try:
            return data.dissoc_strict(self.__dict__, (
                'symbolic_rates',
			    '_rates',
			    '_stoich_matrix'))
        except:
            return self.__dict__


    def __setstate__(self, state):
        """Restore instance attributes, recomputing some of them."""
        try:
            self.__dict__.update(state)
            self._populate_derivative_and_jacobian()
        except:
            return 0

    def stoich_matrix(self):
        '''
        Builds stoichiometry matrix
        Rows: molecules
        Columns: reactions
        Values: reaction stoichiometry
        '''
        shape = (self._stoichMatrixI.max() + 1, self._stoichMatrixJ.max() + 1)
        out = np.zeros(shape, np.float64)
        out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixV
        return out

    def mass_matrix(self):
        '''
        Builds stoichiometry mass matrix
        Rows: molecules
        Columns: reactions
        Values: molecular mass
        '''
        shape = (self._stoichMatrixI.max() + 1, self._stoichMatrixJ.max() + 1)
        out = np.zeros(shape, np.float64)
        out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoich_matrix_mass
        return out

    def mass_balance(self):
        '''
        Sum along the columns of the massBalance matrix to check for reaction mass balance
        '''
        return np.sum(self.balance_matrix, axis=0)

    def _populate_derivative_and_jacobian(self):
        '''Compile callable functions for computing the derivative.'''
        self._make_derivative()

        self._rates = build_ode.derivatives(self.symbolic_rates)
        self._stoich_matrix = self.stoich_matrix()  # Matrix is small and can be cached for derivatives


    def _make_y_dy(self):
        '''This is building the right hand side of the ODE system using MM kinetics.'''

        substrates = self.substrates
        yStrings = ["y[%d]" % x for x in range(len(substrates))]
        y = sp.symbols(yStrings)
        rates = []
        for colIdx in range(len(substrates)):
            flux = self.kcats[colIdx] * y[colIdx]/ (self.kms[colIdx] + y[colIdx])
            rates.append(flux)
        return y, rates

    def _make_derivative(self):
        '''
        Creates symbolic representation of the ordinary differential equations
        and the Jacobian. Used during simulations.
        '''
        y, rates = self._make_y_dy()

        rates = sp.Matrix(rates)
        self.symbolic_rates = rates


    def molecules_to_next_time_step(self, moleculeDict, enzymeDict, cellVolume,
                                    nAvogadro, timeStepSec, method="LSODA"):

        """
        Calculates the changes in the counts of molecules in the next timestep
        by solving an initial value ODE problem.

        Args:
            moleculeDict (dict): name and current counts of molecules
                involved in the ODE
            enzymeDict (dict): dictionary of enzyme counts
            cellVolume (float): current volume of cell
            nAvogadro (float): Avogadro's number
            timeStepSec (float): current length of timestep in seconds
            method (str): name of the ODE method to use
            jit (bool): if True, use the jit compiled version of derivatives
                functions


        Returns:
            moleculesNeeded (1d ndarray, ints): counts of molecules that need
                to be consumed
            allMoleculesChanges (1d ndarray, ints): expected changes in
                molecule counts after timestep
        """

        subCounts, subNames, enzymeCounts, enzymeNames = [], [], [], []
        # get the counts of the substrates.
        # we assume that each reaction has only one substrate
        for sub in self.substrates:
            subCounts.append(moleculeDict[sub[0]])
            subNames.append(sub[0])

        # get the counts of the enzymes
        # we assume that each reaction has only one enzyme
        for reaction in self.rxn_ids:
            # if the reaction has an enzyme
            if self.enzymes[reaction]:
                enzymeCounts.append(enzymeDict[self.enzymes[reaction][0]])
                enzymeNames.append(self.enzymes[reaction][0])
            #If the reaction is spontaneous, we set the enzyme counts to a high value to avoid creaitng a bottleneck.
            #This might need to be modelled differently
            else:
                enzymeCounts.append(100000)
                enzymeNames.append('None')


        enzymeConc = np.array(enzymeCounts) / (cellVolume * nAvogadro)
        sub_init = np.array(subCounts) / (cellVolume * nAvogadro)

        moleculeCounts = list(moleculeDict.values())
        y_init = np.array(moleculeCounts) / (cellVolume * nAvogadro)
        derivatives = self.derivatives

        sol = scipy.integrate.solve_ivp(
            lambda t, y: derivatives([0, timeStepSec], sub_init, enzymeConc), [0, timeStepSec], y_init,
            method=method, t_eval=[0, timeStepSec], atol=1e-8)
        y = sol.y.T
        y[y < 0] = 0
        yMolecules = (y * (cellVolume * nAvogadro)) #.asNumber(COUNTS_UNITS/VOLUME_UNITS)

        # Calculate changes in molecule counts for all molecules
        allMoleculesChanges = yMolecules[-1, :] - yMolecules[0, :]

        # Molecules needed are the change but only for the molecules corresponding to the negative terms in the stoich matrix
        moleculesNeeded = np.negative(allMoleculesChanges).clip(min=0)

        return moleculesNeeded, allMoleculesChanges, self.flux

    def derivatives(self, t, y, enzymeC):
        """
        Calculate derivatives from stoichiometry and rates with argument order
        for solve_ivp.
        """
        self.flux =  np.multiply(enzymeC, self._rates[0](y, t))
        return self._stoich_matrix.dot(np.transpose(self.flux))
