from wholecell.utils import units
import numpy as np
from wholecell.utils import build_ode
import scipy
import scipy.integrate
import sympy as sp
from wholecell.utils import data



class MetabolismVio(object):
    def __init__(self, raw_data, sim_data):
        # Build the abstractions needed for the metabolism of violacein
        molecules = []  # list of all molecules involved in the pathway

        ratesFwd = []  # rate of reaction
        ratesRev = []
        enzymes = []
        rxnIds = []  # ID tied to each rxn equation

        stoichMatrixI = []  # Molecule indices
        stoichMatrixJ = []  # Reaction indices
        stoichMatrixV = []  # Stoichometric coefficients

        stoichMatrixMass = []  # molecular mass of molecules in stoichMatrixI
        miscrnas_with_singleton_tus = sim_data.getter.get_miscrnas_with_singleton_tus()

        compartment_ids_to_abbreviations = {
            comp['id']: comp['abbrev'] for comp in raw_data.compartments
        }
        for reactionIndex, reaction in enumerate(raw_data.new_gene_data.vioAE.metabolic_reactions_new):
            reactionName = reaction["id"]
            if reactionName not in rxnIds:
                rxnIds.append(reactionName)
                ratesFwd.append(reaction["forward_rate"])
                ratesRev.append(reaction["reverse_rate"])

                try:
                    enzymes.append(reaction["catalyzed_by"][0])
                except:
                    enzymes.append(None)
                reactionIndex = len(rxnIds) - 1

            #Build stoichiometry matrix
            for mol_id, coeff in reaction["stoichiometry"].items():
                # Replace miscRNA subunit IDs with TU IDs
                if mol_id in miscrnas_with_singleton_tus:
                    mol_id = sim_data.getter.get_singleton_tu_id(mol_id)

                mol_id_with_compartment = "{}[{}]".format(
                    mol_id.split('[')[0],
                    compartment_ids_to_abbreviations[mol_id.split('[')[1][:-1]]
                )

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
                molecularMass = sim_data.getter.get_mass(mol_id_with_compartment.split('[')[0]).asNumber(units.g / units.mol)
                stoichMatrixMass.append(molecularMass)

        self.molecule_names = np.array(molecules, dtype='U')
        self.rxn_ids = rxnIds
        self.enzymes = enzymes
        self.rates_fwd = np.array(ratesFwd)
        self.rates_rev = np.array(ratesRev)

        self._stoichMatrixI = np.array(stoichMatrixI)
        self._stoichMatrixJ = np.array(stoichMatrixJ)
        self._stoichMatrixV = np.array(stoichMatrixV)

        self._stoich_matrix_mass = np.array(stoichMatrixMass)
        self.balance_matrix = self.stoich_matrix() * self.mass_matrix()

        # Find the mass balance of each equation in the balanceMatrix

        massBalanceArray = self.mass_balance()


        # The stoichometric matrix should balance out to numerical zero. This gives an error with smaller number.
        assert np.max(np.absolute(massBalanceArray)) < 1e-2

        # Build matrices
        self._populate_derivative()

    def __getstate__(self):
        """Return the state to pickle, omitting derived attributes that
        __setstate__() will recompute, esp. the ode_derivatives
        that don't pickle.
        """
        return data.dissoc_strict(self.__dict__, (
            'symbolic_rates',
            '_rates',
            '_stoich_matrix'))

    def __setstate__(self, state):
        """Restore instance attributes, recomputing some of them."""
        self.__dict__.update(state)
        self._populate_derivative()

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

    def _populate_derivative(self):
        '''Compile callable functions for computing the derivative.'''
        self._make_derivative()

        self._rates = build_ode.derivatives(self.symbolic_rates)
        self._stoich_matrix = self.stoich_matrix()  # Matrix is small and can be cached for derivatives
        # WORKAROUND: Avoid Numba LoweringError JIT-compiling these functions:

    def _make_y_dy(self):
        '''This is building the right hand side of the ODE system using mass action kinetics.'''

        S = self.stoich_matrix()

        yStrings = ["y[%d]" % x for x in range(S.shape[0])]

        y = sp.symbols(yStrings)

        rates = []
        for colIdx in range(S.shape[1]):
            negIdxs = np.where(S[:, colIdx] < 0)[0]
            posIdxs = np.where(S[:, colIdx] > 0)[0]

            reactantFlux = self.rates_fwd[colIdx]
            for negIdx in negIdxs:
                stoich = -S[negIdx, colIdx]
                if stoich == 1:
                    if self.enzymes[colIdx]!=[]:
                        reactantFlux *= y[negIdx]
                else:
                    reactantFlux *= y[negIdx] ** stoich

            productFlux = self.rates_rev[colIdx]
            for posIdx in posIdxs:
                stoich = S[posIdx, colIdx]
                if stoich == 1:
                    productFlux *= y[posIdx]
                else:
                    productFlux *= y[posIdx] ** stoich

            rates.append(reactantFlux - productFlux)
        return y, rates

    def molecules_to_next_time_step(self, moleculeCounts, enzymeNames, enzymeDict, cellVolume,
                                    nAvogadro, timeStepSec, random_state, method, min_time_step=None,
                                    jit=True):
        """
        Calculates the changes in the counts of molecules in the next timestep
        by solving an initial value ODE problem.

        Args:
            moleculeCounts (1d ndarray, ints): current counts of molecules
                involved in the ODE
            cellVolume (float): current volume of cell
            nAvogadro (float): Avogadro's number
            timeStepSec (float): current length of timestep in seconds
            random_state (RandomState object): process random state
            method (str): name of the ODE method to use
            min_time_step (int): if not None, timeStepSec will be scaled down until
                it is below min_time_step if negative counts are encountered
            jit (bool): if True, use the jit compiled version of derivatives
                functions
            methods_tried (Optional[Set[str]]): methods for the solver that have
                already been tried

        Returns:
            moleculesNeeded (1d ndarray, ints): counts of molecules that need
                to be consumed
            allMoleculesChanges (1d ndarray, ints): expected changes in
                molecule counts after timestep
        """
        #enzymeC = []
#
        #for colIdx in range(self._stoich_matrix.shape[1]):
        #    if enzymeNames[colIdx] is None:
        #        enzymeC.append(1)
        #    else:
        #        enzymeC.append(enzymeCounts[ind])
        #print(enzymeC)

        enzymeC = []
        for e in enzymeNames:
            if e is None:
                enzymeC.append(float(1))
            else:
                enzymeC.append(float(enzymeDict[e+'[c]']))

        y_init = moleculeCounts / (cellVolume * nAvogadro)
        # TO DO Ioana: Fix the jit version. At the moment it is identical to the non-jit
        if jit:
            derivatives = self.derivatives_jit
        else:
            derivatives = self.derivatives

        sol = scipy.integrate.solve_ivp(
            lambda t, y: derivatives([0, timeStepSec], y_init, enzymeC), [0, timeStepSec], y_init,
            method=method, t_eval=[0, timeStepSec], atol=1e-8)
        y = sol.y.T

        y[y < 0] = 0
        yMolecules = y * (cellVolume * nAvogadro)
        moleculesNeeded = abs(np.sum(self.stoich_matrix().clip(max=0), axis=1))

        # Calculate changes in molecule counts for all molecules
        allMoleculesChanges = yMolecules[-1, :]

        return moleculesNeeded, allMoleculesChanges

    def derivatives(self, t, y, enzymeC):
        """
        Calculate derivatives from stoichiometry and rates with argument order
        for solve_ivp.
        """

        return self._stoich_matrix.dot(np.multiply(enzymeC, self._rates[0](y, t)))


    def derivatives_jit(self, t, y, enzymeC):
        """
        Calculate derivatives from stoichiometry and rates with argument order
        for solve_ivp.
        """
        #print('Current conc ', y)
        #print('S ', self._stoich_matrix)
        #print('enzyme ', enzymeC)
        print('rates ', self._rates[0](y, t))
        return self._stoich_matrix.dot(np.multiply(enzymeC, self._rates[0](y, t)))

    def _make_derivative(self):
        '''
        Creates symbolic representation of the ordinary differential equations.
         Used during simulations.
        '''
        y, rates = self._make_y_dy()
        rates = sp.Matrix(rates)
        self.symbolic_rates = rates
