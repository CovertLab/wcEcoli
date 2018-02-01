
"""
Flux Balance Analysis for simplified core network

"""

import numpy as np
from wholecell.utils.modular_fba import FluxBalanceAnalysis

reactionStoich={
	# Metabolic reactions
	'R1': {'A': -1, 'ATP': -1, 'B': 1},
	'R2a': {'B': -1, 'ATP': 2, 'NADH': 2, 'C': 1},
	'R2b': {'C': -1, 'ATP': -2, 'NADH': -2, 'B': 1},
	'R3': {'B': -1, 'F': -1},
	'R4': {'C': -1, 'G': 1},
	'R5a': {'G': -1, 'C': 0.8, 'NADH': 2},
	'R5b': {'G': -1, 'C': 0.8, 'NADH': 2},
	'R6': {'C': -1, 'ATP': 2, 'D': 3},
	'R7': {'C': -1, 'NADH': -4, 'E': 3},
	'R8a': {'G': -1, 'ATP': -1, 'NADH': -2, 'H': 1},
	'R8b': {'G': 1, 'ATP': 1, 'NADH': 2, 'H': -1},
	'Rres': {'NADH': -1, 'O2': -1, 'ATP': 1},

	# Transport processes
	'Tc1': {'Carbon1': -1, 'A': 1},
	'Tc2': {'Carbon2': -1, 'A': 1},
	'Tf': {'Fext': -1, 'F': 1},
	'Td': {'D': -1, 'Dext': 1},
	'Te': {'E': -1, 'Eext': 1},
	'Th': {'Hext': -1, 'H': 1},
	'To2': {'Oxygen': -1, 'O2': 1},

	# Maintenance and growth processes
	'Growth': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10, 'Biomass': 1}
	}

externalExchangedMolecules=['Carbon1', 'Carbon2', 'Dext', 'Eext', 'Fext', 'Hext', 'Oxygen']

objectiveFunction = {'Biomass': 1} #might be -1

fbaObjectOptions = {
	'reactionStoich' : reactionStoich,
	'externalExchangedMolecules' : externalExchangedMolecules,
	'objective' : objectiveFunction,
	}

fba = FluxBalanceAnalysis(**fbaObjectOptions) # **kwargs allows arbitrary number of arguments to functions
fba.externalMoleculeLevelsIs(0.1 * np.ones(len(externalExchangedMolecules)))


# determine flux across transporters
transportID = 'Tc2'
transportFlux = 0.01


fba.minReactionFluxIs(transportID, transportFlux)
fba.maxReactionFluxIs(transportID, transportFlux)

# results
reactionIDs=fba.reactionIDs()
reactionFluxes = fba.reactionFluxes()
transportFluxes=fba.externalExchangeFluxes()

index = np.where(reactionIDs==transportID)[0][0]
getFlux = reactionFluxes[index]



fba.outputMoleculeLevelsChange()

obj = fba.objectiveValue()





import ipdb; ipdb.set_trace()




