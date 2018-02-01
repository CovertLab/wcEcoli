
"""
Forced Transport Flux Balance Analysis for simplified core network

"""

import numpy as np
import matplotlib.pyplot as plt
import os
from wholecell.utils.modular_fba import FluxBalanceAnalysis

timeTotal = 10
timeStep = 0.5
time = np.linspace(timeStep, timeTotal, num=int(timeTotal/timeStep))

directory = os.path.join('user','forcedFlux')

reactionStoich = {
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

transportID = ['Tc1', 'Tc2', 'Tf', 'Td', 'Te', 'Th', 'To2']
externalExchangedMolecules = ['Carbon1', 'Carbon2', 'Dext', 'Eext', 'Fext', 'Hext', 'Oxygen']
objectiveFunction = {'Biomass': 1} #might be -1

fbaObjectOptions = {
	'reactionStoich' : reactionStoich,
	'externalExchangedMolecules' : externalExchangedMolecules,
	'objective' : objectiveFunction,
	}
fba = FluxBalanceAnalysis(**fbaObjectOptions) # **kwargs allows arbitrary number of arguments to functions

# initialize external molecule levels
fba.externalMoleculeLevelsIs(0.1 * np.ones(len(externalExchangedMolecules)))


# TODO -- iterate FBA over time with timesteps, updating mass and environment.
# TODO -- plot output


objvec = np.empty(len(time))


for i, t in enumerate(time):

	# update transport fluxes
	transportDict = {ID: None for ID in transportID}
	transportDict['Tc1'] = 0.02 #TODO -- this should be michaelis-menten equation

	# constrain transport fluxes
	for ID, value in transportDict.iteritems():
		if value:
			fba.minReactionFluxIs(ID, value)
			fba.maxReactionFluxIs(ID, value)

	# results
	reactionIDs=fba.reactionIDs()
	reactionFluxes = fba.reactionFluxes()
	# flowRates=fba.externalExchangeFluxes() # flux of external nutrients

	indices = [np.where(reactionIDs==id)[0][0] for id in transportDict.keys()]
	transportFluxes = reactionFluxes[indices] # flux across transport reactions

	# exFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * fba.externalExchangeFluxes() / coefficient).asNumber(units.mmol / units.g / units.h)

	fba.outputMoleculeLevelsChange()
	obj = fba.objectiveValue()

	# save values
	objvec[i] = obj






plt.plot(objvec)
plt.xlabel('time')
plt.savefig(os.path.join(directory, 'objvec.png'))



# import ipdb; ipdb.set_trace()



# TODO -- create michaelis-menten function, for which you can initialize kinetic parameters and then get flux given input concentrations


