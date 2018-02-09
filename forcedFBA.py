
"""
Forced Transport Flux Balance Analysis for simplified core network

"""

import numpy as np
import scipy.constants
import matplotlib.pyplot as plt
import os
from wholecell.utils.modular_fba import FluxBalanceAnalysis
# from wholecell.utils import units

NO_TRANSPORT = False
TRANSPORT_CONSTRAINT = True # transporter fluxes are set as constraints on FBA optimization
TRANSPORT_OBJECTIVE = False # transporter fluxes are set as objectives with a penalty
TRANSPORT_PROCESS = False  # transporter fluxes set internal concentrations, which are pulled from by FBA

#michaelis-menten function, gives flux of reaction
def michaelisMenten(Kcat, Km, substrate, enzyme):
	flux = enzyme * substrate * Kcat / (Km + substrate)
	return flux

directory = os.path.join('user','forcedFlux')

timeTotal = 50.0
timeStep = 0.01
time = np.linspace(timeStep, timeTotal, num=int(timeTotal/timeStep))

reactionStoich = {
	#Metabolic reactions
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

	#Transport processes
	'Tc1': {'Carbon1': -1, 'A': 1},
	'Tc2': {'Carbon2': -1, 'A': 1},
	'Tf': {'Fext': -1, 'F': 1},
	'Td': {'D': -1, 'Dext': 1},
	'Te': {'E': -1, 'Eext': 1},
	'Th': {'Hext': -1, 'H': 1},
	'To2': {'Oxygen': -1, 'O2': 1},

	#Maintenance and growth processes
	'Growth': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10, 'Biomass': 1}
	}

transportID = ['Tc1', 'Tc2', 'Tf', 'Td', 'Te', 'Th', 'To2']
metaboliteID = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'O2', 'ATP', 'NADH']
externalExchangedMolecules = ['Carbon1', 'Carbon2', 'Dext', 'Eext', 'Fext', 'Hext', 'Oxygen']

objectiveFunction = {'Biomass': 1} #might be -1

transportKinetics = {
	'Tc1': {'Kcat': 0.5, 'Km': 0.1},
	'Tc2': {'Kcat': 0.8, 'Km': 2},
	'Tf': {'Kcat': 5.0, 'Km': 2.5},
	'Td': {'Kcat': 10.0, 'Km': 2.5},
	'Te': {'Kcat': 10.0, 'Km': 2.5},
	'Th': {'Kcat': 10.0, 'Km': 2.5},
	'To2': {'Kcat': 10.0, 'Km': 2.5},
	}

transportConstraint = ['Tc1','Tc2'] # list of constrained transporters TODO -- put this in transportKinetics.

fbaObjectOptions = {
	'reactionStoich' : reactionStoich,
	'externalExchangedMolecules' : externalExchangedMolecules,
	'objective' : objectiveFunction,
	}

#initialize fba
fba = FluxBalanceAnalysis(**fbaObjectOptions)

#initialize state
cellDensity = 1100.0 #g/L
initdryMass = 400.0 # fg = 1e-15 g
initnutrientLevels = np.full((len(externalExchangedMolecules)), 100) + np.random.normal(0,2,len(externalExchangedMolecules)) # add gaussian noise







#transport kinetics are not forced
if NO_TRANSPORT:

	dryMass = initdryMass
	nutrientLevels = initnutrientLevels

	fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

	objVec = np.empty(len(time))
	bioMassVec = np.empty(len(time))
	nutrientVec = np.empty([len(time),len(externalExchangedMolecules)])

	# iterate FBA over time, updating mass and environment. fba.outputMoleculeLevelsChange() MAKES FBA UPDATE!
	for i, t in enumerate(time):
		fba.externalMoleculeLevelsIs(nutrientLevels.tolist())
		dryMassOld = dryMass
		bXout = fba.outputMoleculeLevelsChange() # this is flux of biomass
		# fba.biomassReactionFlux() # is this biomass exchange flux?

		# when bXout approaches 0, it becomes nan. 
		if bXout == bXout and bXout != 0:
			dryMass = dryMassOld * np.exp(bXout * (timeStep/360))
			nutrientLevels = nutrientLevels + fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld)

		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = nutrientLevels.tolist()


		assert sum(dryMass-dryMassOld) == sum(nutrientLevels-nutrientLevelsOld) # these should be equal.. right?
		# import ipdb; ipdb.set_trace()

	plt.figure(4)
	plt.plot(time, objVec)
	plt.xlabel('time (s)')
	plt.title('objective')
	plt.savefig(os.path.join(directory, 'noforce_objective.png'))

	plt.figure(5)
	plt.plot(time, bioMassVec)
	plt.xlabel('time (s)')
	plt.title('biomass')
	plt.savefig(os.path.join(directory, 'noforce_biomass.png'))

	plt.figure(6)
	for n in range(len(nutrientLevels)):
		plt.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	plt.xlabel('time (s)')
	plt.ylabel('concentration')
	plt.legend()
	plt.title('external nutrient levels')
	plt.savefig(os.path.join(directory, 'noforce_externalNutrients.png'))










#transport kinetics are set as constraints
if TRANSPORT_CONSTRAINT:

	dryMass = initdryMass
	nutrientLevels = initnutrientLevels

	fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

	objVec = np.empty(len(time))
	bioMassVec = np.empty(len(time))
	nutrientVec = np.empty([len(time),len(externalExchangedMolecules)])

	# iterate FBA over time, updating mass and environment. fba.outputMoleculeLevelsChange() MAKES FBA UPDATE!
	for i, t in enumerate(time):

		fba.externalMoleculeLevelsIs(nutrientLevels.tolist())


		#constrain transport flux
		transportDict = {ID: None for ID in transportID}
		for index, ID in enumerate(transportConstraint): #
			Kcat = transportKinetics[ID]['Kcat']
			Km = transportKinetics[ID]['Km']

			#get this transporter's substrate (value = -1) from reactionStoic
			substrateName = reactionStoich[ID].keys()[reactionStoich[ID].values().index(-1)]
			nutrientIndex = externalExchangedMolecules.index(substrateName)
			substrateConc = nutrientLevels[nutrientIndex]

			# TODO -- add vector for transporter concentrations
			transportDict[ID] = michaelisMenten(Kcat, Km, substrateConc, 1) #transporter concentration set to 1

		# Maybe I should set the following:
		# fba.setMinReactionFluxes
		# fba.setMaxReactionFluxes
		for ID, value in transportDict.iteritems():
			if value:
				fba.minReactionFluxIs(ID, value)
				fba.maxReactionFluxIs(ID, value)

		# import ipdb; ipdb.set_trace()


		dryMassOld = dryMass
		bXout = fba.outputMoleculeLevelsChange() # this is flux of biomass
		
		# when bXout approaches 0, it becomes nan. 
		if bXout == bXout and bXout != 0:
			dryMass = dryMassOld * np.exp(bXout * (timeStep/360))
			nutrientLevels = nutrientLevels + fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld) 



		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = nutrientLevels.tolist()

	#plot
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	#wspace = width between subplots, hspace = height between subplots
	f.subplots_adjust(wspace=0.4, hspace=0.5)
	ax1.plot(time, objVec)
	ax1.set(xlabel='time (s)')
	ax1.set_title('objective')

	ax2.plot(time, bioMassVec)
	ax2.set(xlabel='time (s)')
	ax2.set_title('biomass')

	for n in range(len(nutrientLevels)):
		ax3.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	ax3.set(xlabel='time (s)', ylabel='concentration')
	ax3.set_title('external nutrient levels')
	ax3.legend(prop={'size': 6})

	substrate=np.linspace(0.0, 50.0, num=1001)
	for ID, kinetics in transportKinetics.iteritems():
		nflux=np.empty_like(substrate)
		for ind, sConc in enumerate(substrate):
			Kcat = kinetics['Kcat']
			Km = kinetics['Km']
			nflux[ind]=michaelisMenten(Kcat, Km, sConc, 1)
		ax4.plot(substrate, nflux, label=ID)
	ax4.set(xlabel='substrate concentration', ylabel='flux')
	ax4.set_title('transporters flux range')
	ax4.legend(prop={'size': 6})
	f.savefig(os.path.join(directory, 'constraint_transport.png'), dpi=200)





#transport kinetics are set as objectives
if TRANSPORT_OBJECTIVE:

	#enable kinetic targets
	fba.enableKineticTargets()
	# TODO -- how are these Targets set?

	dryMass = initdryMass
	nutrientLevels = initnutrientLevels

	fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

	objVec = np.empty(len(time))
	bioMassVec = np.empty(len(time))
	nutrientVec = np.empty([len(time),len(externalExchangedMolecules)])

	# iterate FBA over time, updating mass and environment. fba.outputMoleculeLevelsChange() MAKES FBA UPDATE!
	for i, t in enumerate(time):

		fba.externalMoleculeLevelsIs(nutrientLevels.tolist())


		#transport flux targets
		transportDict = {ID: None for ID in transportID}
		for index, ID in enumerate(transportConstraint): #
			Kcat = transportKinetics[ID]['Kcat']
			Km = transportKinetics[ID]['Km']

			#get this transporter's substrate (value = -1) from reactionStoic
			substrateName = reactionStoich[ID].keys()[reactionStoich[ID].values().index(-1)]
			nutrientIndex = externalExchangedMolecules.index(substrateName)
			substrateConc = nutrientLevels[nutrientIndex]

			# TODO -- add vector for transporter concentrations
			transportDict[ID] = michaelisMenten(Kcat, Km, substrateConc, 0.1) #transporter concentration set to 1

			targets = np.ones(len(transportConstraint)) # this should be a numpy array of len(transportConstraint) with target values




		'''
		TODO -- set kinetic flux targets according to determined transport fluxes
		'''
		for ID, value in transportDict.iteritems():
			if value:
				# fba.kineticTargetFluxTargets
				# self.fba.setKineticTarget(self.kineticsConstrainedReactions, targets, raiseForReversible = False)
				fba.setKineticTarget(transportConstraint, targets)


		dryMassOld = dryMass
		bXout = fba.outputMoleculeLevelsChange() # this is flux of biomass
		
		# when bXout approaches 0, it becomes nan. 
		if bXout == bXout and bXout != 0:
			dryMass = dryMassOld * np.exp(bXout * (timeStep/360))
			nutrientLevels += fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld) #* (timeStep/360) # TODO -- is the timestep here correct?

		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = nutrientLevels.tolist()


	#disable kinetic targets once simulation is done.
	fba.disableKineticTargets()

	plt.figure(7)
	plt.plot(time, objVec)
	plt.xlabel('time (s)')
	plt.title('objective')
	plt.savefig(os.path.join(directory, 'obj_objective.png'))

	plt.figure(8)
	plt.plot(time, bioMassVec)
	plt.xlabel('time (s)')
	plt.title('biomass')
	plt.savefig(os.path.join(directory, 'obj_biomass.png'))

	plt.figure(9)
	for n in range(len(nutrientLevels)):
		plt.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	plt.xlabel('time (s)')
	plt.ylabel('concentration')
	plt.legend()
	plt.title('external nutrient levels')
	plt.savefig(os.path.join(directory, 'obj_externalNutrients.png'))










#transport kinetics are set as objectives
if TRANSPORT_PROCESS:

	#initialize stores for added step
	externalMolecules = ['Carbon1_out', 'Carbon2_out', 'Dext_out', 'Eext_out', 'Fext_out', 'Hext_out', 'Oxygen_out']
	dryMass = initdryMass
	externalNutrientLevels = initnutrientLevels
	internalNutrientLevels = np.empty_like(initnutrientLevels)

	objVec = np.empty(len(time))
	bioMassVec = np.empty(len(time))
	nutrientVec = np.empty([len(time),len(externalExchangedMolecules)])
	externalVec = np.empty([len(time),len(externalExchangedMolecules)])

	# iterate FBA over time, updating mass and environment. fba.outputMoleculeLevelsChange() MAKES FBA UPDATE!
	for i, t in enumerate(time):

		#transport processes to set downstream limits
		transportDict = {ID: None for ID in transportID}
		for index, ID in enumerate(transportConstraint): #
			Kcat = transportKinetics[ID]['Kcat']
			Km = transportKinetics[ID]['Km']

			#get this transporter's substrate (value = -1) from reactionStoic
			substrateName = reactionStoich[ID].keys()[reactionStoich[ID].values().index(-1)]
			nutrientIndex = externalExchangedMolecules.index(substrateName)
			substrateConc = externalNutrientLevels[nutrientIndex]

			# TODO -- add vector for transporter concentrations
			transportDict[ID] = michaelisMenten(Kcat, Km, substrateConc, 0.1) #transporter concentration set to 1

		# 

		# determine internal nutrient levels based on transport
		# TODO -- if value is "None" (no transport), need to use externalNutrient
		deltaNutrients = transportDict.values()

		for index, val in enumerate(deltaNutrients):
			if val is not None:
				externalNutrientLevels[index] -= deltaNutrients[index]
				internalNutrientLevels[index] += deltaNutrients[index]
			else:
				# TODO -- this should maintain mass... 
				internalNutrientLevels[index] = externalNutrientLevels[index]
		

		# set internal nutrients based on availability from transport
		fba.externalMoleculeLevelsIs(internalNutrientLevels.tolist())

		# import ipdb; ipdb.set_trace()

		dryMassOld = dryMass
		bXout = fba.outputMoleculeLevelsChange() # this is flux of biomass
		
		# when bXout approaches 0, it becomes nan. 
		if bXout == bXout and bXout != 0:
			dryMass = dryMassOld * np.exp(bXout * (timeStep/360))
			internalNutrientLevels += fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld) #* (timeStep/360) # TODO -- is the timestep here correct?

		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = internalNutrientLevels.tolist()
		externalVec[i] = externalNutrientLevels.tolist()


	#disable kinetic targets once simulation is done.
	fba.disableKineticTargets()

	plt.figure(10)
	plt.plot(time, objVec)
	plt.xlabel('time (s)')
	plt.title('objective')
	plt.savefig(os.path.join(directory, 'process_objective.png'))

	plt.figure(11)
	plt.plot(time, bioMassVec)
	plt.xlabel('time (s)')
	plt.title('biomass')
	plt.savefig(os.path.join(directory, 'process_biomass.png'))

	plt.figure(12)
	for n in range(len(nutrientLevels)):
		plt.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	plt.xlabel('time (s)')
	plt.ylabel('concentration')
	plt.legend()
	plt.title('internal nutrient levels')
	plt.savefig(os.path.join(directory, 'process_internalNutrients.png'))


	plt.figure(13)
	for n in range(len(externalNutrientLevels)):
		plt.plot(time, externalVec[:, n], label=externalMolecules[n])
	plt.xlabel('time (s)')
	plt.ylabel('concentration')
	plt.legend()
	plt.title('external nutrient levels')
	plt.savefig(os.path.join(directory, 'process_externalNutrients.png'))






