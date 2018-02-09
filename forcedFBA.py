
"""
Forced Transport Flux Balance Analysis for simplified core network

"""

import numpy as np
import scipy.constants
import matplotlib.pyplot as plt
import os
from wholecell.utils.modular_fba import FluxBalanceAnalysis
# from wholecell.utils import units

NO_TRANSPORT = True
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
	'Tc1': {'Kcat': 0.5, 'Km': 0.25},
	'Tc2': {'Kcat': 0.5, 'Km': 0.25},
	'Tf': {'Kcat': 10.0, 'Km': 2.5},
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
initnutrientLevels = np.full((len(externalExchangedMolecules)), 100) + np.random.normal(0,1.5,len(externalExchangedMolecules)) # add gaussian noise







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

		fba.biomassReactionFlux() # is this biomass exchange flux?

		# when bXout approaches 0, it becomes nan. 
		if bXout == bXout and bXout != 0:
			dryMass = dryMassOld * np.exp(bXout * (timeStep/360))
			nutrientLevels = nutrientLevels + fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld)

		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = nutrientLevels.tolist()


		# assert sum(dryMass-dryMassOld) == sum(nutrientLevels-nutrientLevelsOld) # these should be equal.. right?
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
	plt.title('external molecule levels')
	plt.savefig(os.path.join(directory, 'noforce_externalMolecules.png'))

















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
			transportDict[ID] = michaelisMenten(Kcat, Km, substrateConc, 0.1) #transporter concentration set to 1

		import ipdb; ipdb.set_trace()

		for ID, value in transportDict.iteritems():
			if value:
				fba.minReactionFluxIs(ID, value)
				fba.maxReactionFluxIs(ID, value)

		


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

	plt.figure(1)
	plt.plot(time, objVec)
	plt.xlabel('time (s)')
	plt.title('objective')
	plt.savefig(os.path.join(directory, 'constraint_objective.png'))

	plt.figure(2)
	plt.plot(time, bioMassVec)
	plt.xlabel('time (s)')
	plt.title('biomass')
	plt.savefig(os.path.join(directory, 'constraint_biomass.png'))

	plt.figure(3)
	for n in range(len(nutrientLevels)):
		plt.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	plt.xlabel('time (s)')
	plt.ylabel('concentration')
	plt.legend()
	plt.title('external molecule levels')
	plt.savefig(os.path.join(directory, 'constraint_externalMolecules.png'))










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
	plt.title('external molecule levels')
	plt.savefig(os.path.join(directory, 'obj_externalMolecules.png'))










#transport kinetics are set as objectives
if TRANSPORT_PROCESS:

	# #adjust stoichiometry to include additional variables
	# reactionStoich2 = reactionStoich
	# reactionStoich2 = {
	# #transport processes
	# 'Tc1': {'Carbon1': -1, 'Carbon1in': 1},
	# 'Tc2': {'Carbon2': -1, 'Carbon2in': 1},
	# 'Tf': {'Fext': -1, 'Fin': 1},
	# 'Td': {'Din': -1, 'Dext': 1},
	# 'Te': {'Ein': -1, 'Eext': 1},
	# 'Th': {'Hext': -1, 'Hin': 1},
	# 'To2': {'Oxygen': -1, 'Oxygenin': 1},

	# #internal processes
	# 'Tc1_in': {'Carbon1in': -1, 'A': 1},
	# 'Tc2_in': {'Carbon2in': -1, 'A': 1},
	# 'Tf_in': {'Fin': -1, 'F': 1},
	# 'Td_in': {'D': -1, 'Din': 1},
	# 'Te_in': {'E': -1, 'Ein': 1},
	# 'Th_in': {'Hin': -1, 'H': 1},
	# 'To2_in': {'Oxygenin': -1, 'O2': 1},
	# }
	outsideMolecules = ['Carbon1_out', 'Carbon2_out', 'Dext_out', 'Eext_out', 'Fext_out', 'Hext_out', 'Oxygen_out']








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


		#transport processes to set downstream limits
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

		'''
		TODO -- determine internal concentrations, use them to set a new "externalmolecules" to draw from, which is actually internal
		Tc1 and Tc2 produce A
		Tf produces F
		To2 produces O2
		Th produces H
		
			external --> new internal molecules 
		Carbon1 --> Carbon1in
		Carbon2 --> Carbon2in
		Fext-->Fin

		"outsideMolecules" are actually variables that track how much was drawn and how much was pulled and how much left.
		'''


		# for ID, value in transportDict.iteritems():
		# 	if value:
		# 		fba.maxReactionFluxIs(ID, value)


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
	plt.title('external molecule levels')
	plt.savefig(os.path.join(directory, 'process_externalMolecules.png'))








