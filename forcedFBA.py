
"""
Forced Transport Flux Balance Analysis for simplified core network

"""

import numpy as np
import scipy.constants
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from wholecell.utils.modular_fba import FluxBalanceAnalysis

#FBA types
NO_TRANSPORT = False
TRANSPORT_CONSTRAINT = True # transporter fluxes are set as constraints on FBA optimization
TRANSPORT_OBJECTIVE = False # transporter fluxes are set as objectives with a penalty
TRANSPORT_PROCESS = False  # transporter fluxes set internal concentrations, which are pulled from by FBA

#update conditions
UPDATE_NUTRIENTS = True
UPDATE_TRANSPORTERS = False

transportSD = 0.1 #standard deviation for gaussian transport update

#michaelis-menten function, gives flux of reaction
def michaelisMenten(Kcat, Km, substrate, enzyme):
	flux = enzyme * substrate * Kcat / (Km + substrate)
	return flux

directory = os.path.join('user','forcedFlux')

timeTotal = 300.0 #93.0 #seconds
timeStep = 0.05  #seconds
time = np.linspace(timeStep, timeTotal, num=int(timeTotal/timeStep))

kineticObjectiveWeight = 1e-07

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
forcedMolecules = ['Carbon1', 'Carbon2', 'D', 'E', 'Fext', 'Hext', 'Oxygen'] # these are the forced molecules, including some that are internal
externalExchangedMolecules = ['Carbon1', 'Carbon2', 'Dext', 'Eext', 'Fext', 'Hext', 'Oxygen']

objectiveFunction = {'Biomass': 1} #might be -1

transportKinetics = {
	'Tc1': {'Kcat': 0.5, 'Km': 0.2, 'active': True},
	'Tc2': {'Kcat': 1, 'Km': 10, 'active': True},
	'Tf': {'Kcat': 0.2, 'Km': 0.1, 'active': False},
	'Td': {'Kcat': 10.0, 'Km': 2.5, 'active': False},  # this reaction expells material
	'Te': {'Kcat': 10.0, 'Km': 2.5, 'active': False},  # this reaction expells material
	'Th': {'Kcat': 10.0, 'Km': 2.5, 'active': False},
	'To2': {'Kcat': 10.0, 'Km': 2.5, 'active': False},
	}

#initial state
cellDensity = 1100.0 #g/L
initdryMass = 400.0 #fg = 1e-15 g
envVolume = 10.0 #L
initnutrientLevels = np.full((len(externalExchangedMolecules)), 25) #+ np.random.normal(0,2,len(externalExchangedMolecules)) # add gaussian noise
inittransportLevels = np.full((len(transportID)), 0.5)









#transport kinetics are not forced
if NO_TRANSPORT:
	#initialize fba
	fbaObjectOptions = {
		'reactionStoich' : reactionStoich,
		'externalExchangedMolecules' : externalExchangedMolecules,
		'objective' : objectiveFunction,
		}
	fba = FluxBalanceAnalysis(**fbaObjectOptions)

	#initialize states
	dryMass = initdryMass
	nutrientLevels = initnutrientLevels

	fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

	#initialize vectors for recording data
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
			if UPDATE_NUTRIENTS:
				nutrientLevels = nutrientLevels + fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld)

		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = nutrientLevels.tolist()

		# assert sum(dryMass-dryMassOld) == sum(nutrientLevels-nutrientLevelsOld) # these should be equal.. right?

	#plot
	f = plt.figure(figsize=(18,10))

	ax1 = f.add_subplot(2, 2, 1)
	ax1.plot(time, objVec)
	ax1.set(xlabel='time (s)')
	ax1.set_title('objective')

	ax2 = f.add_subplot(2, 2, 2)
	ax2.plot(time, bioMassVec)
	ax2.set(xlabel='time (s)')
	ax2.set_title('biomass')

	ax3 = f.add_subplot(2, 2, 3)
	for n in range(len(nutrientLevels)):
		ax3.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	ax3.set(xlabel='time (s)', ylabel='concentration')
	ax3.set_title('external nutrient levels')
	ax3.legend(prop={'size': 6})

	f.subplots_adjust(wspace=0.4, hspace=0.5) #wspace = width between subplots, hspace = height between subplots
	f.savefig(os.path.join(directory, 'noForced_transport.png'), dpi=200)








#transport kinetics are set as constraints
if TRANSPORT_CONSTRAINT:
	#initialize fba
	fbaObjectOptions = {
		'reactionStoich' : reactionStoich,
		'externalExchangedMolecules' : externalExchangedMolecules,
		'objective' : objectiveFunction,
		}
	fba = FluxBalanceAnalysis(**fbaObjectOptions)

	dryMass = initdryMass
	nutrientLevels = initnutrientLevels
	transportLevels = inittransportLevels

	fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

	objVec = np.empty(len(time))
	bioMassVec = np.empty(len(time))
	nutrientVec = np.empty([len(time),len(externalExchangedMolecules)])
	transportVec = np.empty([len(time),len(transportLevels)])
	fluxVec = np.empty([len(time),len(externalExchangedMolecules)])	
	# transportVec = np.empty([len(time),len(transportID)])
	# fluxVec = np.empty([len(time),len(transportID)])

	# iterate FBA over time, updating mass and environment. fba.outputMoleculeLevelsChange() MAKES FBA UPDATE!
	for i, t in enumerate(time):
		fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

		#update transporter levels
		if UPDATE_TRANSPORTERS:
			transportLevels += np.random.normal(0,transportSD,len(transportLevels)) * timeStep
			transportLevels[transportLevels<0] = 0

		#compute transport flux
		# transportDict = {ID: None for ID in transportID}
		for ID, kinetics in transportKinetics.iteritems():
			if kinetics['active'] is True:
				index=transportID.index(ID)
				Kcat = kinetics['Kcat']
				Km = kinetics['Km']

				#get this transporter's substrate (value = -1) from reactionStoic
				substrateName = reactionStoich[ID].keys()[reactionStoich[ID].values().index(-1)]
				nutrientIndex = forcedMolecules.index(substrateName)
				substrateConc = nutrientLevels[nutrientIndex]
				transportConc = transportLevels[index]

				transportFlux = michaelisMenten(Kcat, Km, substrateConc, transportConc) #transporter concentration set to 1
				
				#constrain flux
				if substrateConc>0.1:
					# use transport levels as constraint
					fba.minReactionFluxIs(ID, transportFlux)
					fba.maxReactionFluxIs(ID, transportFlux)
				else:
					#cut nutrients off, to avoid FBA errors
					nutrientLevels[nutrientIndex]=0.0
					fba.minReactionFluxIs(ID, 0)
					fba.maxReactionFluxIs(ID, 0)

		dryMassOld = dryMass
		# bXout = fba.outputMoleculeLevelsChange() # this is flux of biomass
		# bXout = fba.biomassReactionFlux()
		fba.outputMoleculeLevelsChange() 
		growthIndex = np.where(fba.reactionIDs() == 'Growth')[0][0]
		bXout = fba.reactionFluxes()[growthIndex]

		#nutrient concentrations should change based on volume
		#nutrient levels are mmol/L, 
		# envVolume
		# import ipdb; ipdb.set_trace()
		
		if bXout == bXout and bXout != 0: # when bXout approaches 0, it becomes nan. 
			dryMass = dryMassOld * np.exp(bXout * (timeStep/360))
			if UPDATE_NUTRIENTS:
				nutrientLevels += fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld) / envVolume

		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = nutrientLevels.tolist()
		transportVec[i] = transportLevels
		fluxVec[i] = fba.externalExchangeFluxes()

	#plot
	f = plt.figure(figsize=(18,10))

	ax1 = f.add_subplot(2, 3, 1)
	ax1.plot(time, bioMassVec)
	ax1.set(xlabel='time (s)', ylabel='DCW (fg)')
	ax1.set_title('biomass')

	ax2 = f.add_subplot(2, 3, 2, projection='3d')
	substrate=np.linspace(0.0, 50.0, num=21)
	enzyme=np.linspace(0.0, 1.5, num=21)
	substrate, enzyme = np.meshgrid(substrate, enzyme)
	for ID, kinetics in transportKinetics.iteritems():
		if kinetics['active'] is True:
		# if ID is 'Tc2':
			nflux=np.empty_like(substrate)
			for ind, sConc in enumerate(substrate):
				Kcat = kinetics['Kcat']
				Km = kinetics['Km']
				eConc = enzyme[ind]
				nflux[ind]=michaelisMenten(Kcat, Km, sConc, eConc)
			ax2.plot_wireframe(substrate, enzyme, nflux)#, label=ID)
			#put simulated behavior on top of michaelis menten
			ind = [i for i, x in enumerate(transportID) if x == ID] #transportID[n]
			ax2.scatter(nutrientVec[:, ind], transportVec[:, ind], abs(fluxVec[:, ind]), '.', label=ID)
	ax2.set(xlabel='substrate concentration (mMol)', zlabel='flux', ylabel='enzyme concentration (mMol)')
	ax2.set_title('transport observed range')

	ax3 = f.add_subplot(2, 3, 3)
	for ID, kinetics in transportKinetics.iteritems():
		if kinetics['active'] is True:
			nflux=np.empty_like(time)
			ind = [i for i, x in enumerate(transportID) if x == ID] #transportID[n]
			for t in range(len(nflux)):
				Kcat = kinetics['Kcat']
				Km = kinetics['Km']
				sConc = nutrientVec[t, ind]
				eConc = transportVec[t, ind]
				nflux[t] = -1 * michaelisMenten(Kcat, Km, sConc, eConc)
			ax3.plot(time, nflux, label=ID+' m-m')
			ax3.plot(time, fluxVec[:, ind], '--', linewidth=2, label=ID+' sim')
	ax3.set(xlabel='time (s)')
	ax3.set_title('transport flux')
	ax3.legend(prop={'size': 8}, bbox_to_anchor=(1.3, 1.0))


	ax4 = f.add_subplot(2, 3, 4)
	for n in range(len(transportID)):
		ax4.plot(time, transportVec[:, n], label=transportID[n])
	ax4.set(xlabel='time (s)', ylabel='transporter conc (mMol)')
	ax4.set_title('transporter levels')
	ax4.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))

	ax5 = f.add_subplot(2, 3, 5)
	for n in range(len(nutrientLevels)):
		ax5.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	ax5.set(xlabel='time (s)', ylabel='concentration (mMol)')
	ax5.set_title('external nutrient levels')
	ax5.legend(prop={'size': 8}, bbox_to_anchor=(1.3, 1.0))

	ax6 = f.add_subplot(2, 3, 6)
	for n in range(len(transportID)):
		ax6.plot(time, fluxVec[:, n], label=transportID[n])
	ax6.set(xlabel='time (s)')
	ax6.set_title('transport flux')
	ax6.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))
	ax6.set_ylim(-5, 5)

	f.subplots_adjust(wspace=0.4, hspace=0.3) #wspace = width between subplots, hspace = height between subplots
	f.savefig(os.path.join(directory, 'constraint_transport.png'), dpi=600)








#transport kinetics are set as objectives
if TRANSPORT_OBJECTIVE:
	#initialize fba
	externalExchangedMolecules = ['Carbon1', 'Carbon2', 'Dext', 'Eext', 'Fext', 'Hext', 'Oxygen', 'Biomass']

	fbaObjectOptions = {
		'reactionStoich' : reactionStoich,
		'externalExchangedMolecules' : externalExchangedMolecules,
		'objective' : objectiveFunction,
		'objectiveType' : 'kinetic_only',
		'objectiveParameters' : {
			'kineticObjectiveWeight' : kineticObjectiveWeight,
			'reactionRateTargets' : {reaction : 1 for reaction, kinetics in transportKinetics.iteritems() if kinetics['active'] is True},
			'oneSidedReactionTargets' : [],
			}
		}
	fba = FluxBalanceAnalysis(**fbaObjectOptions)
	fba.enableKineticTargets()

	dryMass = initdryMass
	nutrientLevels = np.full((len(externalExchangedMolecules)), 25)
	# nutrientLevels = initnutrientLevels
	transportLevels = inittransportLevels

	fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

	objVec = np.empty(len(time))
	bioMassVec = np.empty(len(time))
	nutrientVec = np.empty([len(time),len(externalExchangedMolecules)])
	transportVec = np.empty([len(time),len(transportID)])
	fluxVec = np.empty([len(time),len(externalExchangedMolecules)])


	# iterate FBA over time, updating mass and environment. fba.outputMoleculeLevelsChange() MAKES FBA UPDATE!
	for i, t in enumerate(time):
		fba.externalMoleculeLevelsIs(nutrientLevels.tolist())

		#update transporter levels
		if UPDATE_TRANSPORTERS:
			transportLevels += np.random.normal(0,transportSD,len(transportLevels)) * timeStep
			transportLevels[transportLevels<0] = 0

		#compute transport flux
		# transportDict = {ID: None for ID in transportID}
		IDs=[]
		targets=[]
		for ID, kinetics in transportKinetics.iteritems():
			if kinetics['active'] is True:
				index=transportID.index(ID)
				Kcat = kinetics['Kcat']
				Km = kinetics['Km']

				#get this transporter's substrate (value = -1) from reactionStoic
				substrateName = reactionStoich[ID].keys()[reactionStoich[ID].values().index(-1)]
				nutrientIndex = forcedMolecules.index(substrateName)
				substrateConc = nutrientLevels[nutrientIndex]
				transportConc = transportLevels[index]

				transportFlux = michaelisMenten(Kcat, Km, substrateConc, transportConc) #transporter concentration set to 1

				#update flux targets
				if substrateConc>0.1:
					IDs.append(ID)
					targets.append(transportFlux)
				else:
					#cut nutrients off, to avoid FBA errors
					nutrientLevels[nutrientIndex]=0.0
					IDs.append(ID)
					targets.append(0.0)
					# import ipdb; ipdb.set_trace()

		# set kinetic flux targets according to determined transport fluxes]
		fba.setKineticTarget(IDs, targets)

		dryMassOld = dryMass
		# bXout = fba.outputMoleculeLevelsChange() # this is flux of biomass
		fba.outputMoleculeLevelsChange() 
		growthIndex = np.where(fba.reactionIDs() == 'Growth')[0][0]
		bXout = fba.reactionFluxes()[growthIndex]

		# import ipdb; ipdb.set_trace()
		# TODO -- write new FBA objective, which includes biomass maximization and kinetic targets.

		if bXout == bXout and bXout != 0: # when bXout approaches 0, it becomes nan. 
			dryMass = dryMassOld * np.exp(bXout * (timeStep/360))
			if UPDATE_NUTRIENTS:
				nutrientLevels = nutrientLevels + fba.externalExchangeFluxes()/bXout * (dryMass - dryMassOld) 

		#save starting values
		objVec[i] = fba.objectiveValue()
		bioMassVec[i] = dryMass
		nutrientVec[i] = nutrientLevels.tolist()
		transportVec[i] = transportLevels
		fluxVec[i] = fba.externalExchangeFluxes()

	#disable kinetic targets once simulation is done.
	fba.disableKineticTargets()

	#plot
	f = plt.figure(figsize=(18,10))

	ax1 = f.add_subplot(2, 3, 1)
	ax1.plot(time, bioMassVec)
	ax1.set(xlabel='time (s)')
	ax1.set_title('biomass')

	ax2 = f.add_subplot(2, 3, 2, projection='3d')
	substrate=np.linspace(0.0, 50.0, num=21)
	enzyme=np.linspace(0.0, 1.5, num=21)
	substrate, enzyme = np.meshgrid(substrate, enzyme)
	for ID, kinetics in transportKinetics.iteritems():
		if kinetics['active'] is True:
		# if ID is 'Tc2':
			nflux=np.empty_like(substrate)
			for ind, sConc in enumerate(substrate):
				Kcat = kinetics['Kcat']
				Km = kinetics['Km']
				eConc = enzyme[ind]
				nflux[ind]=michaelisMenten(Kcat, Km, sConc, eConc)
			ax2.plot_wireframe(substrate, enzyme, nflux)#, label=ID)
			#put simulated behavior on top of michaelis menten
			ind = [i for i, x in enumerate(transportID) if x == ID] #transportID[n]
			ax2.scatter(nutrientVec[:, ind], transportVec[:, ind], abs(fluxVec[:, ind]), '.', label=ID)
	ax2.set(xlabel='substrate concentration', zlabel='flux', ylabel='enzyme concentration')
	ax2.set_title('transport observed range')

	ax3 = f.add_subplot(2, 3, 3)
	for ID, kinetics in transportKinetics.iteritems():
		if kinetics['active'] is True:
			nflux=np.empty_like(time)
			ind = [i for i, x in enumerate(transportID) if x == ID] #transportID[n]
			for t in range(len(nflux)):
				Kcat = kinetics['Kcat']
				Km = kinetics['Km']
				sConc = nutrientVec[t, ind]
				eConc = transportVec[t, ind]
				nflux[t] = -1 * michaelisMenten(Kcat, Km, sConc, eConc)
			ax3.plot(time, nflux, label=ID+' m-m')
			ax3.plot(time, fluxVec[:, ind], '--', linewidth=2, label=ID+' sim')
	ax3.set(xlabel='time (s)')
	ax3.set_title('transport flux')
	ax3.legend(prop={'size': 8}, bbox_to_anchor=(1.3, 1.0))


	ax4 = f.add_subplot(2, 3, 4)
	for n in range(len(transportID)):
		ax4.plot(time, transportVec[:, n], label=transportID[n])
	ax4.set(xlabel='time (s)')
	ax4.set_title('transporter levels')
	ax4.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))

	ax5 = f.add_subplot(2, 3, 5)
	for n in range(len(nutrientLevels)):
		ax5.plot(time, nutrientVec[:, n], label=externalExchangedMolecules[n])
	ax5.set(xlabel='time (s)', ylabel='concentration')
	ax5.set_title('external nutrient levels')
	ax5.legend(prop={'size': 8}, bbox_to_anchor=(1.3, 1.0))

	ax6 = f.add_subplot(2, 3, 6)
	for n in range(len(transportID)):
		ax6.plot(time, fluxVec[:, n], label=transportID[n])
	ax6.set(xlabel='time (s)')
	ax6.set_title('transport flux')
	ax6.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))
	ax6.set_ylim(-5, 5)

	f.subplots_adjust(wspace=0.4, hspace=0.3) #wspace = width between subplots, hspace = height between subplots
	f.savefig(os.path.join(directory, 'objectivebased_transport.png'), dpi=200)





#transport kinetics are set as objectives
if TRANSPORT_PROCESS:
	#initialize fba
	fbaObjectOptions = {
		'reactionStoich' : reactionStoich,
		'externalExchangedMolecules' : externalExchangedMolecules,
		'objective' : objectiveFunction,
		}
	fba = FluxBalanceAnalysis(**fbaObjectOptions)

	#initialize stores for added step
	externalMolecules = ['Carbon1_out', 'Carbon2_out', 'Dext_out', 'Eext_out', 'Fext_out', 'Hext_out', 'Oxygen_out']
	dryMass = initdryMass
	externalNutrientLevels = initnutrientLevels
	internalNutrientLevels = np.empty_like(initnutrientLevels)
	transportLevels = inittransportLevels

	objVec = np.empty(len(time))
	bioMassVec = np.empty(len(time))
	nutrientVec = np.empty([len(time),len(externalExchangedMolecules)])
	externalVec = np.empty([len(time),len(externalExchangedMolecules)])
	transportVec = np.empty([len(time),len(transportID)])
	fluxVec = np.empty([len(time),len(transportID)])

	# iterate FBA over time, updating mass and environment. fba.outputMoleculeLevelsChange() MAKES FBA UPDATE!
	for i, t in enumerate(time):

		#update transporter levels
		transportLevels += np.random.normal(0,0.5,len(transportLevels)) * timeStep
		transportLevels[transportLevels<0] = 0

		#compute transport flux
		transportDict = {ID: None for ID in transportID}
		for ID, kinetics in transportKinetics.iteritems():
			if kinetics['active'] is True:

				index=transportID.index(ID)
				Kcat = kinetics['Kcat']
				Km = kinetics['Km']

				#get this transporter's substrate (value = -1) from reactionStoic
				substrateName = reactionStoich[ID].keys()[reactionStoich[ID].values().index(-1)]
				nutrientIndex = forcedMolecules.index(substrateName)
				substrateConc = nutrientLevels[nutrientIndex]
				transportConc = transportLevels[index]
				transportDict[ID] = michaelisMenten(Kcat, Km, substrateConc, transportConc) #transporter concentration set to 1

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
		transportVec[i] = transportLevels
		fluxVec[i] = fba.externalExchangeFluxes()








	#plot
	f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(18,10))
	f.subplots_adjust(wspace=0.4, hspace=0.3) #wspace = width between subplots, hspace = height between subplots
	ax1.plot(time, bioMassVec)
	ax1.set(xlabel='time (s)')
	ax1.set_title('biomass')

	substrate=np.linspace(0.0, 50.0, num=1001)
	for ID, kinetics in transportKinetics.iteritems():
		nflux=np.empty_like(substrate)
		for ind, sConc in enumerate(substrate):
			Kcat = kinetics['Kcat']
			Km = kinetics['Km']
			nflux[ind]=michaelisMenten(Kcat, Km, sConc, 1)
		ax2.plot(substrate, nflux, label=ID)
	ax2.set(xlabel='substrate concentration', ylabel='flux')
	ax2.set_title('transport Michaelis-Menten')
	ax2.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))

	substrate=np.linspace(0.0, 50.0, num=1001)
	for ID, kinetics in transportKinetics.iteritems():
		nflux=np.empty_like(substrate)
		for ind, sConc in enumerate(substrate):
			Kcat = kinetics['Kcat']
			Km = kinetics['Km']
			nflux[ind]=michaelisMenten(Kcat, Km, sConc, 1)
		ax3.plot(substrate, nflux)#, label=ID)
	#put actual behavior on top of michaelis menten
	for n in transportConstraint:
		ind = [i for i, x in enumerate(transportID) if x == n] #transportID[n]
		ax3.plot(nutrientVec[:, ind], abs(fluxVec[:, ind]), '.', label=n)
	ax3.set(xlabel='substrate concentration', ylabel='flux')
	ax3.set_title('transport observed range')
	ax3.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))

	for n in range(len(transportID)):
		ax4.plot(time, transportVec[:, n], label=transportID[n])
	ax4.set(xlabel='time (s)')
	ax4.set_title('transporter levels')
	ax4.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))

	for n in range(len(externalNutrientLevels)):
		ax5.plot(time, externalVec[:, n], label=externalExchangedMolecules[n])
	ax5.set(xlabel='time (s)', ylabel='concentration')
	ax5.set_title('external nutrient levels')
	ax5.legend(prop={'size': 8}, bbox_to_anchor=(1.3, 1.0))


	for n in range(len(transportID)):
		ax6.plot(time, fluxVec[:, n], label=transportID[n])
	ax6.set(xlabel='time (s)')
	ax6.set_title('transport flux')
	ax6.legend(prop={'size': 8}, bbox_to_anchor=(1.2, 1.0))

	plt.savefig(os.path.join(directory, 'process_based.png'))






