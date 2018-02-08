#!/usr/bin/env python
"""
Track mass flow following metabolism via sankey diagram

@author: Nicole Ferraro, code included from other analysis scripts
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/11/17
"""

import argparse
import cPickle
import json
import os
import shutil

import numpy as np

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

    if not os.path.isdir(simOutDir):
        raise Exception, "simOutDir does not currently exist as a directory"

    if not os.path.exists(plotOutDir):
        os.mkdir(plotOutDir)

    #Extract metabolites produced at each time step
    fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
    time = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
    metaboliteIds = fbaResults.readAttribute("metaboliteNames")
    externalIds = fbaResults.readAttribute("externalMoleculeIDs")
    deltaMetabolites = fbaResults.readColumn("deltaMetabolites")
    #Extract molecules exchanged with external environment at each time step
    exFluxes = fbaResults.readColumn("externalExchangeFluxes")

    fbaResults.close()

    #Obtain levels of building blocks, includes counts of amino acids and nucleotides

    bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
    growth = TableReader(os.path.join(simOutDir, "GrowthLimits"))
    moleculeIds = bulkMolecules.readAttribute("objectNames")
    sim_data = cPickle.load(open(simDataFile))
    validation_data = cPickle.load(open(validationDataFile))

    allMWs = sim_data.state.bulkMolecules.bulkData.struct_array
    MWind = sim_data.molecular_weight_keys.index('metabolite')
    
    #DNTP counts
    dntpIDs = sim_data.moleculeGroups.dNtpIds
    dntpIndexes = np.array([moleculeIds.index(dntpId) for dntpId in dntpIDs], np.int)
    dntpIndexesMet = np.array([metaboliteIds.index(dntpId) for dntpId in dntpIDs], np.int)
    dntpMWInd = [allMWs["id"].tolist().index(dntpId) for dntpId in dntpIDs]
    dntpWeights = allMWs["mass"][dntpMWInd][:,MWind]
    dntpCounts = bulkMolecules.readColumn("counts")[:, dntpIndexes]
    #Since not stored as number of dntps used, need to calculate the change in amount at each time step
    dntpCountsDif = [dntpCounts[x] - dntpCounts[x-1] for x in range(1,len(dntpCounts))]
    dntpCountsDif.insert(0, dntpCounts[0])

    #NTP counts
    ntpIDs = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
    ntpIndexes = np.array([moleculeIds.index(ntpId) for ntpId in ntpIDs], np.int)
    ntpIndexesMet = np.array([metaboliteIds.index(ntpId) for ntpId in ntpIDs], np.int)
    ntpCounts = growth.readColumn("ntpUsed")
    ntpMWInd = [allMWs["id"].tolist().index(ntpId) for ntpId in ntpIDs]
    ntpWeights = allMWs["mass"][ntpMWInd][:,MWind]

    #Amino acid counts
    aaIDs = sim_data.moleculeGroups.aaIDs
    aaIndexes = np.array([moleculeIds.index(aaId) for aaId in aaIDs], np.int)
    aaIndexesMet = np.array([metaboliteIds.index(aaId) for aaId in aaIDs], np.int)
    aaCounts = growth.readColumn("aasUsed")
    aaMWInd = [allMWs["id"].tolist().index(aaId) for aaId in aaIDs]
    aaWeights = allMWs["mass"][aaMWInd][:,MWind]

    bulkMolecules.close()
    growth.close()

    metMWInd = [allMWs["id"].tolist().index(mid) for mid in metaboliteIds]
    metWeights = allMWs["mass"][metMWInd][:,MWind]

    #Multiply by molecular weight to obtain value in g / mol, then convert to fg
    nAvo = sim_data.constants.nAvogadro._value

    dntpMolar = (dntpCountsDif * dntpWeights * 1e15) / nAvo
    ntpMolar = (ntpCounts * ntpWeights * 1e15) / nAvo
    aaMolar = (aaCounts * aaWeights * 1e15) / nAvo
    deltaMetMolar = (deltaMetabolites * metWeights * 1e15) / nAvo

    #Third level, macromolecules, includes protein, RNA, DNA, small molecules, and other
    mass = TableReader(os.path.join(simOutDir, "Mass"))
    dryMass = mass.readColumn("dryMass")

    #Convert from g / g-DCW-hr to fg

    exchangeMasses = {} # avoid double counting duplicates in externalIds
    raw_data = KnowledgeBaseEcoli()
    for metabolite in raw_data.metabolites:
        for mID in externalIds:
            if mID.split("[")[0] == "WATER":
                exchangeMasses[mID] = 18.015 * exFluxes[:,externalIds.index(mID)] * 0.001 * dryMass * time / 3600
            if mID.split("[")[0] == metabolite["id"]:
                exchangeMasses[mID] = metabolite["mw7.2"] * exFluxes[:,externalIds.index(mID)] * 0.001 * dryMass * time / 3600

    #Get averages for each value across all time steps
    exFluxesAvg = {}
    for mkey in exchangeMasses.keys():
        exFluxesAvg[mkey] = np.mean(exchangeMasses[mkey])

    deltaMetabolitesAvg = np.mean(deltaMetMolar, axis=0)
    dntpMolarAvg = np.mean(dntpMolar, axis=0)
    ntpMolarAvg = np.mean(ntpMolar, axis=0)
    aaMolarAvg = np.mean(aaMolar, axis=0)

    #Not used, can check values at end with mass of aa -> prot and so on
    proteinMassAvg = np.mean(proteinMassDif, axis=0)
    rnaMassAvg = np.mean(rnaMassDif, axis=0)
    dnaMassAvg = np.mean(dnaMassDif, axis=0)

    #Get sources, targets, and values for diagram
    meta2aa = sum(deltaMetabolitesAvg[aaIndexesMet])
    meta2dntp = sum(deltaMetabolitesAvg[dntpIndexesMet])
    meta2ntp = sum(deltaMetabolitesAvg[ntpIndexesMet])
    meta2other = sum(deltaMetabolitesAvg) - meta2aa - meta2dntp - meta2ntp
    aa2prot = sum(aaMolarAvg)
    dntp2dna = sum(dntpMolarAvg)
    ntp2rna = sum(ntpMolarAvg)

    gluc2metab = exFluxesAvg['GLC[p]']*-1
    phos2metab = exFluxesAvg['PI[p]']*-1

    other2metab = 0
    metab2out = 0
    aa2metab = 0
    for key,value in exFluxesAvg.items():
        if key not in ['GLC[p]', 'PI[p]', 'WATER[p]'] and value < 0:
            other2metab += value*-1
        if key not in ['ETOH[p]', 'FORMATE[p]', 'ACET[p]', 'CARBON-DIOXIDE[p]', 'WATER[p]'] and value > 0:
            metab2out += value
        newkey = key.split('[')[0] + '[c]'
        if newkey in aaIDs:
            aa2metab += (-1*value)
            other2metab = other2metab - (value*-1)

    metab2etoh = exFluxesAvg['ETOH[p]']
    metab2formate = exFluxesAvg['FORMATE[p]']
    metab2acet = exFluxesAvg['ACET[p]']
    metab2co2 = exFluxesAvg['CARBON-DIOXIDE[p]']
    metab2meta = sum(deltaMetabolitesAvg)

    extra_metab = (gluc2metab + phos2metab + other2metab + aa2metab) - metab2meta - metab2out - metab2co2 - metab2etoh - metab2formate - metab2acet
    metab2outNew = metab2out + extra_metab

    values = [gluc2metab, phos2metab, other2metab, metab2meta, metab2outNew, metab2co2, meta2aa, meta2dntp, meta2ntp, meta2other, aa2prot, dntp2dna, ntp2rna]
    sources = ['Glucose', 'Phosphate', 'Other In', 'Metabolism', 'Metabolism', 'Metabolism', 'Metabolites', 'Metabolites', 'Metabolites', 'Metabolites', 'Amino Acids', 'dNTPs', 'NTPs']
    targets = ['Metabolism', 'Metabolism', 'Metabolism', 'Metabolites', 'Out', 'CO2', 'Amino Acids', 'dNTPs', 'NTPs', 'Other Molecules', 'Protein', 'DNA', 'RNA'] 
    
    #Prevent issues with zero weighted nodes in diagram, only add if non-zero
    if aa2metab > 1e-3:
        values.insert(2,aa2metab)
        sources.insert(2,'External Amino Acids')
        targets.insert(2, 'Metabolism')
    i = 7

    if metab2etoh > 1e-3:
        values.insert(i,metab2etoh)
        sources.insert(i,'Metabolism')
        targets.insert(i,'EtOH')
        i += 1
    if metab2formate > 1e-3:
        values.insert(i,metab2formate)
        sources.insert(i,'Metabolism')
        targets.insert(i,'Formate')
        i += 1
    if metab2acet > 1e-3:
        values.insert(i,metab2acet)
        sources.insert(i,'Metabolism')
        targets.insert(i,'Acetate')

    #Generate new javascript file with these variables for access by the html diagram
    #Original js file stored in analysis directory, new one will be moved to plotOutDir 
    js_data = {'sources':sources, 'targets':targets, 'values':values}
    js = open('sankey.js', 'r')
    js_out = "var diagram_data = JSON.parse('{}');".format(json.dumps(js_data))
    js_out += '\n'
    js_out += js.read()
    with open('sankey_new.js', 'w') as njs:
      njs.write(js_out)
    #Copy d3.js and sankey.html to the analysis folder (for access by new js file), and move sankey_new.js as well
    shutil.move("sankey_new.js", os.path.join(plotOutDir, "sankey.js"))
    #Need to copy instead, keep originals in analysis folder for future use
    shutil.copy2("d3.js", os.path.join(plotOutDir, "d3.js"))
    shutil.copy2("sankey.html", os.path.join(plotOutDir, "sankey.html"))

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)
	
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])

