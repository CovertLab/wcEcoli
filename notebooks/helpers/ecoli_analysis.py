import os
import importlib
from functools import partial
from IPython.display import SVG

class analysisPlot:
    def __init__(self, wcEcoliPath, simOutPath, plotOutPath, simDataPath, validationDataPath):
        self._wcEcoliPath = wcEcoliPath
        self._simOutPath = simOutPath
        self._plotOutPath = plotOutPath
        self._simDataPath = simDataPath
        self._validationDataPath = validationDataPath
        
        files = os.listdir(os.path.join(wcEcoliPath, "models", "ecoli", "analysis", "single"))
        analysisScripts = [x for x in files if x.endswith(".py")]
        
        for analysisScript in analysisScripts:
            plot = analysisScript[:-3]
            setattr(self, plot, partial(self.analyse, plot))
            
    def analyse(self, plot):
        mod = importlib.import_module("models.ecoli.analysis.single." + plot)
        mod.main(self._simOutPath, self._plotOutPath, plot, self._simDataPath, self._validationDataPath, None)
        return SVG(filename=os.path.join(PLOT_OUT_PATH, "svg_plots", "%s.svg" % plot))
        
analysis = analysisPlot(WC_ECOLI_DIRECTORY, SIM_OUT_PATH, PLOT_OUT_PATH, SIM_DATA_PATH, None)

import xml.etree.ElementTree as ET
import requests

def showReactionPathway(reaction):
    resp = requests.get("http://websvc.biocyc.org/getxml?id=ECOLI:%s&detail=full" % reaction)
    root = ET.fromstring(resp.text)

    pathways = root.findall(".//Pathway")
    if len(pathways) == 0:
        print "No pathway found for %s" % reaction
        
    for pathway in pathways:
        pathID = pathway.attrib['frameid']
        display(Image(url="http://biocyc.org/ECOLI/diagram-only?type=PATHWAY&object=%s" % pathID))
        print "Pathway:", "http://ecocyc.org/ECOLI/NEW-IMAGE?type=PATHWAY&object=%s" % pathID        
        
    print "Reaction:", "http://ecocyc.org/ECOLI/NEW-IMAGE?type=NIL&object=%s&redirect=T" % reaction