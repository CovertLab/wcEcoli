from __future__ import division

from models.ecoli.sim.variants.gene_knockout import geneKnockout
from models.ecoli.sim.variants.wildtype import wildtype
from models.ecoli.sim.variants.nutrientTimeSeries import nutrientTimeSeries
from models.ecoli.sim.variants.tf_activity import tfActivity
from models.ecoli.sim.variants.condition import condition
from models.ecoli.sim.variants.transcriptionInitiationShuffleParams import transcriptionInitiationShuffleParams
from models.ecoli.sim.variants.kineticTargetShuffleParams import kineticTargetShuffleParams
from models.ecoli.sim.variants.catalystShuffleParams import catalystShuffleParams
from models.ecoli.sim.variants.translationEfficienciesShuffleParams import translationEfficienciesShuffleParams
from models.ecoli.sim.variants.monomerDegRateShuffleParams import monomerDegRateShuffleParams
from models.ecoli.sim.variants.kineticCatalystShuffleParams import kineticCatalystShuffleParams
from models.ecoli.sim.variants.allShuffleParams import allShuffleParams
from models.ecoli.sim.variants.meneParams import meneParams
from models.ecoli.sim.variants.metabolism_kinetic_objective_weight import metabolism_kinetic_objective_weight

nameToFunctionMapping = {
	"geneKnockout": geneKnockout,
	"wildtype": wildtype,
	"nutrientTimeSeries": nutrientTimeSeries,
	"tfActivity": tfActivity,
	"condition": condition,
	"transcriptionInitiationShuffleParams": transcriptionInitiationShuffleParams,
	"kineticTargetShuffleParams": kineticTargetShuffleParams,
	"catalystShuffleParams": catalystShuffleParams,
	"translationEfficienciesShuffleParams": translationEfficienciesShuffleParams,
	"monomerDegRateShuffleParams": monomerDegRateShuffleParams,
	"kineticCatalystShuffleParams": kineticCatalystShuffleParams,
	"allShuffleParams": allShuffleParams,
	"meneParams": meneParams,
	"metabolism_kinetic_objective_weight": metabolism_kinetic_objective_weight,
}
