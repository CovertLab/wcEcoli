
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle
import scipy.stats

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

NUMERICAL_ZERO = 1e-12

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):
	if metadata["variant"] != "tfActivity":
		print "This plot only runs for the 'tfActivity' variant."
		return

	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	ap = AnalysisPaths(inputDir, variant_plot = True)
	variants = sorted(ap._path_data['variant'].tolist()) # Sorry for accessing private data

	if 0 in variants:
		variants.remove(0)

	if len(variants) == 0:
		return

	all_cells = sorted(ap.get_cells(variant = variants, seed = [0], generation = [0]))

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	expectedProbBound = []
	simulatedProbBound = []
	expectedSynthProb = []
	simulatedSynthProb = []
	targetId = []
	targetCondition = []

	for variant, simDir in zip(variants, all_cells):
		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))
		tfList = ["basal (no TF)"] + sorted(sim_data.tfToActiveInactiveConds)
		simOutDir = os.path.join(simDir, "simOut")
		tf = tfList[(variant + 1) // 2]
		tfStatus = None
		if variant % 2 == 1:
			tfStatus = "active"
		else:
			tfStatus = "inactive"

		bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")

		rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		rnaIds = rnaSynthProbReader.readAttribute("rnaIds")


		for tfTarget in sorted(sim_data.tfToFC[tf]):

			tfTargetBoundId = [tfTarget + "__" + tf]
			tfTargetBoundIndex = np.array([bulkMoleculeIds.index(x) for x in tfTargetBoundId])
			tfTargetBoundCounts = bulkMoleculesReader.readColumn("counts")[:, tfTargetBoundIndex].reshape(-1)

			expectedProbBound.append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
			simulatedProbBound.append(tfTargetBoundCounts[5:].mean())


			tfTargetSynthProbId = [tfTarget + "[c]"]
			tfTargetSynthProbIndex = np.array([rnaIds.index(x) for x in tfTargetSynthProbId])
			tfTargetSynthProb = rnaSynthProbReader.readColumn("rnaSynthProb")[:, tfTargetSynthProbIndex].reshape(-1)

			rnaIdx = np.where(sim_data.process.transcription.rnaData["id"] == tfTarget + "[c]")[0][0]

			expectedSynthProb.append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
			simulatedSynthProb.append(tfTargetSynthProb[5:].mean())

			targetId.append(tfTarget)
			targetCondition.append(tf + "__" + tfStatus)

		bulkMoleculesReader.close()
		rnaSynthProbReader.close()

	expectedProbBound = np.array(expectedProbBound)
	simulatedProbBound = np.array(simulatedProbBound)
	expectedSynthProb = np.array(expectedSynthProb)
	simulatedSynthProb = np.array(simulatedSynthProb)

	regressionResult = scipy.stats.linregress(np.log10(expectedProbBound[expectedProbBound > NUMERICAL_ZERO]), np.log10(simulatedProbBound[expectedProbBound > NUMERICAL_ZERO]))
	regressionResultLargeValues = scipy.stats.linregress(np.log10(expectedProbBound[expectedProbBound > 1e-2]), np.log10(simulatedProbBound[expectedProbBound > 1e-2]))

	ax = plt.subplot(2, 1, 1)
	ax.scatter(np.log10(expectedProbBound), np.log10(simulatedProbBound))
	plt.xlabel("log10(Expected probability bound)", fontsize = 6)
	plt.ylabel("log10(Simulated probability bound)", fontsize = 6)
	plt.title("Slope: %0.3f   Intercept: %0.3e      (Without Small Values:  Slope: %0.3f Intercept: %0.3e)" % (regressionResult.slope, regressionResult.intercept, regressionResultLargeValues.slope, regressionResultLargeValues.intercept), fontsize = 6)
	ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

	regressionResult = scipy.stats.linregress(np.log10(expectedSynthProb[expectedSynthProb > NUMERICAL_ZERO]), np.log10(simulatedSynthProb[expectedSynthProb > NUMERICAL_ZERO]))

	ax = plt.subplot(2, 1, 2)
	ax.scatter(np.log10(expectedSynthProb), np.log10(simulatedSynthProb))
	plt.xlabel("log10(Expected synthesis probability)", fontsize = 6)
	plt.ylabel("log10(Simulated synthesis probability)", fontsize = 6)
	plt.title("Slope: %0.3f   Intercept: %0.3e" % (regressionResult.slope, regressionResult.intercept), fontsize = 6)
	ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

	plt.tight_layout()

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


	# Bokeh
	from bokeh.plotting import figure, output_file, ColumnDataSource, show
	from bokeh.models import HoverTool, BoxZoomTool, LassoSelectTool, PanTool, WheelZoomTool, ResizeTool, UndoTool, RedoTool
	from bokeh.io import vplot

	source1 = ColumnDataSource(data = dict(x = np.log10(expectedProbBound), y = np.log10(simulatedProbBound), ID = targetId, condition = targetCondition))
	hover1 = HoverTool(tooltips = [("ID", "@ID"), ("condition", "@condition")])
	tools1 = [hover1, BoxZoomTool(), LassoSelectTool(), PanTool(), WheelZoomTool(), ResizeTool(),	UndoTool(),	RedoTool(), "reset"]
	s1 = figure(
		x_axis_label = "log10(Expected probability bound)",
		y_axis_label = "log10(Simulated probability bound)",
		width = 800,
		height = 400,
		tools = tools1)
	s1.scatter("x", "y", source = source1)

	if not os.path.exists(os.path.join(plotOutDir, "html_plots")):
		os.makedirs(os.path.join(plotOutDir, "html_plots"))
	import bokeh.io
	bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + "__probBound" + ".html"), title = plotOutFileName, autosave = False)
	bokeh.io.save(s1)


	source2 = ColumnDataSource(data = dict(x = np.log10(expectedSynthProb), y = np.log10(simulatedSynthProb), ID = targetId, condition = targetCondition))
	hover2 = HoverTool(tooltips = [("ID", "@ID"), ("condition", "@condition")])
	tools2 = [hover2, BoxZoomTool(), LassoSelectTool(), PanTool(), WheelZoomTool(), ResizeTool(),	UndoTool(),	RedoTool(), "reset"]
	s2 = figure(
		x_axis_label = "log10(Expected synthesis probability)",
		y_axis_label = "log10(Simulated synthesis probability)",
		width = 800,
		height = 400,
		tools = tools2)
	s2.scatter("x", "y", source = source2)

	bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + "__synthProb" + ".html"), title = plotOutFileName, autosave = False)
	bokeh.io.save(s2)


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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
