#!/usr/bin/env python
"""
Plot reaction max rate over course of the simulation.

@date: Created 1/08/2016
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
import pandas as pd
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import mpld3
from mpld3 import plugins, utils

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

class HighlightLines(plugins.PluginBase):
    """A plugin to highlight lines on hover"""

    JAVASCRIPT = """
    mpld3.register_plugin("linehighlight", LineHighlightPlugin);
    LineHighlightPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    LineHighlightPlugin.prototype.constructor = LineHighlightPlugin;
    LineHighlightPlugin.prototype.requiredProps = ["line_ids", "line_names"];
    LineHighlightPlugin.prototype.defaultProps = {alpha_bg:0.3, alpha_fg:1.0}
    function LineHighlightPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    LineHighlightPlugin.prototype.draw = function(){
      for(var i=0; i<this.props.line_ids.length; i++){
         var obj = mpld3.get_element(this.props.line_ids[i], this.fig),
             alpha_fg = this.props.alpha_fg;
             alpha_bg = this.props.alpha_bg;
             id = this.props.line_names[i];
         obj.elements()
             .on("mouseover", function(d, i){
                            console.log(this);
                            d3.select(this).transition().duration(50)
                              .style("stroke-opacity", alpha_fg);
                              })
             .on("mouseout", function(d, i){
                            d3.select(this).transition().duration(200)
                              .style("stroke-opacity", alpha_bg); });
      }
    };
    """

    def __init__(self, lines, lineLabels):
        self.lines = lines
        self.dict_ = {"type": "linehighlight",
                      "line_ids": [utils.get_id(line) for line in lines],
                      "alpha_bg": lines[0].get_alpha(),
                      "alpha_fg": 1.0,
                      "line_names": [str(label) for label in lineLabels]}

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	
	enzymeKineticsArray = enzymeKineticsdata.readColumn("reactionConstraints")

	reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")
	
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	enzymeKineticsdata.close()

	reactionRateArray = np.transpose(enzymeKineticsArray)

	plt.figure(figsize = (8.5, 11))
	plt.title("Enzyme Kinetics")

	lineLabels = []
	timeCourses = []
	for idx, timeCourse in enumerate(reactionRateArray):
		if (np.amax(timeCourse) < np.inf) and (idx < len(reactionIDs)):
			timeCourses.append(list(timeCourse))
			lineLabels.append(reactionIDs[idx][:15])
	timeCourses = np.array(timeCourses)

	fig, ax = plt.subplots(subplot_kw={'xticks': [], 'yticks': []})
	lines = ax.plot(time/60, timeCourses.T, label= lineLabels, color='blue', lw=4, alpha=0.4)
	plt.xlabel("Time (min)")
	plt.ylabel("Reaction Rate ({counts_units}/{volume_units}.{time_units})".format(counts_units=COUNTS_UNITS.strUnit(), volume_units=VOLUME_UNITS.strUnit(), time_units=TIME_UNITS.strUnit()))
	plt.legend(lineLabels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plugins.connect(fig, HighlightLines(lines, lineLabels))

	# # Define some CSS to control our custom labels
	# css = """
	# table
	# {
	#   border-collapse: collapse;
	# }
	# th
	# {
	#   color: #ffffff;
	#   background-color: #000000;
	# }
	# td
	# {
	#   background-color: #cccccc;
	# }
	# table, th, td
	# {
	#   font-family:Arial, Helvetica, sans-serif;
	#   border: 1px solid black;
	#   text-align: right;
	# }
	# """

	# fig, ax = plt.subplots()
	# ax.grid(True, alpha=0.3)

	# N = 50
	# df = pd.DataFrame(index=range(N))
	# df['x'] = np.random.randn(N)
	# df['y'] = np.random.randn(N)
	# df['z'] = np.random.randn(N)

	# import ipdb; ipdb.set_trace()

	# labels = []
	# for i in range(N):
	#     label = df.ix[[i], :].T
	#     label.columns = ['Row {0}'.format(i)]
	#     # .to_html() is unicode; so make leading 'u' go away with str()
	#     labels.append(str(label.to_html()))

	# points = ax.plot(df.x, df.y, 'o', color='b',
	#                  mec='k', ms=15, mew=1, alpha=.6)

	# ax.set_xlabel('x')
	# ax.set_ylabel('y')
	# ax.set_title('HTML tooltips', size=20)

	# tooltip = plugins.PointHTMLTooltip(points[0], labels,
	#                                    voffset=10, hoffset=10, css=css)
	# plugins.connect(fig, tooltip)


	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

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
