"""
Common code for Cohort analysis plots.

TODO: A method to instantiate and run a list of subclasses in a controlled
order, with unneeded ones commented out, to simplify the Firetask.
"""

from __future__ import absolute_import
from __future__ import division

from models.ecoli.analysis import analysisPlot


class CohortAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for Cohort analysis plots."""
	pass
