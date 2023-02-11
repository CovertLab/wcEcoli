"""
Common code for Comparison analysis plots.
"""

from models.ecoli.analysis import analysisPlot
from runscripts.manual import analysisComparison


# noinspection PyAbstractClass
class ComparisonAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for Comparison analysis plots."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = analysisComparison.AnalysisComparison(analysis_plotter=self)
		script.cli()
