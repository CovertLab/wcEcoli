"""
Complete-lineage variant of subgen_expression_definitions.py.

Runs the identical five-definition subgenerational-expression analysis but only
on seeds whose lineage successfully reached the final generation. Lineages that
died mid-run are a major source of "near-ubiquitous" absences (a dying cell drops
many genes at once in its last generations), so excluding them tests how much of
the apparent subgenerational signal is a dying-lineage artifact versus genuine
expression heterogeneity in healthy, dividing cells.

All outputs, columns, and behavior match subgen_expression_definitions.py; only
the cell set differs. Compare the two runs' _definition_summary.tsv to see how the
subgenerational gene counts shift once incomplete lineages are removed.
"""

from models.ecoli.analysis.cohort import subgen_expression_definitions


class Plot(subgen_expression_definitions.Plot):
	REQUIRE_COMPLETE_LINEAGE = True


if __name__ == '__main__':
	Plot().cli()
