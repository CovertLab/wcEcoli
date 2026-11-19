"""
rpoC at each step from gene to protein, high GFP over no GFP.

# Figure 5G

Four bars, each the mean over surviving high-GFP lineages of the lineage's
window mean divided by the same seed's no-GFP lineage, with a t-based 95% CI:
	copies         rpoC gene copies (RnaSynthProb gene_copy_number)
	transcription  rpoC transcription initiations per second
	               (RnapData rna_init_event_per_cistron)
	translation    rpoC translation initiations per second
	               (RibosomeData ribosome_init_event_per_monomer)
	protein        rpoC monomers, free and inside RNAP (MonomerCounts)

"High GFP" is the design whose mean GFP proteome share is nearest the second
entry of FIG5_PROTEOME_SHARES. Translation is a raw rate here, as in 5D.

The supporting claim: rpoC loses far more at transcription and translation
than gene dosage alone would explain.
"""

from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

BARS = (
	('copies', 'rpoc_gene_copies'),
	('transcription', 'rpoc_transcription_per_s'),
	('translation', 'rpoc_translation_per_s'),
	('protein', 'rpoc_monomers'),
)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		high = gfp_paper_data.nearest_design(data.designs,
			C.FIG5_PROTEOME_SHARES[-1])
		ratios = [gfp_paper_data.seed_matched_ratio(data.lineages,
			int(high['variant']), column) for _, column in BARS]

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.7))
		style.ratio_bars(ax, ['$\\it{rpoC}$\n%s' % label for label, _ in BARS],
			[r['mean'] for r in ratios], [r['ci95'] for r in ratios],
			style.COLORS['rnap'])
		ax.set_ylabel('high GFP / no GFP')

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=min(r['n'] for r in ratios),
			notes='High GFP = variant %d (promoter %g, RBS %g, GFP %.3f of'
				' proteome). Seed-matched ratios over %d lineages: %s.'
				% (high['variant'], high['promoter'], high['rbs'],
				high['gfp_proteome_share'], min(r['n'] for r in ratios),
				', '.join('%s %.3f +/- %.3f' % (label, r['mean'], r['ci95'])
				for (label, _), r in zip(BARS, ratios))),
			**data.provenance([gfp_paper_data.CONTROL_VARIANT,
			int(high['variant'])]))


if __name__ == '__main__':
	Plot().cli()
