"""
Gene copies lost under burden, by gene group, against the terminus.

# Figure 5E

For each group, copies are the summed promoter copies (RnaSynthProb
promoter_copy_number) over the group's transcription units, window means.
Loss = 1 - group copies with high GFP / the same seed's no-GFP copies, mean
over surviving high-GFP lineages with a t-based 95% CI. The baseline is the
no-GFP control, the same seed-matched baseline as every other paper panel.

Groups are defined at the transcription-unit level from sim_data flags:
	RNA polymerase subunits   rna_data includes_RNAP
	rRNA operons              rna_data is_rRNA
	ribosomal proteins        rna_data includes_ribosomal_protein
	GFP construct             units carrying the inserted gene
	terminus reference        gfp_paper_data.terminus_tus (beyond
	                          TERMINUS_MIN_FRACTION, machinery excluded)
RNAP and ribosomal-protein units can overlap where they share an operon.

f under each label is the group's mean replichore fraction (0 = oriC). The
right-hand column is each group's loss over the terminus group's loss.
The GFP row's f follows the insertion site, so it moves with INSERTION_POS.

The supporting claim: machinery genes sit near oriC, so they lose the most
copies when slower growth reduces overlapping replication.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0


def gene_groups(sim_data):
	"""[(label, colour, [TU ids])], in display order."""
	transcription = sim_data.process.transcription
	rna = transcription.rna_data
	ids = rna['id']
	new_gene = set()
	cistrons = transcription.cistron_data
	for cistron in cistrons['id'][cistrons['is_new_gene']]:
		new_gene.update(transcription.cistron_id_to_rna_indexes(cistron))
	return [
		('RNA polymerase subunits', style.COLORS['rnap'],
			list(ids[rna['includes_RNAP']])),
		('rRNA operons', style.COLORS['rrna'], list(ids[rna['is_rRNA']])),
		('ribosomal proteins', style.COLORS['ribosome'],
			list(ids[rna['includes_ribosomal_protein']])),
		('GFP construct', style.COLORS['gfp'],
			[ids[i] for i in sorted(new_gene)]),
		('terminus reference', style.COLORS['control'],
			gfp_paper_data.terminus_tus(sim_data)),
	]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		from scipy import stats
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		high = gfp_paper_data.nearest_design(data.designs,
			C.FIG5_PROTEOME_SHARES[-1])
		control = gfp_paper_data.CONTROL_VARIANT
		profiles = gfp_paper_data.vector_profiles(data, 'promoter_copies',
			[control, int(high['variant'])], plotOutDir, self.cpus)
		base, burden = profiles[control], profiles[int(high['variant'])]
		rna = sim_data.process.transcription.rna_data
		coordinate = dict(zip(rna['id'], rna['replication_coordinate']))

		rows = []
		for label, color, tus in gene_groups(sim_data):
			cols = [base.ids.index(t) for t in tus if t in base.ids]
			seeds = sorted(set(base.per_lineage) & set(burden.per_lineage))
			losses = np.array([1 - burden.per_lineage[s][cols].sum()
				/ base.per_lineage[s][cols].sum() for s in seeds
				if base.per_lineage[s][cols].sum() > 0])
			n = losses.size
			half = (stats.t.ppf(0.975, n - 1) * losses.std(ddof=1)
				/ np.sqrt(n)) if n > 1 else np.nan
			f = gfp_paper_data.replichore_fraction(sim_data,
				[coordinate[t] for t in tus]).mean() if tus else np.nan
			rows.append(dict(label=label, color=color, n_tus=len(cols),
				f=f, loss=losses.mean() if n else np.nan, ci95=half, n=n))
		terminus_loss = rows[-1]['loss']

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
		y = np.arange(len(rows))[::-1]
		for yi, r in zip(y, rows):
			ax.barh(yi, PERCENT * r['loss'], height=0.55, color=r['color'],
				xerr=PERCENT * r['ci95'] if np.isfinite(r['ci95']) else None,
				error_kw=dict(ecolor=style.INK, elinewidth=0.7, capsize=1.5))
			end = PERCENT * (r['loss'] + np.nan_to_num(r['ci95']))
			ax.annotate('%.1f%%' % (PERCENT * r['loss']), xy=(end, yi),
				xytext=(3, 0), textcoords='offset points', va='center',
				fontsize=style.FONT_SIZES['annotation'], color=style.INK)
			ax.annotate('%.2f×' % (r['loss'] / terminus_loss),
				xy=(1.0, yi), xycoords=('axes fraction', 'data'),
				xytext=(28, 0), textcoords='offset points', va='center',
				ha='right', fontsize=style.FONT_SIZES['annotation'],
				color=style.INK)
		ax.annotate('× terminus', xy=(1.0, 1.0),
			xycoords='axes fraction', xytext=(28, 4),
			textcoords='offset points', ha='right', va='bottom',
			fontsize=style.FONT_SIZES['annotation'], color=style.INK)
		ax.set_yticks(y)
		ax.set_yticklabels(['%s\nf = %.3f' % (r['label'], r['f'])
			for r in rows])
		ax.tick_params(axis='y', length=0)
		ax.spines['left'].set_visible(False)
		ax.set_xlim(0, max(30, PERCENT * max(
			r['loss'] + np.nan_to_num(r['ci95']) for r in rows) * 1.2))
		ax.set_xlabel('gene copies lost under burden (%)')

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=min(r['n'] for r in rows),
			notes='High GFP = variant %d (GFP %.3f of proteome); baseline ='
				' no-GFP control, seed-matched. %s.' % (high['variant'],
				high['gfp_proteome_share'], '; '.join(
				'%s: %d TUs, f %.3f, loss %.4f +/- %.4f'
				% (r['label'], r['n_tus'], r['f'], r['loss'], r['ci95'])
				for r in rows)),
			**data.provenance([control, int(high['variant'])]))


if __name__ == '__main__':
	Plot().cli()
