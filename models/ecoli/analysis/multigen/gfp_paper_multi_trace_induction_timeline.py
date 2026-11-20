"""
One lineage through induction: GFP, machinery, mass, replication, dosage and
growth on one time axis.

# Figure 6A

Nine stacked tracks, every time step, from LOOP_TIMELINE_FIRST_GEN to
LOOP_TIMELINE_LAST_GEN: gfp mRNA (physical), GFP protein, RNAP and ribosomes
(totals, active plus free), cell mass, replication origins, rRNA operon
copies (promoter copies summed over the rRNA transcription units), gfp copies
and growth. Each track has its own scale, marked by its range at the right,
since the point is when each one moves, not its size. Thin vertical lines are
divisions; the dashed line and shading mark induction onward.

Multigen, so it runs per lineage; every lineage is written to
<batch>/plotOut/paper_figures/Figure_6A_lineages/ named by variant and seed,
and the representative one is chosen from that set.

Two versions while decision D4 is open: the growth track is the
instantaneous growth rate (plain name), smoothed as every paper trace is
(gfp_paper_data.GROWTH_TRACE_SMOOTH_STEPS, a 5 minute rolling median, because
small-molecule pool refills cause bursts of a minute or two), or each cycle's
doubling time as a flat step (suffix DOUBLING_TIME_PART).

The supporting claim: GFP's burden spreads through machinery, mass,
replication and growth in a traceable order.
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style
from wholecell.io.tablereader import TableReader

SECONDS_PER_HOUR = 3600.0
SUBDIR = 'Figure_6A_lineages'
TRACKS = (
	('gfp_mrna', '$\\it{gfp}$ mRNA', 'gfp'),
	('gfp_protein', 'GFP protein', 'gfp'),
	('rnap', 'RNAP counts', 'rnap'),
	('ribosomes', 'ribosome counts', 'ribosome'),
	('cell_mass', 'cell mass', None),
	('oric', 'replication origins', None),
	('rrna_copies', 'rRNA operon copies', 'rrna'),
	('gfp_copies', '$\\it{gfp}$ copies', 'gfp'),
	('growth', 'growth rate', None),
)


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		gene = gfp_paper_data.gene_info(sim_data)
		machinery = gfp_paper_data.machinery_info(sim_data)
		cells = self.ap.get_cells(generation=np.arange(
			C.LOOP_TIMELINE_FIRST_GEN, C.LOOP_TIMELINE_LAST_GEN + 1))
		cycles = []
		for cell in cells:
			sim_out = os.path.join(cell, 'simOut')
			try:
				time = TableReader(os.path.join(sim_out, 'Main')).readColumn(
					'time')
			except Exception:                               # noqa: BLE001
				continue
			if time.size < 2:
				continue
			series = {name: gfp_paper_data.trace_series(sim_out, name, gene,
				machinery) for name, _, _ in TRACKS if name != 'growth'}
			series['growth_rate'] = gfp_paper_data.trace_series(sim_out,
				'growth_rate', gene, machinery)
			cycles.append(dict(time=time, generation=int(sim_out.split(
				'generation_')[1][:6]), **series))
		if not cycles:
			print('No readable generations in the 6A window; skipping.')
			return
		t0 = cycles[0]['time'][0]
		hours = lambda t: (np.asarray(t) - t0) / SECONDS_PER_HOUR
		induction = next((hours(c['time'][0]) for c in cycles
			if c['generation'] == C.NEW_GENE_INDUCTION_GEN), None)
		end = hours(cycles[-1]['time'][-1])

		batch = os.path.dirname(os.path.dirname(seedOutDir))
		variant = int(os.path.basename(os.path.dirname(seedOutDir))[-6:])
		seed = int(os.path.basename(seedOutDir))
		for part in (None, C.DOUBLING_TIME_PART):
			self._version(cycles, hours, induction, end, part, batch,
				variant, seed, simDataFile)

	@staticmethod
	def _version(cycles, hours, induction, end, part, batch, variant, seed,
			simDataFile):
		style.apply()
		fig, axes = plt.subplots(len(TRACKS), 1, sharex=True,
			figsize=style.panel_size('double', 0.45))
		for ax, (name, label, entity) in zip(axes, TRACKS):
			color = style.COLORS[entity] if entity else style.INK
			values = []
			for c in cycles:
				x = hours(c['time'])
				if name == 'growth' and part is not None:
					dt = (c['time'][-1] - c['time'][0]) / 60.0
					y = np.full(x.shape, dt)
				elif name == 'growth':
					y = c['growth_rate']
				else:
					y = c[name]
				if y is None:
					continue
				ax.plot(x, y, color=color, linewidth=0.8)
				values.append(y[np.isfinite(y)])
				ax.axvline(x[-1], color=style.GRID, linewidth=0.5, zorder=0)
			if induction is not None:
				ax.axvspan(induction, end, color=style.COLORS['gfp'],
					alpha=0.06, linewidth=0, zorder=0)
				ax.axvline(induction, color=style.INK_MUTED,
					linestyle=(0, (4, 3)), linewidth=0.8)
			track_label = ('doubling time' if name == 'growth'
				and part is not None else label)
			ax.set_ylabel(track_label, rotation=0, ha='right', va='center',
				fontsize=style.FONT_SIZES['tick'])
			ax.set_yticks([])
			ax.spines['left'].set_visible(False)
			if ax is not axes[-1]:
				ax.tick_params(axis='x', bottom=False)
			if values:
				allv = np.concatenate(values)
				lo, hi = float(np.nanmin(allv)), float(np.nanmax(allv))
				ax.annotate(_range_text(name, lo, hi, part), xy=(1, 0.5),
					xycoords='axes fraction', xytext=(4, 0),
					textcoords='offset points', va='center',
					fontsize=style.FONT_SIZES['annotation'] - 1,
					color=style.INK_MUTED)
		if induction is not None:
			axes[0].annotate('induction', xy=(induction, 1),
				xycoords=('data', 'axes fraction'), xytext=(3, 3),
				textcoords='offset points', va='bottom',
				fontsize=style.FONT_SIZES['annotation'],
				color=style.INK_MUTED)
		axes[-1].set_xlim(0, end)
		axes[-1].set_xlabel('time (h)')
		generations = [c['generation'] for c in cycles]
		style.save_panel(fig, os.path.join(batch, 'plotOut'),
			style.panel_name(__file__, ('v%06d_s%06d' % (variant, seed))
			+ ('' if part is None else '_' + part)), script=__file__,
			batch=batch, sim_data_path=simDataFile, variants=[variant],
			seeds=[seed], generation_window=(generations[0],
			generations[-1]), n=len(cycles), subdir=SUBDIR,
			notes='Each track on its own scale; ranges at the right. rRNA'
				' operon copies = promoter copies summed over the rRNA'
				' transcription units.')


def _range_text(name, lo, hi, part):
	"""The track's range, in its own units."""
	unit = {'cell_mass': ' fg', 'growth': ' /h' if part is None else ' min'}
	fmt = '%.2f' if name == 'growth' and part is None else '{:,.0f}'
	if '{' in fmt:
		return '%s-%s%s' % (fmt.format(lo), fmt.format(hi), unit.get(name, ''))
	return '%s-%s%s' % (fmt % lo, fmt % hi, unit.get(name, ''))


if __name__ == '__main__':
	Plot().cli()
