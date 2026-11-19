"""
One lineage's GFP share, cell mass and replication origins through induction.

# Figure 5A

Three stacked tracks sharing a time axis, from REPLICATION_TIMELINE_FIRST_GEN
(four wildtype generations before induction) to REPLICATION_TIMELINE_LAST_GEN:
	GFP share of proteome (%)   GFP mass / protein mass, every time step
	cell mass (fg)              with a marker at each replication initiation
	replication origins         number of oriC, every time step

A replication initiation is a time step at which the number of origins rises
within a cell; the halving at division is between cells and is not counted.
The dashed line is induction (the first time step of NEW_GENE_INDUCTION_GEN).

This is a multigen analysis, so it runs once per lineage. Every lineage is
written to <batch>/plotOut/paper_figures/Figure_5A_lineages/, named by variant
and seed, and the representative one is chosen from that set; the choice is
a figure decision, not a code one. A lineage that died is drawn up to its
last time step and flagged in the sidecar.

The supporting claim: once GFP slows growth, replication initiates later in
the cell cycle, and origins drop from the 2-4 rhythm to 1-2.
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
PERCENT = 100.0
SUBDIR = 'Figure_5A_lineages'


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		gene = gfp_paper_data.gene_info(sim_data)

		cells = self.ap.get_cells(generation=np.arange(
			C.REPLICATION_TIMELINE_FIRST_GEN,
			C.REPLICATION_TIMELINE_LAST_GEN + 1))
		tracks = [self._read(os.path.join(c, 'simOut'), gene) for c in cells]
		tracks = [t for t in tracks if t is not None]
		if not tracks:
			print('No readable generations in the timeline window; skipping.')
			return

		time = np.concatenate([t['time'] for t in tracks])
		t0 = time[0]
		hours = lambda x: (np.asarray(x) - t0) / SECONDS_PER_HOUR
		induction = next((t['time'][0] for t in tracks
			if t['generation'] == C.NEW_GENE_INDUCTION_GEN), None)

		style.apply()
		fig, axes = plt.subplots(3, 1, sharex=True,
			figsize=style.panel_size('double', 0.5),
			gridspec_kw={'height_ratios': [1, 1, 0.8]})
		share_ax, mass_ax, ori_ax = axes
		for t in tracks:
			x = hours(t['time'])
			share_ax.plot(x, PERCENT * t['gfp_share'],
				color=style.COLORS['gfp'], linewidth=1.0)
			mass_ax.plot(x, t['cell_mass'], color=style.INK, linewidth=0.8)
			ori_ax.fill_between(x, t['oric'], step='post', linewidth=0,
				color=style.COLORS['replication'], alpha=0.85)
		events = [(hours(t['time'][i]), t['cell_mass'][i]) for t in tracks
			for i in t['initiations']]
		if events:
			ex, ey = zip(*events)
			mass_ax.scatter(ex, ey, s=16, color=style.COLORS['replication'],
				edgecolors='white', linewidths=0.6, zorder=3,
				label='replication initiation')

		if induction is not None:
			for ax in axes:
				ax.axvline(hours(induction), color=style.INK_MUTED,
					linestyle=(0, (4, 3)), linewidth=0.8)
			share_ax.annotate('induction', xy=(hours(induction), 1),
				xycoords=('data', 'axes fraction'), xytext=(3, -2),
				textcoords='offset points', va='top',
				fontsize=style.FONT_SIZES['annotation'],
				color=style.INK_MUTED)

		share_ax.set_ylabel('GFP share of\nproteome (%)')
		share_ax.set_ylim(bottom=0)
		mass_ax.set_ylabel('cell mass (fg)')
		mass_ax.set_ylim(bottom=0)
		ori_ax.set_ylabel('replication\norigins')
		top = int(max(t['oric'].max() for t in tracks))
		ori_ax.set_yticks(sorted({1, 2, top} | ({4} if top >= 4 else set())))
		ori_ax.set_ylim(0, top + 0.5)
		ori_ax.set_xlabel('time (h)')
		ori_ax.set_xlim(0, hours(time[-1]))
		if events:
			mass_ax.legend(loc='upper right')

		batch = os.path.dirname(os.path.dirname(seedOutDir))
		variant = int(os.path.basename(os.path.dirname(seedOutDir))[-6:])
		seed = int(os.path.basename(seedOutDir))
		generations = [t['generation'] for t in tracks]
		reached_end = (generations[-1] == C.REPLICATION_TIMELINE_LAST_GEN
			and tracks[-1]['divided'])
		style.save_panel(fig, os.path.join(batch, 'plotOut'),
			style.panel_name(__file__, 'v%06d_s%06d' % (variant, seed)),
			script=__file__, batch=batch, sim_data_path=simDataFile,
			variants=[variant], seeds=[seed],
			generation_window=(generations[0], generations[-1]),
			n=len(events), n_lineages_total=1,
			n_lineages_survived=int(reached_end), subdir=SUBDIR,
			notes='%d replication initiations over generations %d-%d.%s'
				% (len(events), generations[0], generations[-1],
				'' if reached_end else ' Lineage did not reach the end of'
				' the window: it failed or was cut short.'))

	@staticmethod
	def _read(sim_out, gene):
		"""Per-time-step tracks for one cell, or None if it has no data."""
		try:
			time = TableReader(os.path.join(sim_out, 'Main')).readColumn(
				'time')
		except Exception:                                   # noqa: BLE001
			return None
		if time.size == 0:
			return None
		mass = TableReader(os.path.join(sim_out, 'Mass'))
		monomers = TableReader(os.path.join(sim_out, 'MonomerCounts'))
		index = list(monomers.readAttribute('monomerIds')).index(
			gene.monomer_id)
		gfp = monomers.readColumn('monomerCounts')[:, index]
		protein = mass.readColumn('proteinMass')
		oric = TableReader(os.path.join(sim_out, 'ReplicationData')
			).readColumn('numberOfOric')
		generation = int(sim_out.split('generation_')[1][:6])
		with np.errstate(divide='ignore', invalid='ignore'):
			share = np.where(protein > 0, gfp * gene.monomer_mass_fg
				/ protein, np.nan)
		return dict(time=time, generation=generation,
			cell_mass=mass.readColumn('cellMass'), gfp_share=share,
			oric=oric, initiations=np.where(np.diff(oric) > 0)[0] + 1,
			divided=os.path.exists(os.path.join(sim_out,
				'Daughter1_inherited_state.cPickle')))


if __name__ == '__main__':
	Plot().cli()
