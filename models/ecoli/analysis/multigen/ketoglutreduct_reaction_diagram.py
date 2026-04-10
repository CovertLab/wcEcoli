"""
Static pathway diagram showing the KETOGLUTREDUCT-RXN reaction family and how
they are linked through shared metabolites and enzymes.

Single panel using matplotlib text and FancyArrowPatch — no simulation data
needed, but inherits from MultigenAnalysisPlot for framework consistency.

Reactions:
  KETOGLUTREDUCT-RXN  (PGLYCDEHYDROG-CPLX / SerA)
  RXN-16701           (PGLYCDEHYDROG-CPLX / SerA)
  PGLYCDEHYDROG-RXN   (PGLYCDEHYDROG-CPLX / SerA)
  RXN-14932           (CPLX0-9749 / YdiJ)

Metabolites color-coded: blue = homeostatic target, gray = kinetic-only or
neither.
"""

import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure

# Homeostatic targets (blue); all others gray.
HOMEOSTATIC = {'NAD[c]', 'NADH[c]', '2-KETOGLUTARATE[c]'}


def _draw_box(ax, xy, text, width=0.22, height=0.07, fc='#e8e8e8',
			  ec='black', fontsize=7, bold=False):
	"""Draw a labeled rounded box and return its center coordinates."""
	x, y = xy
	box = FancyBboxPatch(
		(x - width / 2, y - height / 2), width, height,
		boxstyle='round,pad=0.01', fc=fc, ec=ec, lw=1.0,
		transform=ax.transData, zorder=2)
	ax.add_patch(box)
	weight = 'bold' if bold else 'normal'
	ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
			weight=weight, zorder=3)
	return x, y


def _draw_met(ax, xy, text, fontsize=6.5):
	"""Draw a metabolite label with color based on homeostatic status."""
	x, y = xy
	color = '#1f77b4' if text in HOMEOSTATIC else '#666666'
	suffix = ' (H)' if text in HOMEOSTATIC else ' (K)' if text != 'PROTON[c]' else ''
	ax.text(x, y, text + suffix, ha='center', va='center', fontsize=fontsize,
			color=color, zorder=3,
			bbox=dict(boxstyle='round,pad=0.02', fc='white', ec=color,
					  lw=0.5, alpha=0.8))
	return x, y


def _arrow(ax, start, end, bidirectional=False, color='#444444'):
	"""Draw an arrow (or bidirectional arrow) between two points."""
	style = 'Simple,tail_width=0.5,head_width=4,head_length=3'
	arrow = FancyArrowPatch(
		start, end,
		arrowstyle='<->' if bidirectional else '->',
		color=color, lw=1.0, mutation_scale=10,
		connectionstyle='arc3,rad=0', zorder=1)
	ax.add_patch(arrow)


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		fig, ax = plt.subplots(1, 1, figsize=(12, 8))
		ax.set_xlim(-0.05, 1.05)
		ax.set_ylim(-0.05, 1.05)
		ax.set_aspect('equal')
		ax.axis('off')

		# --- Reaction boxes ---
		rx1 = _draw_box(ax, (0.25, 0.75), 'KETOGLUTREDUCT-RXN', width=0.24,
						fc='#fff3cd', bold=True)
		rx2 = _draw_box(ax, (0.75, 0.75), 'RXN-16701', width=0.18,
						fc='#fff3cd', bold=True)
		rx3 = _draw_box(ax, (0.25, 0.25), 'PGLYCDEHYDROG-RXN', width=0.22,
						fc='#d4edda')
		rx4 = _draw_box(ax, (0.75, 0.25), 'RXN-14932', width=0.16,
						fc='#d4edda')

		# Enzyme labels
		ax.text(0.25, 0.67, 'Enzyme: PGLYCDEHYDROG-CPLX (SerA)',
				ha='center', fontsize=5.5, style='italic', color='#8B4513')
		ax.text(0.75, 0.67, 'Enzyme: PGLYCDEHYDROG-CPLX (SerA)',
				ha='center', fontsize=5.5, style='italic', color='#8B4513')
		ax.text(0.25, 0.17, 'Enzyme: PGLYCDEHYDROG-CPLX (SerA)',
				ha='center', fontsize=5.5, style='italic', color='#8B4513')
		ax.text(0.75, 0.17, 'Enzyme: CPLX0-9749 (YdiJ)',
				ha='center', fontsize=5.5, style='italic', color='#8B4513')

		# --- Metabolite nodes ---
		akg = _draw_met(ax, (0.50, 0.88), '2-KETOGLUTARATE[c]')
		r2hg = _draw_met(ax, (0.15, 0.50), 'R-2-HYDROXYGLUTARATE[c]')
		cpd381 = _draw_met(ax, (0.85, 0.50), 'CPD-381[c]')
		nad = _draw_met(ax, (0.50, 0.58), 'NAD[c]')
		nadh = _draw_met(ax, (0.50, 0.48), 'NADH[c]')
		proton = _draw_met(ax, (0.50, 0.38), 'PROTON[c]')
		g3p = _draw_met(ax, (0.08, 0.25), 'G3P[c]')
		p3hp = _draw_met(ax, (0.42, 0.12), '3-P-HYDROXYPYRUVATE[c]')

		# --- Arrows: KETOGLUTREDUCT-RXN ---
		# 2-KG + NADH + H+ -> R-2-HG + NAD (reversible)
		_arrow(ax, (0.50, 0.85), (0.37, 0.78), bidirectional=True)  # 2-KG <-> rxn1
		_arrow(ax, (0.15, 0.53), (0.18, 0.71), bidirectional=True)  # R-2-HG <-> rxn1
		_arrow(ax, (0.47, 0.56), (0.35, 0.78), bidirectional=True, color='#1f77b4')  # NAD/NADH <-> rxn1

		# --- Arrows: RXN-16701 ---
		# 2-KG + NADH + H+ -> CPD-381 + NAD (reversible)
		_arrow(ax, (0.53, 0.85), (0.70, 0.78), bidirectional=True)  # 2-KG <-> rxn2
		_arrow(ax, (0.85, 0.53), (0.82, 0.71), bidirectional=True)  # CPD-381 <-> rxn2
		_arrow(ax, (0.53, 0.56), (0.68, 0.78), bidirectional=True, color='#1f77b4')  # NAD/NADH <-> rxn2

		# --- Arrows: PGLYCDEHYDROG-RXN ---
		# G3P + NAD -> 3-P-HYDROXYPYRUVATE + NADH + H+
		_arrow(ax, (0.12, 0.25), (0.18, 0.25), bidirectional=True)  # G3P <-> rxn3
		_arrow(ax, (0.39, 0.14), (0.32, 0.22), bidirectional=True)  # 3-PHP <-> rxn3
		_arrow(ax, (0.47, 0.45), (0.32, 0.28), bidirectional=True, color='#1f77b4')  # NAD/NADH <-> rxn3

		# --- Arrows: RXN-14932 ---
		# R-2-HG + Acceptor -> 2-KG + Donor-H2
		_arrow(ax, (0.20, 0.47), (0.68, 0.28))  # R-2-HG -> rxn4
		_arrow(ax, (0.82, 0.28), (0.55, 0.85))  # rxn4 -> 2-KG

		# --- Stoichiometry annotations ---
		ax.text(0.25, 0.82, 'fwd: 2-KG + NADH + H+ → R-2-HG + NAD',
				ha='center', fontsize=5, color='#555555')
		ax.text(0.75, 0.82, 'fwd: 2-KG + NADH + H+ → CPD-381 + NAD',
				ha='center', fontsize=5, color='#555555')
		ax.text(0.25, 0.30, 'G3P + NAD → 3-PHP + NADH + H+',
				ha='center', fontsize=5, color='#555555')
		ax.text(0.75, 0.30, 'R-2-HG + Acceptor → 2-KG + Donor-H2',
				ha='center', fontsize=5, color='#555555')

		# --- Legend ---
		ax.text(0.50, 0.02,
				'Blue text = homeostatic target (H)  |  '
				'Gray text = kinetic-only (K) or neither  |  '
				'Yellow box = futile cycle reactions  |  '
				'Green box = context reactions',
				ha='center', fontsize=6, color='#333333',
				bbox=dict(boxstyle='round', fc='#f0f0f0', ec='gray', lw=0.5))

		ax.set_title('KETOGLUTREDUCT-RXN reaction family — pathway diagram',
					 fontsize=11, pad=10)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
