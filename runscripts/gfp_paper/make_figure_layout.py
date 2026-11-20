#!/usr/bin/env python
"""
Write figure_layout.tsv, the position and size of every panel, from the
locked storyboard deck.

The deck has one page per figure with one picture per panel (alt text
"Panel A", "Panel B", ...). Each figure's panels are scaled together so the
widest extent fills FIGURE_WIDTH_MM, the journal's double column, and
shifted so the figure starts at (0, 0). Only geometry is written: no slide
text, captions or images, so the layout can be committed.

Rerun when the storyboard's arrangement changes; edit the TSV by hand for
small adjustments, and add a `parts` entry for a slot that holds more than
one file side by side (3C).

Usage:
	python runscripts/gfp_paper/make_figure_layout.py <storyboard.pptx>
"""

import argparse
import csv
import os
import re

FIGURE_WIDTH_MM = 183.0
MM_PER_INCH = 25.4
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
	'figure_layout.tsv')
# Slots that hold several panel files side by side, by panel label.
PARTS = {'3C': 'mrna,growth'}


def main():
	from pptx import Presentation
	from pptx.util import Emu
	parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
	parser.add_argument('deck')
	args = parser.parse_args()
	deck = Presentation(args.deck)
	rows = []
	for figure, slide in enumerate(deck.slides, start=1):
		boxes = []
		for shape in slide.shapes:
			if shape.shape_type != 13:              # pictures only
				continue
			name = shape._element.xpath('./p:nvPicPr/p:cNvPr')[0].get(
				'descr') or ''
			match = re.fullmatch(r'Panel ([A-Z])', name.strip())
			if match is None:
				continue
			boxes.append((match.group(1), Emu(shape.left).inches,
				Emu(shape.top).inches, Emu(shape.width).inches,
				Emu(shape.height).inches))
		if not boxes:
			continue
		left = min(b[1] for b in boxes)
		top = min(b[2] for b in boxes)
		right = max(b[1] + b[3] for b in boxes)
		scale = FIGURE_WIDTH_MM / ((right - left) * MM_PER_INCH)
		for letter, x, y, w, h in sorted(boxes):
			label = '%d%s' % (figure, letter)
			rows.append(dict(figure=figure, panel=label,
				x_mm='%.2f' % ((x - left) * MM_PER_INCH * scale),
				y_mm='%.2f' % ((y - top) * MM_PER_INCH * scale),
				w_mm='%.2f' % (w * MM_PER_INCH * scale),
				h_mm='%.2f' % (h * MM_PER_INCH * scale),
				parts=PARTS.get(label, '')))
	with open(OUT, 'w', newline='') as f:
		# Plain LF line endings, as git stores them; csv defaults to CRLF.
		writer = csv.DictWriter(f, delimiter='\t', lineterminator='\n',
			fieldnames=['figure', 'panel', 'x_mm', 'y_mm', 'w_mm', 'h_mm',
			'parts'])
		writer.writeheader()
		writer.writerows(rows)
	print('wrote %s: %d panels across %d figures' % (OUT, len(rows),
		len({r['figure'] for r in rows})))


if __name__ == '__main__':
	main()
