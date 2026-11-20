#!/usr/bin/env python
"""
Assemble the paper's figures from their saved panels: SVG, PNG and PPTX, plus
PDF when cairosvg is installed. Panel letters only, no captions.

Every panel script saves its panel separately (gfp_paper_style.save_panel),
named Figure_<label>_<content>[_<part>].{svg,png} under each batch's
plotOut/paper_figures/. This script finds the right file for each slot in
figure_layout.tsv, fits it into the slot's box, adds the panel letter, and
writes one file per figure in each format:

	SVG   the panels' own vector SVGs placed and scaled, text still editable
	PNG   the 600 dpi panel PNGs composed at 600 dpi
	PPTX  one slide per figure, panel PNGs and letters as separate shapes
	PDF   converted from the SVG, only if cairosvg is importable

Choosing between versions of a panel:
	--growth doubling_time   use each growth panel's doubling-time version
	--version 2D=concentration --version 1E=relative
	                         use a named part instead of the plain panel
	--lineage 5A=v000005_s000003
	                         which lineage to show for the one-lineage panels
	                         (5A, 6A); otherwise the first one found is used
	                         and the choice is reported
Schematics and other hand-made panels are read from --schematics as
Figure_<label>.svg (and .png for the raster outputs). A slot with no file
gets a dashed placeholder box, so a figure can be assembled at any stage.

Each panel is fitted to its box less a band at the top for the letter
(gfp_paper_style.LETTER_BAND_MM). save_panel draws panels at exactly that
size, so the scale is near 1; a panel saved before that, or with content
outside its figure, is scaled to fit, and since that changes its type size the
scale of every panel is printed and written to the provenance JSON, with a
warning outside SCALE_WARN.

Usage:
	python runscripts/gfp_paper/assemble_figures.py \\
		--panels <batch>/plotOut/paper_figures [<dir> ...] \\
		[--schematics DIR] [--figures 1 2 ...] [--growth rate|doubling_time] \\
		[--version LABEL=PART ...] [--lineage LABEL=vNNNNNN_sNNNNNN ...] \\
		[--out DIR]
Output defaults to out/files_for_claude/gfp_paper/figures/<timestamp>/, which
is not tracked by git.
"""

import argparse
import copy
import csv
import glob
import json
import os
import re
from datetime import datetime

from lxml import etree

from wholecell.analysis.gfp_paper_style import LETTER_BAND_MM

HERE = os.path.dirname(os.path.abspath(__file__))
LAYOUT = os.path.join(HERE, 'figure_layout.tsv')
MANIFEST = os.path.join(HERE, 'panels.tsv')
SCRIPT_PREFIX = 'gfp_paper_'
DPI = 600
MM_PER_INCH = 25.4
PT_PER_MM = 72.0 / MM_PER_INCH
LETTER_PT = 9.0
LETTER_FONT = 'Arial'
SCALE_WARN = (0.85, 1.15)
SVG_NS = 'http://www.w3.org/2000/svg'
XLINK_NS = 'http://www.w3.org/1999/xlink'
LINEAGE_PANELS = {'5A', '6A'}


# --------------------------------------------------------------------------
# Which file goes in each slot
def read_layout():
	with open(LAYOUT) as f:
		return list(csv.DictReader(f, delimiter='\t'))


def read_manifest():
	"""{label: (kind, content stem)} from panels.tsv."""
	out = {}
	with open(MANIFEST) as f:
		for row in csv.DictReader(f, delimiter='\t'):
			stem = os.path.splitext(os.path.basename(row['script']))[0]
			if stem.startswith(SCRIPT_PREFIX):
				stem = stem[len(SCRIPT_PREFIX):]
			out[row['panel']] = (row['kind'], stem)
	return out


def find(dirs, name, subdir=None):
	"""The first directory's copy of name.svg (or None), and its PNG twin."""
	for d in dirs:
		base = os.path.join(d, subdir, name) if subdir else os.path.join(
			d, name)
		if os.path.exists(base + '.svg') or os.path.exists(base + '.png'):
			svg = base + '.svg' if os.path.exists(base + '.svg') else None
			png = base + '.png' if os.path.exists(base + '.png') else None
			return svg, png
	return None, None


def resolve(label, part, manifest, dirs, schematics, growth, lineages,
		notes):
	"""(svg path, png path, chosen name) for one slot or slot part."""
	kind, stem = manifest.get(label, ('', ''))
	if kind in ('schematic', 'experimental') or not stem or stem.startswith(
			'('):
		svg, png = find([schematics] if schematics else [], 'Figure_%s'
			% label)
		return svg, png, 'Figure_%s (%s)' % (label, kind or 'hand-made')
	base = 'Figure_%s_%s' % (label, stem)
	subdir = None
	if label in LINEAGE_PANELS:
		subdir = 'Figure_%s_lineages' % label
		chosen = lineages.get(label)
		if chosen is None:
			found = sorted({os.path.basename(p)[len(base) + 1:].split(
				'_doubling_time')[0].rsplit('.', 1)[0] for d in dirs for p in
				glob.glob(os.path.join(d, subdir, base + '_v*.svg'))})
			if not found:
				return None, None, base + ' (no lineage files)'
			chosen = found[0]
			notes.append('%s: no --lineage given, used %s of %d' % (label,
				chosen, len(found)))
		part = chosen if part is None else '%s_%s' % (chosen, part)
	names = []
	want = base if not part else '%s_%s' % (base, part)
	if growth == 'doubling_time':
		if part == 'growth':
			names.append('%s_doubling_time' % base)
		names.append('%s_doubling_time' % want)
	names.append(want)
	for name in names:
		svg, png = find(dirs, name, subdir)
		if svg or png:
			return svg, png, name
	return None, None, want + ' (not found)'


# --------------------------------------------------------------------------
# SVG
def svg_size_pt(root):
	"""(width, height) of a panel SVG in its own user units (points)."""
	box = root.get('viewBox')
	if box:
		_, _, w, h = (float(v) for v in box.replace(',', ' ').split())
		return w, h
	def pt(value):
		value = value.strip()
		for unit, factor in (('pt', 1.0), ('px', 0.75), ('mm', PT_PER_MM),
				('in', 72.0)):
			if value.endswith(unit):
				return float(value[:-len(unit)]) * factor
		return float(value)
	return pt(root.get('width')), pt(root.get('height'))


def prefixed(root, prefix):
	"""A copy of a panel's content with every id made unique to it."""
	root = copy.deepcopy(root)
	ids = {el.get('id') for el in root.iter() if el.get('id')}
	href = '{%s}href' % XLINK_NS
	for el in root.iter():
		if not isinstance(el.tag, str):
			continue
		if el.get('id'):
			el.set('id', prefix + el.get('id'))
		for key, value in list(el.attrib.items()):
			if 'url(#' in value:
				el.set(key, re.sub(r'url\(#([^)]+)\)', lambda m: 'url(#%s%s)'
					% (prefix if m.group(1) in ids else '', m.group(1)),
					value))
			elif key in (href, 'href') and value.startswith('#') and (
					value[1:] in ids):
				el.set(key, '#' + prefix + value[1:])
	for el in root.iter('{%s}metadata' % SVG_NS):
		el.getparent().remove(el)
	return root


def letter_svg(parent, x_pt, y_pt, letter):
	text = etree.SubElement(parent, '{%s}text' % SVG_NS, x='%.2f' % x_pt,
		y='%.2f' % (y_pt + LETTER_PT), style='font-family:%s;font-size:%gpx;'
		'font-weight:bold' % (LETTER_FONT, LETTER_PT))
	text.text = letter


def placeholder_svg(parent, x, y, w, h, label):
	etree.SubElement(parent, '{%s}rect' % SVG_NS, x='%.2f' % x,
		y='%.2f' % y, width='%.2f' % w, height='%.2f' % h,
		style='fill:none;stroke:#9AA39F;stroke-width:0.6;'
		'stroke-dasharray:3,2')
	text = etree.SubElement(parent, '{%s}text' % SVG_NS, x='%.2f' % (
		x + w / 2), y='%.2f' % (y + h / 2), style='font-family:%s;'
		'font-size:6px;fill:#9AA39F;text-anchor:middle' % LETTER_FONT)
	text.text = label


# --------------------------------------------------------------------------
# Assembly
def slots_for(figure, layout):
	return [r for r in layout if int(r['figure']) == figure]


def assemble(figure, layout, manifest, args, out_dir):
	"""Write one figure in every format; returns its provenance."""
	from PIL import Image, ImageDraw, ImageFont
	slots = slots_for(figure, layout)
	width_mm = max(float(s['x_mm']) + float(s['w_mm']) for s in slots)
	height_mm = max(float(s['y_mm']) + float(s['h_mm']) for s in slots)
	svg = etree.Element('{%s}svg' % SVG_NS, nsmap={None: SVG_NS,
		'xlink': XLINK_NS}, width='%.2fmm' % width_mm, height='%.2fmm'
		% height_mm, viewBox='0 0 %.2f %.2f' % (width_mm * PT_PER_MM,
		height_mm * PT_PER_MM), version='1.1')
	px = lambda mm: int(round(mm / MM_PER_INCH * DPI))
	canvas = Image.new('RGB', (px(width_mm), px(height_mm)), 'white')
	draw = ImageDraw.Draw(canvas)
	try:
		font = ImageFont.truetype('/System/Library/Fonts/Supplemental/'
			'Arial Bold.ttf', int(LETTER_PT / 72 * DPI))
		small = ImageFont.truetype('/System/Library/Fonts/Supplemental/'
			'Arial.ttf', int(6 / 72 * DPI))
	except OSError:
		font = small = ImageFont.load_default()
	placed, notes = [], []

	for slot in slots:
		label = slot['panel']
		x, y = float(slot['x_mm']), float(slot['y_mm'])
		w, h = float(slot['w_mm']), float(slot['h_mm'])
		letter_y = y
		y, h = y + LETTER_BAND_MM, h - LETTER_BAND_MM
		override = args.version.get(label)
		parts = ([override] if override else
			[p for p in slot['parts'].split(',') if p] or [None])
		part_w = w / len(parts)
		for k, part in enumerate(parts):
			px0 = x + k * part_w
			svg_path, png_path, name = resolve(label, part, manifest,
				args.panels, args.schematics, args.growth, args.lineage,
				notes)
			record = dict(panel=label, part=part, file=name, svg=svg_path,
				png=png_path, box_mm=[round(px0, 2), round(y, 2),
				round(part_w, 2), round(h, 2)])
			if svg_path:
				panel = etree.parse(svg_path).getroot()
				pw, ph = svg_size_pt(panel)
				scale = min(part_w * PT_PER_MM / pw, h * PT_PER_MM / ph)
				group = etree.SubElement(svg, '{%s}g' % SVG_NS,
					transform='translate(%.3f %.3f) scale(%.5f)' % (
					px0 * PT_PER_MM, y * PT_PER_MM, scale))
				content = prefixed(panel, 'f%s_%s_' % (label, k))
				for child in content:
					group.append(child)
				record['scale'] = round(scale, 3)
				if not SCALE_WARN[0] <= scale <= SCALE_WARN[1]:
					notes.append('%s%s drawn at %.0f%% of its saved size, so'
						' its type is %.1f pt not %.1f pt; resize the panel'
						' or its box' % (label, '' if part is None else ' '
						+ part, 100 * scale, 7 * scale, 7.0))
			else:
				placeholder_svg(svg, px0 * PT_PER_MM, y * PT_PER_MM,
					part_w * PT_PER_MM, h * PT_PER_MM, name)
			if png_path:
				image = Image.open(png_path).convert('RGB')
				iw_mm = image.width / DPI * MM_PER_INCH
				ih_mm = image.height / DPI * MM_PER_INCH
				s = min(part_w / iw_mm, h / ih_mm)
				image = image.resize((max(1, px(iw_mm * s)),
					max(1, px(ih_mm * s))), Image.LANCZOS)
				canvas.paste(image, (px(px0), px(y)))
				record['png_size_mm'] = [round(iw_mm * s, 2),
					round(ih_mm * s, 2)]
			else:
				draw.rectangle([px(px0), px(y), px(px0 + part_w),
					px(y + h)], outline='#9AA39F', width=4)
				draw.text((px(px0 + part_w / 2), px(y + h / 2)), name,
					fill='#9AA39F', font=small, anchor='mm')
			placed.append(record)
		letter = label[len(str(figure)):]
		letter_svg(svg, x * PT_PER_MM, letter_y * PT_PER_MM, letter)
		draw.text((px(x), px(letter_y)), letter, fill='black', font=font)

	stem = os.path.join(out_dir, 'Figure_%d' % figure)
	etree.ElementTree(svg).write(stem + '.svg', xml_declaration=True,
		encoding='utf-8', pretty_print=False)
	canvas.save(stem + '.png', dpi=(DPI, DPI))
	try:
		import cairosvg
		cairosvg.svg2pdf(url=stem + '.svg', write_to=stem + '.pdf')
		pdf = True
	except ImportError:
		pdf = False
	provenance = dict(figure=figure, written=datetime.now().isoformat(
		timespec='seconds'), size_mm=[round(width_mm, 2),
		round(height_mm, 2)], growth=args.growth, versions=args.version,
		lineages=args.lineage, panels=placed, notes=notes, pdf=pdf)
	with open(stem + '.json', 'w') as f:
		json.dump(provenance, f, indent=1)
	return provenance


def write_pptx(figures, provenances, out_dir):
	"""One slide per figure, sized to the largest figure."""
	from pptx import Presentation
	from pptx.util import Mm, Pt
	deck = Presentation()
	deck.slide_width = Mm(max(p['size_mm'][0] for p in provenances))
	deck.slide_height = Mm(max(p['size_mm'][1] for p in provenances))
	blank = deck.slide_layouts[6]
	for figure, provenance in zip(figures, provenances):
		slide = deck.slides.add_slide(blank)
		letters = {}
		for record in provenance['panels']:
			x, y, w, h = record['box_mm']
			if record['png']:
				iw, ih = record['png_size_mm']
				slide.shapes.add_picture(record['png'], Mm(x), Mm(y),
					Mm(iw), Mm(ih))
			letters.setdefault(record['panel'], (x, y - LETTER_BAND_MM))
		for label, (x, y) in letters.items():
			box = slide.shapes.add_textbox(Mm(x), Mm(y), Mm(8), Mm(6))
			frame = box.text_frame
			frame.margin_left = frame.margin_top = 0
			run = frame.paragraphs[0].add_run()
			run.text = label[len(str(figure)):]
			run.font.size = Pt(LETTER_PT)
			run.font.bold = True
			run.font.name = LETTER_FONT
	path = os.path.join(out_dir, 'figures.pptx')
	deck.save(path)
	return path


def pairs(values):
	out = {}
	for value in values or []:
		key, _, chosen = value.partition('=')
		if not chosen:
			raise SystemExit('expected LABEL=VALUE, got %r' % value)
		out[key] = chosen
	return out


def main():
	parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
	parser.add_argument('--panels', nargs='+', required=True,
		help='directories holding saved panels, searched in order')
	parser.add_argument('--schematics', default=None)
	parser.add_argument('--figures', nargs='*', type=int, default=None)
	parser.add_argument('--growth', choices=('rate', 'doubling_time'),
		default='rate')
	parser.add_argument('--version', action='append', default=[])
	parser.add_argument('--lineage', action='append', default=[])
	parser.add_argument('--out', default=None)
	args = parser.parse_args()
	args.version = pairs(args.version)
	args.lineage = pairs(args.lineage)
	layout = read_layout()
	manifest = read_manifest()
	figures = args.figures or sorted({int(r['figure']) for r in layout})
	out_dir = args.out or os.path.join(HERE, '..', '..', 'out',
		'files_for_claude', 'gfp_paper', 'figures',
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	os.makedirs(out_dir, exist_ok=True)
	provenances = []
	for figure in figures:
		provenance = assemble(figure, layout, manifest, args, out_dir)
		provenances.append(provenance)
		missing = [r['file'] for r in provenance['panels']
			if not (r['svg'] or r['png'])]
		print('Figure %d: %d slots, %d missing%s' % (figure,
			len(provenance['panels']), len(missing), (': ' + ', '.join(
			missing)) if missing else ''))
		for note in provenance['notes']:
			print('  note: ' + note)
	print('wrote', write_pptx(figures, provenances, out_dir))
	print('output in', os.path.abspath(out_dir))


if __name__ == '__main__':
	main()
