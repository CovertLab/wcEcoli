"""
Interactive protein-memory scatter + outlier analysis.

Joins the per-gene protein-memory table (subgen_protein_memory.py output) to each
protein's translation efficiency and degradation rate from sim_data, then:

  * writes a self-contained interactive HTML scatter (hover a point for its gene
    name, transcript-off frequency, protein-absence rate, def-5 rate, translation
    efficiency, and half-life; optional color-by half-life / translation
    efficiency), and
  * analyzes the outliers -- especially genes ABOVE the y = x line (protein absent
    a larger fraction of time than the gene is transcriptionally off, i.e. the
    OPPOSITE of memory) -- testing whether they are explained by short protein
    half-life (fast degradation) and/or low translation efficiency.

Axes (both per-seed rates, from subgen_protein_memory_pergene.tsv):
  x = transcript-off frequency (fraction of cell cycles with 0 completed transcripts)
  y = protein-absence rate (fraction of cell-cycle time with 0 protein copies)
  memory_gap = x - y  (positive = memory: protein present more often than transcribed)

Inputs
------
  pergene_memory.tsv : subgen_protein_memory_pergene.tsv (x, y, def5_rate per gene).
  sim_data.cPickle   : for translation_efficiencies_by_monomer and
      monomer_data['deg_rate'] (1/s) and gene common names.

Usage
-----
  python subgen_protein_memory_interactive.py <pergene_memory.tsv> \
      <sim_data.cPickle> <output_dir>

Self-contained: numpy + scipy + the sim_data pickle.
"""

import os
import csv
import json
import math
import pickle
import argparse

import numpy as np
from scipy import stats

from wholecell.utils import units

ACCENT = '#1667B8'
OUTLIER = '#E8A33D'


def load_memory(path):
	rows = []
	with open(path) as f:
		for r in csv.DictReader(f, delimiter='\t'):
			rows.append({
				'gene_id': r['gene_id'],
				'cistron_id': r['cistron_id'],
				'x': float(r['transcript_off_freq_mean']),
				'y': float(r['protein_zero_frac_mean']),
				'def5_rate': float(r['def5_rate_mean']),
				})
	return rows


def annotate(rows, sim_data):
	"""Attach monomer_id, translation efficiency, deg_rate, half-life, and gene
	name to each row (matched on cistron_id)."""
	tl = sim_data.process.translation
	md = tl.monomer_data
	mono_cistron = np.asarray(md['cistron_id'])
	mono_id = np.asarray(md['id'])
	deg_rate = md['deg_rate'].asNumber(1 / units.s)
	trl_eff = np.asarray(tl.translation_efficiencies_by_monomer)
	cis_to_i = {c: i for i, c in enumerate(mono_cistron)}
	cn = sim_data.common_names
	out = []
	for r in rows:
		i = cis_to_i.get(r['cistron_id'])
		if i is None:
			continue
		dr = float(deg_rate[i])
		hl_min = math.log(2) / dr / 60 if dr > 0 else float('inf')
		try:
			name = cn.get_common_name(r['gene_id']) or r['gene_id']
		except Exception:
			name = r['gene_id']
		r = dict(r)
		r.update({'monomer_id': str(mono_id[i]), 'trl_eff': float(trl_eff[i]),
			'deg_rate': dr, 'half_life_min': hl_min, 'name': name})
		out.append(r)
	return out


# ------------------------------------------------------------------ analysis

def analyze_outliers(rows):
	x = np.array([r['x'] for r in rows])
	y = np.array([r['y'] for r in rows])
	gap = x - y                                   # + = memory
	hl = np.array([r['half_life_min'] for r in rows])
	te = np.array([r['trl_eff'] for r in rows])
	finite_hl = np.isfinite(hl)
	log_hl = np.log10(np.where(finite_hl, hl, np.nan))
	log_te = np.log10(te)

	def corr(a, b, mask):
		m = mask & np.isfinite(a) & np.isfinite(b)
		return {'pearson_r': float(stats.pearsonr(a[m], b[m])[0]),
			'spearman_r': float(stats.spearmanr(a[m], b[m])[0]), 'n': int(m.sum())}

	all_mask = np.ones(len(rows), dtype=bool)
	res = {
		'n_genes': len(rows),
		'n_above_line': int((y > x).sum()),
		'median_gap': float(np.median(gap)),
		'corr_gap_vs_log_halflife': corr(gap, log_hl, all_mask),
		'corr_gap_vs_log_trleff': corr(gap, log_te, all_mask),
		'corr_y_vs_log_halflife': corr(y, log_hl, all_mask),
		'corr_y_vs_log_trleff': corr(y, log_te, all_mask),
		}

	# Standardized regression: gap ~ log_halflife + log_trleff (relative weight).
	m = finite_hl
	X = np.column_stack([_z(log_hl[m]), _z(log_te[m]), np.ones(m.sum())])
	beta, *_ = np.linalg.lstsq(X, _z(gap[m]), rcond=None)
	res['std_regression_gap'] = {'beta_log_halflife': float(beta[0]),
		'beta_log_trleff': float(beta[1]), 'n': int(m.sum())}

	# Population medians (context for the outlier characterization).
	res['population_median_halflife_min'] = float(np.nanmedian(hl[finite_hl]))
	res['population_median_trleff'] = float(np.median(te))

	# Memory gap by translation-efficiency band (fixed thresholds -- quantiles are
	# degenerate because most genes share the default efficiency). Tests whether
	# very low translation efficiency shrinks the memory gap / flips genes above.
	bands = [(0.0, 0.1), (0.1, 0.5), (0.5, 1.5), (1.5, np.inf)]
	res['gap_by_trleff_band'] = []
	for lo, hi in bands:
		mq = (te >= lo) & (te < hi)
		if not mq.any():
			continue
		res['gap_by_trleff_band'].append({
			'trleff_range': [lo, (None if np.isinf(hi) else hi)],
			'n': int(mq.sum()), 'mean_gap': round(float(gap[mq].mean()), 3),
			'mean_protein_absence': round(float(y[mq].mean()), 3),
			'n_above_line': int((y[mq] > x[mq]).sum())})

	def te_pctile(v):
		return round(float(100 * (te < v).mean()), 1)

	# The above-line genes and the least-memory tail, characterized.
	order = np.argsort(gap)                        # smallest gap first
	def describe(i):
		r = rows[i]
		return {'name': r['name'], 'gene_id': r['gene_id'], 'x': round(r['x'], 3),
			'y': round(r['y'], 3), 'gap': round(float(gap[i]), 3),
			'def5_rate': round(r['def5_rate'], 4),
			'trl_eff': round(r['trl_eff'], 3),
			'trl_eff_percentile': te_pctile(r['trl_eff']),
			'half_life_min': (round(r['half_life_min'], 1)
				if np.isfinite(r['half_life_min']) else None)}
	res['above_line_genes'] = [describe(i) for i in np.where(y > x)[0]]
	res['least_memory_20'] = [describe(i) for i in order[:20]]
	res['most_memory_20'] = [describe(i) for i in order[::-1][:20]]
	return res


def _z(a):
	return (a - np.nanmean(a)) / np.nanstd(a)


# ------------------------------------------------------------------ interactive HTML

def write_html(rows, path):
	data = [{'n': r['name'], 'g': r['gene_id'], 'x': round(r['x'], 4),
		'y': round(r['y'], 4), 'r': round(r['def5_rate'], 4),
		't': round(r['trl_eff'], 3),
		'h': (round(r['half_life_min'], 1)
			if np.isfinite(r['half_life_min']) else None)} for r in rows]
	payload = json.dumps(data, separators=(',', ':'))
	html = _HTML_TEMPLATE.replace('__DATA__', payload) \
		.replace('__ACCENT__', ACCENT).replace('__OUTLIER__', OUTLIER)
	with open(path, 'w') as f:
		f.write(html)


_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Protein memory: transcription vs protein absence</title>
<style>
  :root { --accent:__ACCENT__; --outlier:__OUTLIER__; --ink:#1b2530; --muted:#5b6672; }
  body { font-family:'DejaVu Sans',-apple-system,Arial,sans-serif; color:var(--ink);
    margin:0; padding:24px; background:#fff; }
  h1 { font-size:18px; margin:0 0 2px; }
  p.sub { color:var(--muted); font-size:13px; margin:0 0 14px; }
  .controls { font-size:13px; margin-bottom:8px; color:var(--muted); }
  select { font-size:13px; padding:2px 4px; }
  #wrap { position:relative; width:640px; max-width:100%; }
  svg { overflow:visible; }
  .dot { cursor:pointer; }
  .tip { position:absolute; pointer-events:none; background:#fff; border:1px solid #ccc;
    border-radius:6px; padding:8px 10px; font-size:12px; box-shadow:0 2px 8px rgba(0,0,0,.15);
    opacity:0; transition:opacity .08s; max-width:240px; z-index:5; }
  .tip b { color:var(--ink); } .tip .k { color:var(--muted); }
  .legend { font-size:12px; color:var(--muted); margin-top:6px; }
  .sw { display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:4px;
    vertical-align:middle; }
</style></head>
<body>
<h1>Sub-generational transcription vs. protein absence</h1>
<p class="sub">Each point is a sub-generational gene (n=<span id="ng"></span>). Below the dashed
  y=x line, protein is present more often than it is transcribed &mdash; expression "memory."</p>
<div class="controls">Color by:
  <select id="colorby">
    <option value="none">none (highlight anti-memory outliers)</option>
    <option value="h">protein half-life</option>
    <option value="t">translation efficiency</option>
  </select>
</div>
<div id="wrap"><svg id="plot" width="640" height="620"></svg><div class="tip" id="tip"></div></div>
<div class="legend" id="legend"></div>
<script>
const DATA = __DATA__;
const ACCENT="__ACCENT__", OUTLIER="__OUTLIER__";
document.getElementById('ng').textContent = DATA.length;
const svg=document.getElementById('plot'), tip=document.getElementById('tip'),
      wrap=document.getElementById('wrap'), NS='http://www.w3.org/2000/svg';
const W=640,H=620,m={l:70,r:20,t:20,b:66}, pw=W-m.l-m.r, ph=H-m.t-m.b;
const sx=v=>m.l+v*pw, sy=v=>m.t+(1-v)*ph;
function el(t,a){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}
// axes box + gridlines
for(let g=0;g<=1.0001;g+=0.2){
  svg.appendChild(el('line',{x1:sx(g),y1:sy(0),x2:sx(g),y2:sy(1),stroke:'#eef1f4'}));
  svg.appendChild(el('line',{x1:sx(0),y1:sy(g),x2:sx(1),y2:sy(g),stroke:'#eef1f4'}));
  const tx=el('text',{x:sx(g),y:sy(0)+18,'text-anchor':'middle',fill:'#5b6672','font-size':11});
  tx.textContent=g.toFixed(1); svg.appendChild(tx);
  const ty=el('text',{x:m.l-8,y:sy(g)+4,'text-anchor':'end',fill:'#5b6672','font-size':11});
  ty.textContent=g.toFixed(1); svg.appendChild(ty);
}
svg.appendChild(el('line',{x1:sx(0),y1:sy(0),x2:sx(0),y2:sy(1),stroke:'#999'}));
svg.appendChild(el('line',{x1:sx(0),y1:sy(0),x2:sx(1),y2:sy(0),stroke:'#999'}));
svg.appendChild(el('line',{x1:sx(0),y1:sy(0),x2:sx(1),y2:sy(1),stroke:'#5b6672',
  'stroke-dasharray':'5 4'}));
// axis labels
let xl=el('text',{x:m.l+pw/2,y:H-8,'text-anchor':'middle',fill:'#1b2530','font-size':13});
xl.textContent='Fraction of cells with no transcription event'; svg.appendChild(xl);
let yl=el('text',{x:16,y:m.t+ph/2,'text-anchor':'middle',fill:'#1b2530','font-size':13,
  transform:`rotate(-90 16 ${m.t+ph/2})`});
yl.textContent='Fraction of cell-cycle time with absent protein'; svg.appendChild(yl);
// sequential ramp (light->dark accent) for color-by
function hex2rgb(h){return [parseInt(h.slice(1,3),16),parseInt(h.slice(3,5),16),parseInt(h.slice(5,7),16)];}
const A=hex2rgb(ACCENT);
function ramp(f){f=Math.max(0,Math.min(1,f));const r=Math.round(232+(A[0]-232)*f),
  g=Math.round(238+(A[1]-238)*f),b=Math.round(244+(A[2]-244)*f);return `rgb(${r},${g},${b})`;}
function logspan(key){const v=DATA.map(d=>d[key]).filter(x=>x!=null&&x>0).map(Math.log10);
  return [Math.min(...v),Math.max(...v)];}
const dots=[];
DATA.forEach((d,i)=>{
  const above=d.y>d.x;
  const c=el(above?'rect':'circle', above
    ? {x:sx(d.x)-4,y:sy(d.y)-4,width:8,height:8,transform:`rotate(45 ${sx(d.x)} ${sy(d.y)})`}
    : {cx:sx(d.x),cy:sy(d.y),r:3.4});
  c.setAttribute('class','dot');
  c.setAttribute('fill', above?OUTLIER:ACCENT);
  c.setAttribute('fill-opacity', above?0.95:0.42);
  if(above){c.setAttribute('stroke','#8a5a00');c.setAttribute('stroke-width',0.6);}
  c.addEventListener('mousemove',ev=>showTip(ev,d,above));
  c.addEventListener('mouseleave',()=>{tip.style.opacity=0;});
  svg.appendChild(c); dots.push({el:c,d,above});
  if(above){ // label anti-memory outliers
    const t=el('text',{x:sx(d.x)+7,y:sy(d.y)+3,fill:'#8a5a00','font-size':10});
    t.textContent=d.n; svg.appendChild(t);
  }
});
function showTip(ev,d,above){
  tip.innerHTML=`<b>${d.n}</b> <span class="k">(${d.g})</span><br>`+
   `<span class="k">no-transcription cells:</span> ${d.x.toFixed(3)}<br>`+
   `<span class="k">protein absent (time):</span> ${d.y.toFixed(3)}<br>`+
   `<span class="k">def-5 rate:</span> ${d.r} transcripts/gen<br>`+
   `<span class="k">translation efficiency:</span> ${d.t}<br>`+
   `<span class="k">protein half-life:</span> ${d.h==null?'&gt; cell cycle':d.h+' min'}`+
   (above?'<br><b style="color:#8a5a00">above line: anti-memory</b>':'');
  const rect=wrap.getBoundingClientRect();
  let px=ev.clientX-rect.left+12, py=ev.clientY-rect.top+12;
  if(px>rect.width-180)px=ev.clientX-rect.left-180;
  tip.style.left=px+'px'; tip.style.top=py+'px'; tip.style.opacity=1;
}
const legend=document.getElementById('legend');
function setColor(mode){
  let span=null;
  if(mode!=='none') span=logspan(mode);
  dots.forEach(o=>{
    if(o.above) return;                       // outliers keep their color+shape
    if(mode==='none'){o.el.setAttribute('fill',ACCENT);o.el.setAttribute('fill-opacity',0.42);}
    else{const v=o.d[mode];
      if(v==null||v<=0){o.el.setAttribute('fill','#ccc');}
      else{o.el.setAttribute('fill',ramp((Math.log10(v)-span[0])/(span[1]-span[0])));}
      o.el.setAttribute('fill-opacity',0.85);}
  });
  legend.innerHTML = mode==='none'
    ? `<span class="sw" style="background:${ACCENT}"></span>subgen gene`+
      `&nbsp;&nbsp;<span class="sw" style="background:${OUTLIER};border-radius:0;transform:rotate(45deg)"></span>anti-memory (protein absent more than transcript)`
    : (mode==='h'?'protein half-life':'translation efficiency')+
      ` (light = low, dark = high, log scale); &#9670; = anti-memory outlier`;
}
document.getElementById('colorby').addEventListener('change',e=>setColor(e.target.value));
setColor('none');
</script>
</body></html>
"""


# ------------------------------------------------------------------ main

def main():
	ap = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument('pergene_memory', help='subgen_protein_memory_pergene.tsv')
	ap.add_argument('sim_data', help='kb/simData.cPickle')
	ap.add_argument('output_dir')
	args = ap.parse_args()
	os.makedirs(args.output_dir, exist_ok=True)
	prefix = os.path.join(args.output_dir, 'subgen_protein_memory_interactive')

	rows = load_memory(args.pergene_memory)
	with open(args.sim_data, 'rb') as f:
		sim_data = pickle.load(f)
	rows = annotate(rows, sim_data)
	print('annotated %d genes' % len(rows))

	# Enriched per-gene TSV.
	cols = ['gene_id', 'name', 'cistron_id', 'monomer_id', 'x', 'y', 'def5_rate',
		'trl_eff', 'deg_rate_per_s', 'half_life_min']
	with open(prefix + '_pergene.tsv', 'w') as f:
		w = csv.writer(f, delimiter='\t')
		w.writerow(cols)
		for r in rows:
			w.writerow([r['gene_id'], r['name'], r['cistron_id'], r['monomer_id'],
				'%.5g' % r['x'], '%.5g' % r['y'], '%.5g' % r['def5_rate'],
				'%.5g' % r['trl_eff'], '%.5g' % r['deg_rate'],
				('%.5g' % r['half_life_min']) if np.isfinite(r['half_life_min'])
					else 'inf'])

	res = analyze_outliers(rows)
	with open(prefix + '_outlier_analysis.json', 'w') as f:
		json.dump(res, f, indent=2)
	write_html(rows, prefix + '.html')

	print('\n== outlier analysis ==')
	print('genes above the y=x line (anti-memory): %d' % res['n_above_line'])
	print('memory_gap vs log10(half-life): r=%.3f (spearman %.3f)'
		% (res['corr_gap_vs_log_halflife']['pearson_r'],
			res['corr_gap_vs_log_halflife']['spearman_r']))
	print('memory_gap vs log10(trl_eff):   r=%.3f (spearman %.3f)'
		% (res['corr_gap_vs_log_trleff']['pearson_r'],
			res['corr_gap_vs_log_trleff']['spearman_r']))
	print('std regression gap ~ log_hl + log_te: beta_hl=%.3f beta_te=%.3f'
		% (res['std_regression_gap']['beta_log_halflife'],
			res['std_regression_gap']['beta_log_trleff']))
	print('population median half-life=%.1f min, trl_eff=%.2f'
		% (res['population_median_halflife_min'], res['population_median_trleff']))
	print('\nabove-line genes:')
	for g in res['above_line_genes']:
		print('  %-8s x=%.2f y=%.2f gap=%.2f  trl_eff=%.2f  half_life=%s min'
			% (g['name'], g['x'], g['y'], g['gap'], g['trl_eff'],
				g['half_life_min']))
	print('\nWrote HTML + enriched TSV + outlier JSON to %s' % args.output_dir)


if __name__ == '__main__':
	main()
