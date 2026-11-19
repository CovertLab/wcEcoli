# GFP burden paper: simulations and figures

Everything needed to reproduce the paper's simulations and figure panels.
Modelled on `runscripts/operon_paper`.

| File | What it is |
| --- | --- |
| `paper_runs.sh` | Every simulation batch, in submission order |
| `paper_analysis.sh` | Every paper analysis, per batch |
| `panels.tsv` | The manifest: panel label to script, batch and status |
| `batches.tsv` | The batch registry: id, size, paths, git and sim_data hashes |

Shared code lives outside this directory, so analyses can import it:

| Module | What it holds |
| --- | --- |
| `wholecell/analysis/gfp_paper_constants.py` | Generation window, grid definition, new gene id lookup |
| `wholecell/analysis/gfp_paper_style.py` | Colours, fonts, panel sizes, `save_panel()`, failed-lineage rule, stats |

## Conventions

**Panel labels are not script names.** The storyboard may still be reordered,
so each script is named for what it plots and the `Figure_3F_` prefix is
attached at save time from `panels.tsv`. Reordering the paper is then an edit
to one file instead of a round of renaming. Each script also carries a
`# Figure 3F` comment, and `panels.tsv` is the source of truth if the two
disagree.

**Every panel is saved three ways plus metadata.** `save_panel()` writes PNG,
PDF and SVG with editable text, and a JSON sidecar naming the git commit,
sim_data hash, batch, variants, seeds, generation window and sample size. Any
figure can be traced back to the code and data that made it.

**Everything lands in one folder per batch**, `<batch>/plotOut/paper_figures/`,
including output from multigen analyses, which would otherwise scatter panels
across generation directories. One download gets the whole set. Panels with
one file per variant or seed go in a subfolder named after the panel.

**Averages use generations 16-23**, the last eight. Generations 0-7 are
wildtype steady state and 8-15 are the induction transient. Take the window
from `gfp_paper_constants`, never a local constant.

**Designs whose lineages died stay in the figures, drawn hollow**, via
`gfp_paper_style.design_scatter`. The survivor count goes in the panel
metadata, because a hollow point averages only the lineages that survived.

**Growth is reported as growth rate relative to the seed-matched no-GFP
control**, not doubling time, except where a formula needs the period
(Cooper-Helmstetter).

## Running

Simulations go to the `owners` partition, pinned to CPU generations with
AVX-512 so results are reproducible across nodes; analyses go to `mcovert`.
**Only one simulation batch runs at a time**, so `paper_runs.sh` is a sequence,
not a parallel launch. Analyses may run on `mcovert` while a batch is going.

Submission order, and why each is where it is:

1. `smoke` - proves the setup before 6,000 sims depend on it
2. `B1` scan - everything downstream waits on it
3. `W` wildtype - independent, and fills the slot while B1 is analysed
4. `B2` focus - needs the design chosen from B1
5. `B1deep` - needs the subset chosen from B1
6. `M1` - needs the kcat table from W
7. `M2` - needs the kcat table and the below-the-line set from B2

Record each batch in `batches.tsv` as it is submitted and again when it
finishes, including the resolved CPU generation. That file is also the
manifest for the eventual data deposit.

## Status

Skeleton only. `paper_runs.sh` and `paper_analysis.sh` are filled in once the
ParCa and smoke test have passed; `panels.tsv` lists all 44 panels of the
locked 2026.11.13 storyboard with their scripts still to be written.
