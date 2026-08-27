# TissueMosaic on the JHH IMC cohort

Trains a self-supervised model (DINO) on the spatial arrangement of cell types in
the JHH IMC PDAC ROIs, then embeds patches of tissue so micro-environments can be
compared across ROIs and patients.

Run everything from the **repo root** with the `tissue-mosaic-260819` conda env.
The upstream TissueMosaic checkout is expected at `~/TissueMosaic`; override with
`TISSUEMOSAIC_HOME`.

## Layout

```
tissuemosaic/
├── imc_tm.py                  all shared logic; notebooks and scripts import from here
├── tm_compat.py               the one compatibility shim, applied on import of imc_tm
├── config_dino_ssl_imc.yaml   DINO config for this cohort
├── make_imc_anndata.py        headless CLI around imc_tm.build_anndata()
├── analyze_imc_features.py    featurize every ROI, PCA + UMAP, sanity figure
├── run_*.py                   training / featurisation drivers -- see Pipeline order
├── ecm_*.csv                  small result tables, committed
├── notebooks/                 one per experiment
│   └── archive/               superseded, outputs stripped
├── imc_anndata*/              generated, ignored -- see Generated directories
└── runs/                      training output, ignored -- ~84 GB
```

### Drivers

| script | role |
|--------|------|
| `run_train_imc.py` | shim + handoff to upstream `main_1_train_ssl.py` |
| `run_featurize_imc.py` | shim + handoff to upstream `main_2_featurize.py` |
| `run_featurize_screen.py` | featurize a whole screen and record `_provenance.json` |
| `run_window_sweep.py` | trains the six S1 models (long-running, run from a terminal) |
| `run_final_sweep.py` | trains the two final configs, 8 ROIs / 1000 epochs (`--parallel`) |
| `run_ecm_screen.py` | the 11-channel runs: `--parallel`, `--long` (1000 epochs), `--control` |

### Notebooks

| notebook | what it answers |
|----------|-----------------|
| `0.prepare_and_check_data` | CSV -> h5ad, coordinate-frame proof, cohort + JHH369 sanity checks |
| `S1.window_size_sweep` | six-way window size / resolution comparison: UMAP, Leiden, metrics |
| `S1e.compare_A_vs_E` | A (128 um) vs E (64 um) head to head, full downstream analysis |
| `S1f.compare_fov_series` | the 2.0 um/px FOV series 64 -> 96 -> 128 um: smooth climb or a step? |
| `S2.final_sweep` | final sweep evaluation: both filters, per-ROI gain, three-tier transfer |
| `S3.qualitative` | what the two final models produce, read by eye rather than by score |
| `S4.ecm_channel` | does the 11th channel help? matched-epoch vs the no-ECM arm |
| `S4b.ecm_qualitative` | the S3 treatment applied to the 200-epoch ECM runs, 96/128/192 um |
| `S4c.fov_three_way` | window size inside the ECM arm; why raw kappa misleads across FOV |
| `S4d.ecm_qualitative_ep499` | the 1000-epoch runs at epoch 499, plus the arrangement-vs-composition test |

Every notebook resolves `imc_tm` by walking up from the working directory, so it runs
from the repo root, from `tissuemosaic/`, or from `notebooks/`.

### notebooks/archive/

Kept for provenance, **outputs stripped** -- 11 MB of stored figures down to 152 KB.
Nothing here is maintained, and paths in them may reference run directories that no
longer exist.

| notebook | why archived |
|----------|--------------|
| `S1a.inspect_one_run` | single-config inspection; superseded by S1's six-way comparison |
| `S1a.inspect_one_run copy` | editor duplicate of the above |
| `S1b.inspect_opt1_P4` | training-parameter option P4; settled by S1 and S2 |
| `S1c.inspect_opt2_P5` | training-parameter option P5; settled by S1 and S2 |
| `S1d.inspect_opt3_v1proxy` | v1-proxy option; settled by S1 and S2 |
| `S3.qualitative copy` | editor duplicate of `S3.qualitative` |
| `notebook1` | exploratory scratch from before the pipeline existed; no reproducible path |

### Do not add an `__init__.py` to this directory

This directory is named `tissuemosaic`, the same as the installed TissueMosaic
library. Without an `__init__.py` it is only a *namespace-package candidate*, which
loses to the regular package in site-packages -- so `import tissuemosaic` correctly
resolves to the library. **Adding an `__init__.py` makes this directory win instead**,
and `from tissuemosaic.data import AnndataFolderDM` then fails with
`ModuleNotFoundError`. Import these modules flat:

```python
import sys; sys.path.insert(0, "<repo>/tissuemosaic")
import imc_tm
```

For the same reason, keep module names distinctive (`imc_tm`, not `utils`) -- that
`sys.path.insert` puts this directory *first*.

Run notebooks with the `tissue-mosaic-260819` kernel, created from
`conda_env_configs/tissuemosaic_260818.yaml` (the filename dates the file, not the env;
its `name:` field is `tissue-mosaic-260819` and its pins are the ones that match). Notebook 0 **deletes and
regenerates** `imc_anndata/`, so it is safe to re-run and is the canonical way to
rebuild the inputs.

## S1: window size and rasterisation resolution

Window FOV = `pixel_size x global_size`, so the two knobs are not independent --
sweeping one without pinning the other confounds scale with resolution. Six configs
(`imc_tm.WINDOW_SWEEP`), added in two rounds:

| | FOV | pixel_size | global/local | measured train |
|---|---|---|---|---|
| E | 64 um | 2.0 | 32 / 24 | 1.4 h |
| F | 96 um | 2.0 | 48 / 32 | 1.5 h |
| A | 128 um | 2.0 | 64 / 48 | 1.6 h |
| B | 192 um | 2.0 | 96 / 64 | 2.1 h | <- original production config |
| C | 192 um | 3.0 | 64 / 48 | 1.6 h |
| D | 192 um | 1.5 | 128 / 80 | 2.6 h |

Round 1 was A-D, an L hinged on B. It found 128 um beat every 192 um variant -- and
128 um was the *smallest* window tested, so the trend pointed below the sampled range.
Round 2 added E and F to bracket it, making 128 um an interior point rather than an
endpoint. 64 um is the practical floor at 2.0 um/px: the backbone downsamples 16x, so
gs 32 -> 2x2 and local 24 -> 2x2 pre-pool, while an input of 16 would collapse to 1x1.

- **E vs F vs A vs B** -- FOV at fixed resolution, a 3x range. *How much tissue?*
- **C vs B vs D** -- resolution at fixed FOV. *How finely rasterised?*
- **A vs C** -- identical network input (64/48) and compute, differing only in physical
  scale, so FOV is isolated from model capacity.

```bash
python tissuemosaic/run_window_sweep.py --dry-run       # plan only
python tissuemosaic/run_window_sweep.py                 # ~10.8 h, JHH369 only
python tissuemosaic/run_window_sweep.py --skip-trained  # resume
```

Then run `S1.window_size_sweep.ipynb`.

### Two filters, and why both must scale

`n_element_min_for_crop` (10) is the **training** floor, and `cropper_test` reuses it
(`datamodule.py:331`). On top of it the analysis applies `imc_tm.fov_scaled_threshold`,
= 0.5 x the cells a window of that FOV is expected to hold at JHH369 density, giving
**13 / 29 / 52 / 117** cells at 64 / 96 / 128 / 192 um.

A *fixed* threshold is itself FOV-dependent -- the exact bias it exists to remove.
Round 1 used a flat 50, which is 48% of an expected 128 um window but only 21% of a
192 um one, so it stripped sparse border tiles at small FOV and spared them at large
(retained border share 32% at 128 um vs 51% at 192). Under the scaled rule the six sit
at 13.6 / 26.8 / 15.3 / 28.2 / 16.9 / 19.4 % (A-F), and retained *coverage* equalises
to 26-28k cells for every config, an 8% spread.

### Categorical colours

**Every categorical variable is coloured from `colorcet.glasbey_dark`**, each family on
its own slice so that a cluster, a cell type and an ROI never share a colour -- several
figures put those legends side by side:

| family | helper | slice | min dLab | replaced |
|---|---|---|---|---|
| clusters | `imc_tm.cluster_palette()` | `[0:n]` | 23.7 | tab20, 16.6 |
| cell types | `imc_tm.palette()` | `[16:26]` | 37.5 | tab10, 27.7 |
| ROIs | `imc_tm.roi_palette()` | `[40:44]` | 49.3 | tab10, 40.3 |

Offsets live in `imc_tm.GLASBEY_OFFSET`. They cost nothing: the cell-type slice at 16
measures *better* than glasbey's own first ten (37.5 vs 36.2), so cross-family
separation is free.

Why it mattered: Leiden gives 6-15 clusters here and tab20 is seven light/dark hue
PAIRS, so 15 clusters read as roughly 8 colours and neighbouring clusters were
indistinguishable in the tissue maps. Glasbey picks each successive colour to be as far
as possible from every colour already chosen, so distinctness degrades gracefully.

The cluster row above is measured in CIELAB, which *understates* colorcet's palette --
it optimises in CAM02-UCS. Both beat tab20 on the same yardstick, which is the
comparison that matters. `cluster_palette_source()` reports which palette is live, and
there is a CIELAB Glasbey fallback if colorcet is ever absent.

### Tissue maps use data coordinates

`imc_tm.plot_window_map()` draws each window as a real rectangle of its physical size
rather than a points-sized scatter marker. A points-sized square only tiles the tissue
if you guess a size matching the figure geometry; guess low and the map becomes a dot
grid with white gutters, which hides the one thing the panel exists to show. With data
coordinates a 128 um window is 128 um wide at any figure size, zoom or dpi, so adjacent
windows always share an edge and the FOV rows are directly comparable. Remaining white
space is then meaningful: those are windows the FOV-scaled filter dropped.

### Rank on kappa, not on the z-score

`spatial_coherence` returns both. `coherence_kappa` = (obs - null) / (1 - null) is a
probability difference and does not grow with N. The z-score divides by a null SD that
shrinks as ~1/sqrt(n_pairs), so **z grows as ~sqrt(N) at a fixed effect size**.

Equalising coverage does not equalise count -- a smaller window tiles the same tissue
more times, so the six configs retain 118 to 1078 windows, a 9x spread and a 3x swing
in z from arithmetic alone. Ranking on z would hand the win to the smallest FOV
automatically. Use z only to ask whether a config beats chance at all. (Within a
single geometry, as in the P1-P5 training probes, N is matched and z is fine.)

Kappa is not the end of it either. `imc_tm.coherence_ceiling` scores synthetic
ground-truth domains of a fixed physical size on each config's real window grid, and
a PERFECT model scores 3.9x higher at 64 um than at 192 um for 250 um domains -- a
192 um window is roughly one domain, so neighbours rarely share one, while sixteen
64 um windows tile the same domain. Report `kappa / ceiling`, at several assumed
domain sizes, and check the ranking is stable across them.

### Schedule

`max_epochs` 1000 -> 200 is the deliberate screen budget. `warm_up_epochs` 100 -> 20
with `warm_down_epochs` at 100 preserves TissueMosaic's schedule *shape*: the LR
milestones are `(0, warm_up, max_epochs - warm_down, max_epochs)` = (0, 20, 100, 200),
a short warm-up then a 50% cosine decay. Note that leaving both at 100 would collapse
the peak-LR plateau to zero width **without erroring**, because the assert is only
`warm_up + warm_down <= max_epochs`, i.e. `200 <= 200`.

That budget is 6,400 steps per config: a **comparative screen**, not converged models.
The ranking is fair because every config gets an identical schedule; retrain the
winner properly, and confirm on the full cohort, before it becomes production.

`global_size` / `local_size` are multiples of 16. The backbone downsamples 16x, giving
clean 4x4 / 6x6 / 8x8 pre-pool maps. Other sizes do not error (it ends in
`AdaptiveAvgPool2d(1)`) but downsample raggedly.

**Caching hazard:** `AnndataFolderDM` writes `train_dataset.pt` / `test_dataset.pt`
next to the h5ad, and those rasters are built for ONE `pixel_size`. Reusing them
across sweep configs would silently train on the wrong resolution. `run_window_sweep.py`
deletes them before every run.

## S2: the final sweep

Two configs, 8 ROIs, 1000 epochs. `run_final_sweep.py --parallel`.

```
                      final_fov96      final_fov128
pixel_size                    2.0              2.0
global / local            48 / 32          64 / 48
FOV                         96 um           128 um
max_epochs                   1000             1000
warm_up / warm_down       100 / 100        100 / 100
n_element_min_for_crop         10               10     <- train AND inference
n_crops_for_tissue_train     1024             1024
weight decay              0.0 / 0.0        0.0 / 0.0
batch_size_per_gpu            128              128
```

`warm_up`/`warm_down` 100/100 at 1000 epochs is **upstream's shape verbatim** (their
config, not their argparse defaults). The S1 screen used 20/100 instead, because
100/100 at `max_epochs` 200 collapses the peak-LR plateau to zero width.

### What this run is actually for

Not 96-vs-128 -- the screen already settled that, with 128 um gaining 3.03x over a
cell-count vector against 96 um's 1.23x. The open question is **transfer**: a model
trained on JHH369 alone scored gain **0.33-0.78x on unseen patients**, i.e. WORSE than
clustering cell counts. Five patients in training is the test of whether that was a
data-diversity problem or something the objective cannot fix.

### The 8 ROIs

`imc_tm.FINAL_ROIS`. The four JHH369 ROIs for continuity, plus two high-quality and two
low-quality, where the low pair is low on **orthogonal** axes so the failure modes stay
separable:

| ROI | density | biology % | role |
|---|---|---|---|
| JHH387ROI1 | 10,590 | 48.4% | high: dense and signal-rich |
| JHH380ROI1 | 8,262 | 62.1% | high: richest classified signal |
| JHH357ROI5 | 7,396 | **2.2%** | low: dense but signal-poor |
| JHH317ROI2 | **2,643** | 43.7% | low: sparse but signal-rich |

61,103 cells, 5 patients, statuses 0/2/6/8/10. Density spread 4.0x, biology 28x.

### Three tiers of generalisation

Four training patients contribute other ROIs to the cohort, so the 40 untrained ROIs
split into two different tests rather than one:

| tier | ROIs | patients | isolates |
|---|---|---|---|
| training | 8 | 5 | in-sample fit |
| held-out ROI, seen patient | 11 | 4 | a new section of a known patient |
| held-out ROI, unseen patient | 29 | 9 | a genuinely new patient |

Gain surviving tier 2 but dying at tier 3 is a patient-level failure, not section-level.

### Two operational gotchas

**One h5ad folder per config.** `AnndataFolderDM` caches its raster to
`<data_folder>/train_dataset.pt`, and both configs share `pixel_size 2.0`, so a shared
folder has them racing to write the same file under `--parallel`. `final_h5ad_dir()`
gives each its own (~38 MB).

**Launch the local `run_train_imc.py` shim, never upstream's `main_1_train_ssl.py`
directly.** `tm_compat` is applied by importing `imc_tm`, and a subprocess inherits none
of the parent's imports -- launching upstream directly leaves `weights_only=True` and
TissueMosaic dies reading back its own `train_dataset.pt` (`SparseImage` is not an
allowed global). This cost one failed launch.

A corollary that has misled me twice: because the shim runs upstream via
`runpy.run_path`, **`main_1_train_ssl.py` never appears in the process argv**. Grep for
`run_train_imc` when checking whether training is alive, or for `main_1` and conclude
the run has died when it is fine. The same trap applies to any kill pattern.

### Runtime

8 ROIs x 1024 crops / batch 128 = **64 steps/epoch**, 64,000 steps per config.

Measured ~26 h each, run concurrently, not the ~17 h first estimated. The error was
scaling *steps* with ROI count while assuming the single-config per-step rate would
hold: instantaneous rate does hold (~1.05 it/s), but roughly **30 s/epoch of
non-step overhead** -- the cropper re-drawing 1024 crops per tissue every epoch --
also doubled with the ROI count, and it is single-threaded so it neither overlaps GPU
compute nor parallelises. Sequential would be ~52 h; parallel contention is minimal
(GPU ~57%, 5.9 GB of 16).

Evaluate with `S2.final_sweep.ipynb`, which accepts an `EPOCH` parameter and runs
against periodic checkpoints mid-flight. Note the cosine decay runs epochs 900-1000, so
anything read before 900 has not consolidated.

## S4: the ECM channel

`tissuemosaic/runs/ecm_screen/` — two 200-epoch runs, `ecm_fov96` and `ecm_fov128`, identical to the
final sweep except that the raster carries an **11th channel holding the PhysiCell ECM
initial condition**. Built by `imc_tm.build_ecm_anndata`, driven by
`run_ecm_screen.py`, evaluated by `S4.ecm_channel.ipynb`.

### Why ECM needs its own rows but not its own image

ECM is a BioFVM *continuum substrate*, not an agent: a scalar field on a Cartesian voxel
mesh at a **20 µm lattice**, seeded from `<ROI>_ecm_scaled.csv`. Cells are points at
measured positions. Different supports, so they enter the AnnData as different rows —
cell rows are one-hot in channels 0-9 with 0 in channel 10; ECM rows are 0 in 0-9 and
carry `ecm / 10` in channel 10.

Only the *rows* are separate. After rasterisation there is one image of shape
`(11, W, H)`; the model never sees two objects. Splitting the rows is what lets matrix
occupy cell-free space — folding the field onto cell rows instead would discard it
precisely where no cells were segmented, which in a sparse ROI like JHH317ROI2 is where
matrix is highest.

`ECM_DIVISOR = 10` because the source CSV is already normalised **per ROI to its own
95th percentile**. Channel 10 is a within-section relative density with no cross-section
absolute scale.

### The one-hot branch changes what `n_element_min_for_crop` means

Every ECM value is `< 1`, which trips `CropperSparseTensor`'s
`if torch.any(values < 1.0)` branch (`dataset.py:365`). That branch counts **nonzero
entries** instead of summing cell mass. Matrix alone contributes up to 23.5 entries per
96 µm window and 41.3 per 128 µm window, so `n_element_min_for_crop = 20` is satisfied
by matrix by itself at 128 µm and needs only a handful of cells at 96 µm.

**The ECM arm therefore trained on a larger, emptier tile set than the 10-channel arm at
threshold 10.** That is a second difference between the arms beyond the extra channel.
The original code comment estimated ~9 and ~20 entries and was wrong by ~2x; it has been
corrected in place.

### Comparing arms honestly: the shared grid

The tiling origin is `torch.randint(low=-crop_size, high=0)` seeded per ROI and the
stride follows from `global_size` alone, so **the candidate window positions are
identical between an ECM run and its no-ECM counterpart** regardless of channel count.
Only the `n_element_min` cull differs. So featurise both with
`n_element_min_for_crop=1`, join on `imc_tm.window_key` (ROI + centre to 1 nm) via
`imc_tm.align_to_grid`, then apply one filter to both. Measured: the base grid is a
strict subset of the ECM grid — 1161 of 1171 at 96 µm, 700 of 704 at 128 µm, the extras
being tiles holding matrix but zero cells — and centres and cell counts agree exactly.

### SmartPca is unseeded upstream, and it made the baseline irreproducible

Found while building S4. Rebuilding the graph on identical input gave bit-identical kappa
for the 512-dimensional DINO features but drew anywhere in 0.147-0.160 for the 10-channel
composition baseline — a ±4% wobble on the denominator of every gain, in a comparison
whose smaller effect is +8%.

The cause is `SmartPca.fit` in `tissuemosaic/utils/validation_util.py`, which regularises
the covariance matrix with an **unseeded** draw before taking its SVD:

```python
eps = 1E-4 * torch.randn(p, dtype=cov.dtype, device=cov.device)
cov += torch.diag(eps)
```

Negligible against 512 dimensions of DINO features. Against a 10-dimensional composition
vector — z-scored, so the covariance diagonal is ~1, with many near-duplicate windows — it
perturbs the retained subspace enough to flip kNN ties, and Leiden moves with them.

`imc_tm.umap_graph` and `imc_tm.cluster_embedding_tm` now call `torch.manual_seed(seed)`
before `SmartPca`. Verified reproducible across repeats, call orders and processes; S4
asserts it over all 72 (fov, arm, resolution, seed) settings. **Numbers produced before
this fix carry that variance**: any notebook that clusters through `cluster_embedding_tm`
(S2, S3) has a composition baseline good to roughly ±4%, so gains quoted there can move by
about that much on re-run. The DINO arms were never affected — the jitter is negligible at
512 dimensions — and ratios between two DINO arms are immune either way, because the
baseline cancels. `cluster_embedding` (the scanpy path) goes through `sc.pp.pca` and was
never affected. Worth reporting upstream alongside the `weights_only` issue.

### Results, epoch 199, 8 ROIs, gain over the 10-channel composition baseline

| arm | 96 µm | 128 µm |
|---|---|---|
| composition, 10 ch (reference) | 1.00x | 1.00x |
| composition, 11 ch (+ECM) | 1.09x | 1.25x |
| DINO 10 ch @199 | 1.27x | 2.37x |
| **DINO 11 ch (+ECM) @199** | **1.81x** | **2.56x** |
| DINO 10 ch, best of 8 checkpoints | 1.49x (ep 499) | 3.06x (ep 599) |

**At 96 µm the channel helps and the result is robust.** ECM beats the 10-channel arm at
every Leiden resolution 0.6-1.5 (by 1.23-1.43x) and beats the base model's best checkpoint,
not just the matched one.

**At 128 µm it is not established.** The +8% over `base@199` shrinks to +1% at resolution
0.8, and the base model reaches 3.06x by epoch 599 — above anything the 200-epoch ECM run
produced. Whether an ECM run trained as long would exceed it is untested.

The 11-channel *composition* baseline (1.09x / 1.25x) shows some matrix signal is trivially
recoverable without a model, but the 11-channel DINO model clears it by 1.66x at 96 µm and
2.06x at 128 µm, so it is not just reading the field off.

### What is missing

**No matched no-ECM control was trained** — the design was two ECM runs only. The
comparison borrows the final sweep's checkpoints, which share ROIs, pixel size, crop
geometry, batch size and weight decay but *not* the LR protocol: `(0, 100, 900, 1000)`
against `(0, 20, 180, 200)`, so the base is at peak LR at epoch 199 while the ECM runs
have annealed. **That confound favours ECM**, which is why `S4` compares against the base
model's whole trajectory rather than one epoch. A 200-epoch 10-channel run on the
identical schedule is the experiment that would settle it:
`python tissuemosaic/run_ecm_screen.py --control --parallel`, about 5.5 h.

Nothing in S4 is a transfer test: eight ROIs, five patients, all training data.

### Runtime

64 steps/epoch, 200 epochs, both concurrent: **5.6 h** wall clock (20:49 to 02:28), about
2.5 min/epoch early and faster once the crop cache warmed. Roughly 22% slower per epoch
than the 10-channel sweep at the same step count — the extra channel plus ~22.6k
additional sparse entries per pass.

### Three quantities get called "ECM content", and only one of them means it

This bit us twice. Channel 10 is a continuous field on a 20 um lattice, and "how much
matrix is in this window" has three plausible answers that behave very differently.

| | definition | what it ALSO tracks |
|---|---|---|
| `n_ecm` / `ecm_mass` | sum of channel-10 values in the window | **window occupancy** -- r = +0.77 with the voxel count alone, and it scales with window AREA |
| `frac_ecm` | `n_ecm / total window mass` | **cell density** -- it is a ratio against cell mass |
| `ncv_ecm` | `patch_ncv[:, 10]`, an ENTRY COUNT fraction | almost purely inverse cell density, r = -0.62 with cell mass |
| **`ecm_density`** | `n_ecm / n_ecm_voxels` (`imc_tm.ecm_density`) | nothing -- geometry divided out |

Measured medians, which show which ones are area-invariant:

| | 96 um | 128 um | 192 um |
|---|---|---|---|
| `ecm_mass` | 12.26 | 20.15 | 38.57 |
| `ecm_voxels` | 20 | 36 | 81 |
| **`ecm_density`** | **0.580** | **0.578** | **0.582** |
| `frac_ecm` | 0.174 | 0.172 | 0.175 |

`ecm_mass` runs 1.00 / 1.64 / 3.15 against window areas of 1.00 / 1.78 / 4.00, so it must
never be compared across FOV. `ecm_density` is flat, as a density should be.

**The control that settles it is `n_ecm_voxels`** -- how many lattice sites the window
covers, containing no matrix values whatsoever. Clustering eta^2 against each measure:

| | `ecm_density` | `frac_ecm` | `ecm_voxels` (control) |
|---|---|---|---|
| ecm_fov96 | 0.430 | 0.683 | 0.540 |
| ecm_fov128 | 0.322 | 0.661 | **0.670** |
| ecm_fov192 | 0.323 | 0.667 | 0.515 |

At 128 um the pure-geometry control scores HIGHER than the real measure. Any eta^2 quoted
on mass or on fraction is substantially "the model resolved where the tissue edge is".

**S4's headline was affected.** It originally reported the base -> ECM jump on `n_ecm`:
0.153 -> 0.546 at 96 um and 0.411 -> 0.539 at 128 um. On `ecm_density` the same comparison
is **0.145 -> 0.381 (2.6x)** and **0.257 -> 0.388 (1.5x)**. The conclusion -- the
11-channel model genuinely separates on matrix -- survives; the absolute numbers were
inflated. At 128 um the correction makes the ECM effect relatively LARGER, because the
10-channel model's 0.411 was mostly occupancy. S4 and S4b now plot all three measures with
the voxel-count control beside them.

Note S4c is unaffected: the FOV comparison never quantifies ECM at all -- its baseline is
the ten cell channels with matrix dropped.

### S4c: window size inside the ECM arm (96 vs 128 vs 192 um)

`ecm_fov192` was added to test S4's reading -- that ECM earns its +43% at 96 um by
supplying large-scale context a small window lacks. If so, a bigger window needs matrix
less and the optimum should not move upward. Geometry is `B_fov192_px2.0`'s verbatim
(global 96 / local 64), so the arm stays comparable to S1's round-1 192 um result.
Evaluated by `S4c.fov_three_way.ipynb`.

**Raw kappa is not comparable across FOV; gain is.** S4 compared two models at the SAME
fov on a shared grid, where raw kappa is fine. Across FOVs the grids differ by construction
(~1050 / ~640 / ~330 windows of different physical size) and bigger windows are more like
their neighbours whether or not the model improved. MEASURED: **192 um has the highest raw
kappa of any config in the screen (0.443) and only the second-best gain**, because its
composition baseline rose just as much (0.214 against 0.148 at 128 um). Scoring each config
against a baseline clustered on its own windows cancels the geometry.

The baseline is the ten CELL channels, renormalised, with ECM dropped -- same quantity at
every FOV, and it does not let matrix into the reference. Verified: dropping channel 10 and
renormalising reproduces the 10-channel run's own composition to **0.000000**.

Two filters, since `n_cells >= 10` admits a 4x emptier window at 192 um than at 96 um:
`flat` (10/10/10, S4's operating point) and `scaled` (10/18/40, constant minimum density).

| arm | 96 um | 128 um | 192 um |
|---|---|---|---|
| 10 channels, `flat` | 1.27x | 2.37x | *no run* |
| **11 channels, `flat`** | 1.81x | **2.56x** | 2.07x |
| 10 channels, `scaled` | 1.27x | 2.23x | *no run* |
| **11 channels, `scaled`** | 1.81x | **2.48x** | 2.17x |

**128 um wins under both filters, and the resolution bands for 128 and 192 do not overlap**
(128: 2.48-3.01 flat / 2.41-2.75 scaled; 192: 1.88-2.09 flat / 2.01-2.35 scaled). The
optimum did not move when the channel was added, and the FOV response stays an interior
optimum rather than a monotone climb -- the same shape S1 found without ECM.

192 um also resolves less: 8 clusters over 337 windows against 12 over 661 at 128 um, and
in S4b's tissue maps its windows are visibly too coarse to follow tissue structure.

**The ECM uplift shrinks with window size** -- +43% at 96 um, +8% (`flat`) to +11%
(`scaled`) at 128 um -- which is what the "matrix supplies missing context" reading
predicts, and means the channel is worth least exactly where you would operate. It cannot
be measured at 192 um: there is no `final_fov192`, and building one is only worthwhile if
192 um were competitive, which it is not.

Caveats: 192 um is the noisiest arm (286-325 windows after filtering against 595-1053);
`n_element_min_for_crop = 20` is fully vacuous there, since the ECM lattice alone puts ~92
voxels in a window and the cropper counts entries, so that arm trained on the emptiest tile
set of the three; and the S4 LR-schedule confound still inflates the ECM uplift, though it
does not touch the FOV comparison, which is ECM-vs-ECM throughout.

### S4b: the qualitative read

`S4b.ecm_qualitative.ipynb` — the S3 treatment applied to the two ECM models: merged
window-level and cell-level AnnData, UMAPs by cluster / by ROI / by ECM content, per-ROI
tissue maps, composition heatmaps, and the eight boxplot grids. Data comes from
`run_featurize_screen.py --screen ecm --epoch 199`, which drives upstream
`main_2_featurize` and writes `_provenance.json` so a notebook can refuse to compare two
configs featurized from different epochs.

Three things had to change for eleven channels:

- **`imc_tm.channels_for(n)`** returns the channel-name list matching an n-column block,
  and `merge_featurized` / `merge_cells` now use it. They previously zipped
  `len(CELL_TYPES)` names onto whatever width the data had, which on ECM data drops
  channel 10 from every plot *without erroring* — the worst possible failure mode. The
  cell-level `cell_type` categorical had the same bug in reverse: with `categories=
  CELL_TYPES` every ECM row silently becomes NaN.
- **`patch_ncv` counts entries, not mass.** `sparse_image.py:1026` builds it as
  `count_nonzero(channel == code) / total`. For a one-hot cell channel a count *is* a
  mass, so S1-S3 are unaffected. For a continuous channel it is not: on JHH369ROI1's
  densest 128 um window `patch_ncv` puts `ecm` at **0.205** while matrix is **0.017** of
  the window's actual mass. The count version is `n_ecm_voxels / (n_cells + n_ecm_voxels)`
  and, since ECM voxels sit at fixed 20 um lattice density, it is close to an inverse
  CELL-DENSITY map -- and those voxels being ~20% of the entries diluted every cell-type
  fraction with the same geometry term. `imc_tm.window_mass` computes a proper mass basis
  and S4b swaps `obsm['composition']` to it before any analysis, keeping the original as
  `obsm['composition_ncv']` for one comparison figure.
- **`palette(ECM_CELL_TYPES)`** appends `ecm` after the eight biology types in the glasbey
  cell-type slice, so it gets its own colour (`#5901a3`) and every existing biology colour
  is unchanged.

Note `main_2_featurize` rebuilds its datamodule from `model._hparams`, so the 11-channel
`categories_to_channels`, `pixel_size` and `global_size` all come from the checkpoint —
nothing about channel count is passed on the command line.

### Window membership must be tested in pixels

`imc_tm.window_mass` assigns rows to windows by integer `patch_xywh` in raster pixels,
which is the cropper's own test (`dataset.py:361`), not by a micron box around the window
centre. A micron box round-trips through a `floor()` and disagrees at the boundary by up
to one pixel. On densely packed cells that is a few percent and invisible; on ECM, which
sits on a 20 um lattice, one pixel of slop flips a whole row of voxels in or out, and a
micron-space alignment check reports a spurious 0.09-0.19 disagreement. Residual after the
pixel-exact test is pixel collisions only -- two rows in one 2 um pixel merge into one
raster entry -- so it scales as 1/occupancy: ~0.007 median on 100-200 entry windows.

Featurized output is ~5 GB per config. `.gitignore` matches `tissuemosaic/imc_anndata*/`
for that reason; the previous per-folder entries would not have caught the new screens.

## Spatial coordinates

The ICS CSVs are in **microns**, origin at the **ROI centre**, with **y increasing
upward**. QuPath exports `Centroid X/Y µm` already micron-calibrated from the OME
metadata, in the image frame (origin top-left, y down); `assemble_initial_conditions.py`
then applies

```python
cells['x'] = cells['x'] - img_w / 2
cells['y'] = img_h / 2 - cells['y']      # note the sign -- y is FLIPPED
```

Reconstructing `JHH369ROI1` from its detection file with `img_w = img_h = 1100`
reproduces the ICS coordinates to `max|diff| = 0.0` on both axes; dropping the flip
gives `64.32 µm`. Notebook 0 runs this proof.

**So plot `y` as-is.** Matplotlib's default (y up) reproduces the tissue in its
original anatomical orientation -- never `invert_yaxis()` or negate `y`. For rasterised
images use `origin='lower'`.

ROI domains are not uniform: 25 are 1000x1000 um, 22 are 1100x1100 um, one is
1000x580 um. Harmless -- `SparseImage` derives its own bounding box from the data plus
padding. A handful of edge-touching cells have centroids up to ~6 um outside the
declared domain. `volume` is the sphere volume from QuPath's `Max diameter µm`,
i.e. `(pi/6)*d^3`, in um^3.

## Pipeline order

**The inputs are committed, but they live in the `PhysiCell` submodule.** Run
`git submodule update --init --recursive` from the repo root first, or nothing below has
anything to read.

| # | stage | made by | reads | writes |
|---|-------|---------|-------|--------|
| 0 | inputs | `prep_imc_spatial/` (its own pipeline) | `qupath_detections/`, OME-TIFFs | `PhysiCell/.../config/ics/JHH_IMC/*.csv` (50 cell ICs) and `.../ics/substrates/*_ecm_scaled.csv` (51 ECM fields) -- 41 MB, tracked in the submodule |
| 1 | base AnnData | `make_imc_anndata.py` | those ICS CSVs + `prep_imc_spatial/imc_spatial_roi_specs.csv` | `tissuemosaic/imc_anndata/` -- 49 h5ad, 52 MB |
| 2 | per-config rasters | whichever driver needs them | `imc_anndata/` | `tissuemosaic/imc_anndata_<config>/` -- 9-126 MB each |
| 3 | training | `run_window_sweep.py`, `run_final_sweep.py`, `run_ecm_screen.py` | the per-config folder | `tissuemosaic/runs/<screen>/<config>/` -- 334 MB per checkpoint, up to 40 per run |
| 4 | featurisation | `run_featurize_screen.py` | h5ad + one checkpoint | `tissuemosaic/imc_anndata_feat_<config>/` -- ~5.2 GB per config, plus `_provenance.json` |
| 5 | evaluation | `notebooks/` | featurized h5ad, or checkpoints directly | figures in-notebook, `tissuemosaic/ecm_*.csv` |

```bash
conda activate tissue-mosaic-260819                    # conda_env_configs/tissuemosaic_260818.yaml

python tissuemosaic/make_imc_anndata.py                # stage 1

python tissuemosaic/run_ecm_screen.py --parallel       # stages 2+3, 200 epochs, ~5.6 h
python tissuemosaic/run_ecm_screen.py --long --parallel  # 1000 epochs, ~28 h

python tissuemosaic/run_featurize_screen.py --screen ecm --epoch 199   # stage 4
jupyter lab tissuemosaic/notebooks/                    # stage 5
```

Every driver takes `--dry-run`, which writes the configs and prints the plan without
training. Stage 2 is implicit: each driver builds its own h5ad folder if it is missing.

### Why one h5ad folder per config

`AnndataFolderDM` caches its rasterisation to `<data_folder>/train_dataset.pt`. The cache
is built for **one** `pixel_size` and channel count, so a stale one trains on the wrong
image silently, and two configs sharing a folder would race to write it. Hence
`imc_anndata_ecm_fov96/` and `imc_anndata_ecm_fov128/` rather than one shared directory,
and hence the drivers clearing the cache before every real launch (but not on `--dry-run`).

## Generated directories

**None of this is tracked** -- ~121 GB, regenerable from stage 0.

| path | what | size |
|------|------|------|
| `tissuemosaic/imc_anndata/` | stage 1 output, all 49 ROIs | 52 MB |
| `tissuemosaic/imc_anndata_<config>/` | the per-config subset + its raster cache | 9-126 MB each |
| `tissuemosaic/imc_anndata_feat_<config>/` | featurized, carries the full sparse-image state | 5.2 GB each, ~36 GB |
| `tissuemosaic/runs/<screen>/<config>/` | checkpoints under `.neptune/offline-name/OFFLINE/checkpoints/` | ~84 GB |

Training output sits under `tissuemosaic/`, not at the repo root, so this directory holds
the code and everything it generates -- and because `runs/<label>/` already means something
else in `slurm/` (the `pdac-spatial-pipeline` repo on the cluster, not this one).

## What the data looks like

50 IC CSVs, of which 48 are used. `JHH387ROI1_withductfiller.csv` and
`..._not_movable.csv` are skipped: both are variants of an ROI already in the set,
both hold the only NaN coordinate rows (1386 each), and they are the only source
of `epithelial_normal`, which no other ROI has.

- 48 ROIs, 14 patients, 357,274 cells
- 2.6k-13.7k cells per ROI, each ROI ~1000 x 1000 um (two are 1600 x 1600)
- median nearest-neighbour spacing 6.9 um
- 10 cell types, in the fixed channel order set by `CELL_TYPES` in
  `make_imc_anndata.py`, which **must** match `categories_to_channels` in the yaml

Composition is heavily skewed. `other_tissue` averages 57% of cells per ROI and
the hex-packed synthetic `duct_filler` another 14%, so 71% of the signal is
background and modelling artifact. The antigen-presentation biology sits in the
tail: `apCAF` 3.3%, `epithelial_tumor_class1_class2` 1.2%,
`mesenchymal_tumor_class1_class2` 0.3%. All 10 are kept as channels so that the
same model can also featurize PhysiCell simulation output, which contains
`duct_filler` and `other_tissue` too.

## Why the config differs from the testis defaults

`config_dino_ssl_imc.yaml` is a fork of `TissueMosaic/run/config_dino_ssl_testis.yaml`.
Every changed line is marked `# IMC`. The substantive one is spatial scale:

| pixel_size | global_size | field of view | non-overlapping crops in the cohort |
|---|---|---|---|
| 4.0 (testis default) | 96 | 384 um | 300 |
| 2.5 | 96 | 240 um | 760 |
| **2.0 (used here)** | **96** | **192 um** | **1190** |
| 1.7 | 64 | 109 um | 4270 |

The testis defaults were tuned for Slide-seq pucks and leave too little tissue per
ROI. 2.0 um is about a third of the 6.9 um cell spacing, which is the rule of
thumb in `SparseImage.from_anndata`, and a 192 um crop covers a tumour nest plus
its stromal cuff. Training draws *overlapping* random crops
(`n_crops_for_tissue_train` per ROI per epoch), so the table is a floor on
distinct content, not the number of samples seen.

`batch_size_per_gpu` is 128 rather than 256 to fit a 16 GB card.

Note `main_1_train_ssl.py` resolves parameters in the order **yaml > CLI >
defaults**, so anything present in the yaml cannot be overridden from the command
line. To vary a parameter, copy the yaml.

## tm_compat.py

The rebuilt environment (python 3.10, numpy 1.26.4, pytorch-lightning 1.9.4,
neptune-client 1.14.0) resolves three of the four API breakages natively. Only one
remains, and it cannot be fixed by version choice: torch >= 2.6 flipped
`torch.load`'s `weights_only` default to True, while the GPU floors torch at >= 2.7
(RTX 5060 Ti is sm_120 Blackwell; the oldest conda-forge pytorch with a CUDA 12.8+
build is 2.7.1).

`tm_compat` fixes that **without disabling the check globally**:

1. `add_safe_globals` for numpy's scalar/dtype types -- purely additive, so Lightning
   checkpoints load with `weights_only=True` still enforced.
2. `weights_only=False` only for `train_dataset.pt` / `test_dataset.pt`, the two files
   TissueMosaic pickles itself mid-run. Their allowlist tail runs through `slice`,
   pandas and anndata internals and does not converge, so relaxation is unavoidable
   there -- and is scoped to exactly those names.

Anything else you `torch.load` is still refused. Scoping is by basename, so any file
with those names is relaxed wherever it lives; that is the intended trade.

`tissuemosaic` is installed as a **non-editable snapshot** (`direct_url.json` has
`"dir_info": {}`, and the site-packages copy has a different inode from
`~/TissueMosaic/src/...`), so patching the upstream source would have no effect until
a reinstall and would silently revert on any env rebuild. The wrapper is version
controlled here instead. Worth reporting upstream regardless -- TissueMosaic is broken
on torch >= 2.6 for everyone.

## Not applicable

`TissueMosaic/notebooks/notebook3.ipynb` (gene regression) needs an expression
matrix. The IC CSVs carry only position, type and volume, so there is nothing to
regress. The raw marker intensities do exist -- `qupath_detections/*.txt` has 54
channels including HLA-DR, PD-L1, PD-1, FOXP3 and GZMB -- and `category_key`
accepts continuous weights, so a marker-intensity variant of this pipeline is
possible. It would also recover the labels `TYPE_MAP` collapses: `Treg`,
`macrophage`, `B cell`, and HLA-DR status on T cells.

## Generated directories

`imc_anndata/` (step 1 output, plus the `train_dataset.pt` / `test_dataset.pt`
caches step 2 writes beside it) is regenerable and not tracked.
