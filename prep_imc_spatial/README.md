# IMC spatial initial-condition preparation

Everything that turns raw IMC imaging into the PCMM inputs `slurm/run_imc_spatial.jl`
consumes. Run from the **repo root**, e.g. `python3 prep_imc_spatial/batch_assemble_ics.py`
-- the scripts anchor paths with `Path(__file__).parent.parent`.

The funnel narrows at each stage: 56 detections -> 51 substrates -> 50 cell ICs ->
48 ROI configs -> 48 PCMM input pairs.

## Pipeline order

| # | script | in | out |
|---|--------|----|-----|
| 1 | `scale_biwt.py` | `ics/substrates/*_ecm.csv` | `*_ecm_scaled.csv` (95th-pct clip) |
| 2 | `batch_assemble_ics.py` -> `assemble_initial_conditions.py` | `qupath_detections/*.txt`, `ics/ductcoordinates/*.geojson`, `*_ecm_scaled.csv` | `ics/JHH_IMC/<ROI>.csv` |
| 3 | `generate_roi_configs.py` | cell ICs + substrates | `PhysiCell_settings_<ROI>.xml` x48 |
| 4 | `make_imc_count_summary.py` | cell ICs | `prep_imc_spatial/assignmentsummary_JHH_IMC.csv` |
| 5 | `setup_imc_spatial_pcmm.py` | the 48 ROI configs | `data/inputs/` + `prep_imc_spatial/imc_spatial_roi_specs.csv` |
| 6 | `make_imc_wellmixed_ics.py` | stage 5's ICs, substrates and spec table | `data/inputs/ics/{cells,substrates}/<ROI>_wellmixed*` + `imc_wellmixed_roi_specs.csv`, `imc_wellmixed_layouts.csv` |

Stage 2 is where antigen presentation is encoded: `TYPE_MAP` in
`assemble_initial_conditions.py` maps QuPath classes onto PhysiCell cell types
(`tumor_epithelial: HLA-DR` -> `epithelial_tumor_class1_class2`,
`CD4 T cell: FOXP3` -> `Treg`, ...), and duct polygons are hex-packed with
`duct_filler` cells.

Stage 2 writes **two** copies of each cell IC: the repo-root `PhysiCell/config/ics/JHH_IMC`
(the `make load` copy) and `PhysiCell/user_projects/antigen_presentation/config/ics/JHH_IMC`
(tracked, canonical). Stage 5 reads the canonical one.

### QC / inspection

`plot_initial_counts_imc.py`, `plot_apcaf_vs_tcells.py`, `plot_cell_size_by_type.py`
-- these write into `analysis/`, which is not created automatically.

## Ordering constraint

If the ICs are rebuilt, **stage 5 must be re-run**: the PCMM `ic_cell`/`ic_substrate`
folders are copies one step downstream, and `run_imc_spatial.jl` cannot detect that
they have gone stale.

## IMC well-mixed (stage 6)

`slurm/run_imc_wellmixed.jl` runs the IMC spatial model (same config, rules,
custom code and per-ROI volumes) with each ROI's cells placed at random.
`make_imc_wellmixed_ics.py` builds its inputs; run it in the
`physicell-sim-260606` environment, which has the numpy and scipy it needs.

- **Cells:** each ROI's IC rows minus the structural `other_tissue` and
  `duct_filler`, keeping type, `is_movable` and measured volume.
- **Disk:** placed in a disk centred in the +/-800 um domain (HTAN well-mixed's).
  Its area is the ROI's image area minus the area the structural cells cover (a
  1 um raster of their discs), so the biological cells keep the free space they
  had in the spatial ROI. No two cells overlap.
- **ECM:** uniform at the mean of the ROI's ECM image.
- **Layouts:** one random layout per replicate, `<ROI>_wellmixed_r<k>`, each
  seeded from (base seed, ROI, k). Re-running rebuilds identical files, and
  raising `n_layouts` adds layouts without changing existing ones; then raise
  `N_LAYOUTS` in the runner. Each layout is its own PCMM IC folder, so PCMM
  lists it as its own monad rather than a replicate of the others.

If stage 5 is re-run (so the spatial ICs change), re-run stage 6 too.

## Edge rows are clipped into the domain

Each ROI's domain is the extent of its ECM image, but QuPath places the centroids
of cells cut by the image border up to ~6 µm outside that frame (1,809 rows across
the 48 ROIs, 312 of them biological cells). PhysiCell freezes any cell it loads
outside the domain: it never moves, divides or dies, yet is written out as live.
Stage 5 therefore moves each such row `CLIP_INSET_UM` (0.5 µm) inside the edge as it
copies the IC, and stops if any row would move more than `CLIP_MAX_MOVE_UM`
(20 µm). Every other row is copied byte for byte, and the tracked ICs under
`user_projects/antigen_presentation/config/ics/JHH_IMC` stay as assembled.

## Known issue: stale root IC directory

`generate_roi_configs.py`, `plot_initial_counts_imc.py` and `plot_apcaf_vs_tcells.py`
read the repo-root `PhysiCell/config/ics/JHH_IMC`, which is **stale** -- an older
generation with a different schema, not merely an unprocessed copy:

|          | rows | columns                                    | volume     |
|----------|------|--------------------------------------------|------------|
| root     | 5202 | 9 (incl. `cycle entry`, `custom:GFP`, ...) | all NaN    |
| canonical| 6992 | 6                                          | populated  |

Of the 11 ROIs present in both, zero are identical, and the root copy holds only 22
CSVs against the canonical 50. Repoint these to the canonical directory before
relying on them; `setup_imc_spatial_pcmm.py` already does (it takes the *filename*
from each ROI config but resolves it against the canonical path).

## archive/

Superseded or one-off scripts from this pipeline, kept for provenance. None are
part of the current flow.

| script | why archived |
|--------|--------------|
| `scaleimage.py` | TIFF -> `_ecm.csv`; hardcoded `/Users/...` path, and writes to `ics/substrate` (singular), which nothing reads |
| `writequpathcoords2physicell.py` | superseded by `assemble_initial_conditions.py`; hardcoded OneDrive paths |
| `hex_mask_workflow.py` | hex-grid prototype; the packing now lives inside `assemble_initial_conditions.py` |
| `hex_mask_workflow (1).py` | byte-identical duplicate of the above (download artifact) |
| `winsorize_ics_volumes.py` | self-described one-off; targets the stale root IC directory |
| `fix_is_movable.py` | targets the stale root IC directory; the assembler now sets `is_movable` directly |
| `IMCSegmentationSteinbock.py` | upstream segmentation, run once off-cluster; hardcoded paths |
| `imctrying.py` | exploratory scratch; hardcoded paths |

Stage 1's true input (`_ecm.csv` from the ECM channel) is **not currently
reproducible in this repo** -- `scaleimage.py` is the only TIFF -> CSV script and it
is archived for the reasons above. Only its committed outputs exist.
