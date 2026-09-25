"""
make_imc_wellmixed_ics.py
-------------------------
Builds the PCMM initial conditions for IMC well-mixed: each ROI's cells from
the IMC spatial IC, placed at random, one layout per replicate.

IMC well-mixed runs the IMC spatial model (same config, rules, custom code and
per-ROI cell volumes); only the starting arrangement differs. For each ROI:

  - Cells: the rows of data/inputs/ics/cells/<ROI>/cells.csv minus the
    structural types (other_tissue, duct_filler), which HTAN well-mixed has no
    counterpart for. Type, is_movable and measured volume are kept; only the
    positions are new.
  - Disk: centred at the origin of the +/-800 um domain (HTAN well-mixed's),
    with area equal to the ROI's image area minus the area its structural cells
    cover -- the union of their discs on a 1 um raster, so overlaps are not
    counted twice. The biological cells then have the free space they had in
    the spatial ROI.
  - Placement: uniform random start in the disk, then overlapping pairs are
    pushed apart until no two cells overlap, every cell kept inside the disk.
    A cell's radius is PhysiCell's, from its volume.
  - ECM: uniform at the mean of the ROI's ECM image, on the domain's voxel
    grid. ECM never changes during a run, so this is its value throughout.

Each layout gets its own seed from (base seed, ROI, layout number), so
re-running this rebuilds the same files and raising n_layouts adds layouts
without touching existing ones. Each layout is its own PCMM IC folder:
PCMM counts runs as replicates only when they share input folders, and its own
per-run random placement (cells.xml patches) cannot carry per-cell volume or
is_movable, or keep cells apart.

Outputs
  data/inputs/ics/cells/<ROI>_wellmixed_r<k>/cells.csv     k = 1..n_layouts, x48
  data/inputs/ics/substrates/<ROI>_wellmixed/substrates.csv  x48
  prep_imc_spatial/imc_wellmixed_roi_specs.csv   one row per ROI: disk, areas, ECM
  prep_imc_spatial/imc_wellmixed_layouts.csv     one row per layout: seed, SHA-1

Run from the repo root, in the physicell-sim-260606 environment (it needs only
numpy and scipy, both pinned there), after setup_imc_spatial_pcmm.py, whose IMC
spatial ICs and spec table this reads:
    conda activate physicell-sim-260606
    python prep_imc_spatial/make_imc_wellmixed_ics.py
slurm/run_imc_wellmixed.jl then runs every layout; its N_LAYOUTS must not
exceed the n_layouts built here.
"""

import csv
import hashlib
import math
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

SETTINGS = {
    # Seeds: one per (base_seed, ROI, layout number).
    "base_seed": 20260925,
    # Layouts (replicates) per ROI. Raise to add replicates; existing layouts
    # are rebuilt unchanged.
    "n_layouts": 1,
    # Structural rows of the spatial IC, left out of the well-mixed runs.
    "excluded_types": ("other_tissue", "duct_filler"),
    # Raster step for the area the structural cells cover (um).
    "raster_um": 1.0,
    # Domain half-width (um), as HTAN well-mixed, and the voxel size of the
    # IMC config (um); the ECM file covers every voxel of this domain.
    "domain_half_width_um": 800.0,
    "voxel_um": 20.0,
    # Placement: how far apart a pushed pair is set beyond touching (um), and
    # when to give up.
    "push_margin_um": 1e-4,
    "max_iterations": 20000,
}

HERE = Path(__file__).parent
BASE = HERE.parent
ICS = BASE / "data" / "inputs" / "ics"
SPATIAL_SPECS = HERE / "imc_spatial_roi_specs.csv"
ROI_SPECS_OUT = HERE / "imc_wellmixed_roi_specs.csv"
LAYOUTS_OUT = HERE / "imc_wellmixed_layouts.csv"


def radius(volume):
    """PhysiCell's cell radius for a volume (um^3 -> um)."""
    return (3.0 * np.asarray(volume, float) / (4.0 * math.pi)) ** (1.0 / 3.0)


def read_rows(path):
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        return header, [row for row in reader if row]


def covered_area(xy, radii, bounds, step):
    """Area of the union of discs inside the rectangle `bounds`, on a raster."""
    x_min, x_max, y_min, y_max = bounds
    nx = int(round((x_max - x_min) / step))
    ny = int(round((y_max - y_min) / step))
    mask = np.zeros((ny, nx), bool)
    for (x, y), r in zip(xy, radii):
        i0 = max(int((x - r - x_min) / step), 0)
        i1 = min(int((x + r - x_min) / step) + 1, nx)
        j0 = max(int((y - r - y_min) / step), 0)
        j1 = min(int((y + r - y_min) / step) + 1, ny)
        if i0 >= i1 or j0 >= j1:
            continue
        px = x_min + (np.arange(i0, i1) + 0.5) * step
        py = y_min + (np.arange(j0, j1) + 0.5) * step
        mask[j0:j1, i0:i1] |= (px[None, :] - x) ** 2 + (py[:, None] - y) ** 2 <= r * r
    return mask.sum() * step * step


def overlapping_pairs(pos, radii):
    pairs = cKDTree(pos).query_pairs(2.0 * radii.max(), output_type="ndarray")
    if len(pairs) == 0:
        return pairs, np.empty(0), np.empty(0)
    dist = np.linalg.norm(pos[pairs[:, 1]] - pos[pairs[:, 0]], axis=1)
    overlap = radii[pairs[:, 0]] + radii[pairs[:, 1]] - dist
    keep = overlap > 0
    return pairs[keep], dist[keep], overlap[keep]


def place(radii, disk_radius, rng):
    """Non-overlapping positions for discs of `radii` inside a disk; returns (pos, iterations)."""
    limit = disk_radius - radii
    if (limit <= 0).any():
        raise ValueError("a cell is wider than the disk")
    u, theta = rng.random(len(radii)), rng.random(len(radii)) * 2.0 * math.pi
    pos = np.c_[np.sqrt(u) * limit * np.cos(theta), np.sqrt(u) * limit * np.sin(theta)]
    margin = SETTINGS["push_margin_um"]
    for iteration in range(SETTINGS["max_iterations"]):
        pairs, dist, overlap = overlapping_pairs(pos, radii)
        if len(pairs) == 0:
            return pos, iteration
        i, j = pairs[:, 0], pairs[:, 1]
        direction = (pos[j] - pos[i]) / np.maximum(dist, 1e-9)[:, None]
        shift = (overlap / 2.0 + margin)[:, None] * direction
        np.add.at(pos, i, -shift)
        np.add.at(pos, j, shift)
        r = np.linalg.norm(pos, axis=1)
        out = r > limit
        pos[out] *= (limit[out] / r[out])[:, None]
    raise RuntimeError(f"cells still overlap after {SETTINGS['max_iterations']} iterations")


def write_uniform_ecm(path, value):
    half, dx = SETTINGS["domain_half_width_um"], SETTINGS["voxel_um"]
    centres = np.arange(-half + dx / 2.0, half, dx)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write("x,y,z,ecm\n")
        for y in centres:
            for x in centres:
                fh.write(f"{x:g},{y:g},0,{value!r}\n")


def main():
    with open(SPATIAL_SPECS, newline="") as fh:
        specs = list(csv.DictReader(fh))
    if not specs:
        raise SystemExit(f"No ROIs in {SPATIAL_SPECS}")

    excluded = set(SETTINGS["excluded_types"])
    roi_rows, layout_rows = [], []
    for spec in specs:
        roi = spec["roi"]
        bounds = tuple(float(spec[k]) for k in ("x_min", "x_max", "y_min", "y_max"))
        image_area = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2])

        header, rows = read_rows(ICS / "cells" / roi / "cells.csv")
        ix, iy, it, iv = (header.index(k) for k in ("x", "y", "type", "volume"))
        structural = [r for r in rows if r[it] in excluded]
        cells = [r for r in rows if r[it] not in excluded]

        s_xy = np.array([[float(r[ix]), float(r[iy])] for r in structural]).reshape(-1, 2)
        s_area = covered_area(s_xy, radius([float(r[iv]) for r in structural]), bounds,
                              SETTINGS["raster_um"]) if structural else 0.0
        free_area = image_area - s_area
        disk_radius = math.sqrt(free_area / math.pi)
        radii = radius([float(r[iv]) for r in cells])
        coverage = float((math.pi * radii ** 2).sum() / free_area)

        _, ecm_rows = read_rows(ICS / "substrates" / roi / "substrates.csv")
        ecm_mean = float(np.mean([float(r[3]) for r in ecm_rows]))
        write_uniform_ecm(ICS / "substrates" / f"{roi}_wellmixed" / "substrates.csv", ecm_mean)

        roi_crc = zlib.crc32(roi.encode())
        for k in range(1, SETTINGS["n_layouts"] + 1):
            rng = np.random.default_rng(np.random.SeedSequence([SETTINGS["base_seed"], roi_crc, k]))
            pos, iterations = place(radii, disk_radius, rng)
            folder = f"{roi}_wellmixed_r{k}"
            out = ICS / "cells" / folder / "cells.csv"
            out.parent.mkdir(parents=True, exist_ok=True)
            with open(out, "w", newline="") as fh:
                w = csv.writer(fh, lineterminator="\n")
                w.writerow(header)
                for row, (x, y) in zip(cells, pos):
                    row = list(row)
                    row[ix], row[iy] = f"{x:.6f}", f"{y:.6f}"
                    w.writerow(row)
            layout_rows.append({
                "roi": roi, "layout": k, "folder": folder,
                "seed_base": SETTINGS["base_seed"], "seed_roi": roi_crc,
                "iterations": iterations,
                "sha1": hashlib.sha1(out.read_bytes()).hexdigest(),
            })

        roi_rows.append({
            "roi": roi, "n_cells": len(cells),
            "image_area_um2": round(image_area, 1),
            "structural_area_um2": round(s_area, 1),
            "free_area_um2": round(free_area, 1),
            "disk_radius_um": round(disk_radius, 3),
            "coverage": round(coverage, 4),
            "ecm_mean": round(ecm_mean, 6),
        })
        print(f"{roi}: {len(cells)} cells in r={disk_radius:.1f} um "
              f"(coverage {coverage:.2f}), ECM {ecm_mean:.2f}")

    for path, rows in ((ROI_SPECS_OUT, roi_rows), (LAYOUTS_OUT, layout_rows)):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    print(f"{len(roi_rows)} ROIs, {len(layout_rows)} layouts -> "
          f"{ROI_SPECS_OUT.relative_to(BASE)}, {LAYOUTS_OUT.relative_to(BASE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
