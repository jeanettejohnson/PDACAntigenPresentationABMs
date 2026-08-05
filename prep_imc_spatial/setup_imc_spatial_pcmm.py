"""
setup_imc_spatial_pcmm.py
-------------------------
Builds the PhysiCellModelManager project inputs for the IMC spatial sweep, so
the 48 per-ROI simulations can run through PCMM instead of the hand-rolled
sbatch loop.

The 48 PhysiCell_settings_<ROI>.xml files remain the source of truth: this
script reads them and derives everything else. Of 4435 leaf elements in those
configs only 33 differ across ROIs, and they split cleanly into

  - 2 initial-condition file references -> PCMM ic_cell / ic_substrate folders
  - 4 domain bounds                     -> config DiscreteVariations
  - 26 cell volumes (13 types x 2)      -> config DiscreteVariations
  - 1 save/folder                       -> irrelevant; PCMM owns output paths

Outputs
  data/inputs/configs/antigen_presentation/          PhysiCell_settings.xml (base)
  data/inputs/rulesets_collections/antigen_presentation/  base_rulesets.csv
  data/inputs/ics/cells/<ROI>/cells.csv              x48
  data/inputs/ics/substrates/<ROI>/substrates.csv    x48
  imc_spatial_roi_specs.csv                          one row per ROI, 30 varying values

Not generated (tracked in git, edit directly):
  data/inputs/custom_codes/antigen_presentation/     main.cpp, Makefile, custom_modules/

Re-runnable: everything it writes is derived and safe to regenerate. Note the
IC folders sit one step downstream of batch_assemble_ics.py -- if that is re-run,
re-run this too so the PCMM inputs match the assembler's current output.

Source-of-truth caveat: the cell-position <folder> in the ROI configs reads
"config/ics/JHH_IMC", which is the make-load copy -- untracked and already
stale. This script takes the *filename* from the config (authoritative for
which file a ROI uses) but resolves it against the git-tracked
user_projects/antigen_presentation/config/ics/JHH_IMC.
"""

import csv
import re
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# Inputs and the PCMM data/ tree resolve against the repo root; the spec table is
# written next to this script, with the rest of the prep pipeline's outputs.
HERE = Path(__file__).parent
BASE = HERE.parent
PROJ = "antigen_presentation"

USER_PROJECT = BASE / "PhysiCell" / "user_projects" / PROJ
CONFIG_DIR = USER_PROJECT / "config"
TRACKED_CELL_IC_DIR = CONFIG_DIR / "ics" / "JHH_IMC"
PHYSICELL_DIR = BASE / "PhysiCell"

INPUTS = BASE / "data" / "inputs"
SPEC_PATH = HERE / "imc_spatial_roi_specs.csv"

# The base config is one ROI's file; every value the 30 variations do not
# override comes from it. Cell types whose volumes are constant across all 48
# ROIs are supplied here rather than varied.
BASE_CONFIG_ROI = "JHH317ROI1"

# The 13 cell types whose volumes vary by ROI (measured per-ROI from the IMC
# data). The other 8 defined types are constant and come from the base config.
VARYING_VOLUME_TYPES = [
    "CAF",
    "CD4_Tcell",
    "CD8_Tcell",
    "CD8_exhausted",
    "Treg",
    "apCAF",
    "epithelial_tumor_class1",
    "epithelial_tumor_class1_class2",
    "epithelial_tumor_class2",
    "mesenchymal_tumor_class1",
    "mesenchymal_tumor_class1_class2",
    "mesenchymal_tumor_class2",
    "other_tissue",
]

DOMAIN_KEYS = ["x_min", "x_max", "y_min", "y_max"]

ROI_FROM_CONFIG = re.compile(r"^PhysiCell_settings_(.+)\.xml$")


def roi_configs():
    """(roi_key, path) for each per-ROI config, sorted by ROI key."""
    out = []
    for path in sorted(CONFIG_DIR.glob("PhysiCell_settings_JHH*ROI*.xml")):
        m = ROI_FROM_CONFIG.match(path.name)
        if m:
            out.append((m.group(1), path))
    return out


def text_at(root, path, attr_name=None):
    node = root.find(path)
    if node is None:
        raise KeyError(f"missing XML path: {path}")
    return (node.text or "").strip()


def read_spec(roi, path):
    """Pull the 30 varying values plus both IC references out of one ROI config."""
    root = ET.parse(path).getroot()

    spec = {"roi": roi}
    for key in DOMAIN_KEYS:
        spec[key] = text_at(root, f"domain/{key}")

    volumes = {c.get("name"): c.find("phenotype/volume") for c in root.find("cell_definitions")}
    for ct in VARYING_VOLUME_TYPES:
        vol = volumes.get(ct)
        if vol is None:
            raise KeyError(f"{roi}: cell type {ct} has no volume block")
        spec[f"vol_{ct}_total"] = (vol.find("total").text or "").strip()
        spec[f"vol_{ct}_nuclear"] = (vol.find("nuclear").text or "").strip()

    # IC references: take the filename from the config (authoritative), but
    # resolve cell positions against the tracked directory rather than the
    # stale make-load copy the <folder> element points at.
    cell_file = text_at(root, "initial_conditions/cell_positions/filename")
    spec["_cell_src"] = TRACKED_CELL_IC_DIR / cell_file

    sub_rel = text_at(root, "microenvironment_setup/options/initial_condition/filename")
    spec["_substrate_src"] = PHYSICELL_DIR / sub_rel

    return spec


def copy_into(src: Path, dest_dir: Path, dest_name: str):
    dest_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest_dir / dest_name)


def main():
    configs = roi_configs()
    if not configs:
        raise SystemExit(f"No ROI configs found in {CONFIG_DIR}")

    specs = []
    missing = []
    for roi, path in configs:
        spec = read_spec(roi, path)
        for key in ("_cell_src", "_substrate_src"):
            if not spec[key].is_file():
                missing.append(f"{roi}: {spec[key]}")
        specs.append(spec)

    if missing:
        print("Missing referenced IC files:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        raise SystemExit(1)

    # --- custom code -------------------------------------------------------
    # NOT generated. data/inputs/custom_codes/antigen_presentation/ is tracked in
    # git, matching how antigen_presentation_htan_singlecell already works.
    #
    # It used to be rmtree'd and re-copied from the user project on every run,
    # which meant any hand edit here -- the SVG palette in custom.cpp, say --
    # was silently destroyed the next time this script ran. Edit the tracked
    # copy directly, and keep it in sync with the user project
    # (tools/check_cell_colors.py checks the palette half of that).
    cc = INPUTS / "custom_codes" / PROJ
    missing = [p for p in ("main.cpp", "Makefile", "custom_modules/custom.cpp")
               if not (cc / p).is_file()]
    if missing:
        raise SystemExit(
            f"Tracked custom code is incomplete under {cc.relative_to(BASE)}: "
            f"missing {', '.join(missing)}. It is version-controlled, not generated -- "
            "restore it from git rather than re-running this script."
        )
    print(f"custom_code   -> {cc.relative_to(BASE)}  (tracked, not generated)")

    # --- base config -------------------------------------------------------
    base_src = CONFIG_DIR / f"PhysiCell_settings_{BASE_CONFIG_ROI}.xml"
    cfg = INPUTS / "configs" / PROJ
    cfg.mkdir(parents=True, exist_ok=True)

    # Repoint cell positions at the tracked IC directory. PCMM passes -i so this
    # element is inert at runtime, but a correct path means a manual run of the
    # base config resolves to tracked data instead of the stale make-load copy.
    tree = ET.parse(base_src)
    folder = tree.getroot().find("initial_conditions/cell_positions/folder")
    if folder is not None:
        folder.text = "user_projects/antigen_presentation/config/ics/JHH_IMC"

    # Run single-threaded, matching the other three simulations and the
    # cpus-per-task=1 that slurm_common.jl requests. The per-ROI configs ask for
    # 8 threads, which only made sense under the old runner's own allocation;
    # one thread per job trades wall time for a 8x smaller core footprint and
    # lets 48 ROIs run in one scheduling batch.
    omp = tree.getroot().find("parallel/omp_num_threads")
    if omp is None:
        raise KeyError("missing XML path: parallel/omp_num_threads")
    omp.text = "1"

    tree.write(cfg / "PhysiCell_settings.xml", encoding="UTF-8", xml_declaration=True)
    print(f"config        -> {(cfg / 'PhysiCell_settings.xml').relative_to(BASE)}  (from {BASE_CONFIG_ROI})")

    # --- rulesets ----------------------------------------------------------
    rules = INPUTS / "rulesets_collections" / PROJ
    rules.mkdir(parents=True, exist_ok=True)
    shutil.copy2(CONFIG_DIR / "cell_rules.csv", rules / "base_rulesets.csv")
    print(f"rulesets      -> {(rules / 'base_rulesets.csv').relative_to(BASE)}")

    # --- per-ROI initial conditions ---------------------------------------
    for spec in specs:
        roi = spec["roi"]
        copy_into(spec["_cell_src"], INPUTS / "ics" / "cells" / roi, "cells.csv")
        copy_into(spec["_substrate_src"], INPUTS / "ics" / "substrates" / roi, "substrates.csv")
    print(f"ic_cell       -> {len(specs)} folders under data/inputs/ics/cells/")
    print(f"ic_substrate  -> {len(specs)} folders under data/inputs/ics/substrates/")

    # --- spec table --------------------------------------------------------
    fields = ["roi"] + DOMAIN_KEYS
    for ct in VARYING_VOLUME_TYPES:
        fields += [f"vol_{ct}_total", f"vol_{ct}_nuclear"]

    with open(SPEC_PATH, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for spec in specs:
            writer.writerow({k: spec[k] for k in fields})
    print(f"spec table    -> {SPEC_PATH.relative_to(BASE)}  ({len(specs)} ROIs x {len(fields) - 1} values)")


if __name__ == "__main__":
    main()
