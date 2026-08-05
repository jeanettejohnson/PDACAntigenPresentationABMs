"""
generate_roi_configs.py
-----------------------
For each ROI that has both a cell ICS CSV and a scaled substrate CSV, produces
a copy of PhysiCell_settings.xml with:
  - domain bounds derived from the substrate voxel grid
  - output folder pointing to outputs/<ROI>
  - cell IC filename updated to the ROI's ICS
  - substrate IC filename updated to the ROI's ecm_scaled CSV
  - per-type target volumes updated from the ICS medians for that ROI
    (only for types with >= MIN_CELLS cells; others keep the global value)
"""

import re
import pandas as pd
import numpy as np
from pathlib import Path

BASE        = Path(__file__).parent
CONFIG_DIR  = BASE / "PhysiCell/user_projects/antigen_presentation/config"
DEPLOY_DIR  = BASE / "PhysiCell/config"   # where PhysiCell actually runs from
ICS_DIR     = BASE / "PhysiCell/config/ics/JHH_IMC"
SUB_DIR     = CONFIG_DIR / "ics/substrates"
TEMPLATE    = CONFIG_DIR / "PhysiCell_settings.xml"

NUCLEAR_FRACTION = 540.0 / 2494.0
DX               = 20.0
MIN_CELLS        = 5   # minimum cells of a type to trust its per-ROI median

PARENT_VOLUME_TYPE = {
    'Treg':                              'CD4_Tcell',
    'CD8_exhausted':                     'CD8_Tcell',
    'apCAF':                             'CAF',
    'epithelial_tumor_class1_class2':    'epithelial_tumor_class1',
    'epithelial_tumor_class2':           'epithelial_tumor_class1',
    'mesenchymal_tumor_class1_class2':   'mesenchymal_tumor_class1',
    'mesenchymal_tumor_class2':          'mesenchymal_tumor_class1',
}

def clean_roi_stem(stem):
    """Strip non-ROI suffixes like _ductfiller to get the canonical ROI name."""
    return re.sub(r'_(ductfiller.*)$', '', stem)

def compute_domain(sub_df):
    """Return (x_min, x_max, y_min, y_max) from substrate voxel centres."""
    half = DX / 2
    return (sub_df.x.min() - half, sub_df.x.max() + half,
            sub_df.y.min() - half, sub_df.y.max() + half)

def compute_type_volumes(ics_df):
    """Median volume per type, honouring parent-type inheritance."""
    if 'volume' not in ics_df.columns:
        return {}
    counts  = ics_df.groupby('type')['volume'].count()
    medians = ics_df.groupby('type')['volume'].median()
    result  = {}
    for t, med in medians.items():
        if counts[t] >= MIN_CELLS:
            result[t] = round(med)
    # induced types inherit parent's median where parent has enough cells
    for child, parent in PARENT_VOLUME_TYPE.items():
        if parent in result:
            result[child] = result[parent]
    return result

def update_volume_block(chunk, total, nuclear):
    vol_match = re.search(r'(<volume>.*?</volume>)', chunk, re.DOTALL)
    if not vol_match:
        return chunk
    vb = vol_match.group(1)
    vb = re.sub(r'(<total[^>]*>)[^<]*(</total>)',   rf'\g<1>{total}\g<2>',   vb, count=1)
    vb = re.sub(r'(<nuclear[^>]*>)[^<]*(</nuclear>)', rf'\g<1>{nuclear}\g<2>', vb, count=1)
    return chunk[:vol_match.start()] + vb + chunk[vol_match.end():]

def apply_volumes(text, type_volumes):
    if not type_volumes:
        return text
    parts = re.split(r'(?=<cell_definition name=")', text)
    out = []
    for part in parts:
        m = re.match(r'<cell_definition name="([^"]+)"', part)
        if m and m.group(1) in type_volumes:
            total   = type_volumes[m.group(1)]
            nuclear = round(total * NUCLEAR_FRACTION)
            part    = update_volume_block(part, total, nuclear)
        out.append(part)
    return ''.join(out)

# ── Find ROIs with both ICS and substrate ──────────────────────────────────────
SKIP_SUFFIXES = {'rectangle', 'withTreg', 'withductfiller', 'ductfiller', 'withTregs', 'not_movable'}

ics_files = {
    p.stem: p for p in ICS_DIR.glob("*.csv")
    if not any(s in p.stem for s in SKIP_SUFFIXES)
}

template_text = TEMPLATE.read_text()
generated = []

for ics_stem, ics_path in sorted(ics_files.items()):
    roi_stem = clean_roi_stem(ics_stem)          # e.g. JHH317ROI3 (strips _ductfiller)
    sub_path = SUB_DIR / f"{roi_stem}_ecm_scaled.csv"

    if not sub_path.exists():
        print(f"  SKIP {roi_stem}: no substrate file ({sub_path.name})")
        continue

    ics_df = pd.read_csv(ics_path)
    sub_df = pd.read_csv(sub_path)

    x_min, x_max, y_min, y_max = compute_domain(sub_df)
    type_volumes = compute_type_volumes(ics_df)

    text = template_text

    # domain
    text = re.sub(r'<x_min>[^<]*</x_min>', f'<x_min>{x_min}</x_min>', text)
    text = re.sub(r'<x_max>[^<]*</x_max>', f'<x_max>{x_max}</x_max>', text)
    text = re.sub(r'<y_min>[^<]*</y_min>', f'<y_min>{y_min}</y_min>', text)
    text = re.sub(r'<y_max>[^<]*</y_max>', f'<y_max>{y_max}</y_max>', text)

    # output folder
    text = re.sub(r'<folder>outputs/[^<]*</folder>',
                  f'<folder>outputs/{roi_stem}</folder>', text, count=1)

    # All paths below are written RELATIVE to PhysiCell/ (the executable's
    # cwd at runtime), not absolute. PhysiCell resolves <folder>/<filename>
    # via plain fstream against cwd, so this makes configs portable across
    # machines (confirmed broken twice already: old laptop -> new Mac ->
    # HPC cluster all have different absolute home paths).
    sub_path_rel = sub_path.relative_to(BASE / "PhysiCell")
    ics_dir_rel = ICS_DIR.relative_to(BASE / "PhysiCell")
    config_dir_rel = CONFIG_DIR.relative_to(BASE / "PhysiCell")

    # substrate IC filename — force enabled="true" regardless of template value
    text = re.sub(
        r'<initial_condition type="csv" enabled="[^"]*">(\s*<filename>)[^<]*(</filename>)',
        rf'<initial_condition type="csv" enabled="true">\g<1>{sub_path_rel}\g<2>',
        text, count=1
    )

    # cell IC filename and folder
    text = re.sub(
        r'(<cell_positions type="csv" enabled="true">\s*<folder>)[^<]*(</folder>\s*<filename>)[^<]*(</filename>)',
        rf'\g<1>{ics_dir_rel}\g<2>{ics_stem}.csv\g<3>',
        text, flags=re.DOTALL, count=1
    )

    # ruleset folder
    text = re.sub(
        r'(<ruleset[^>]*>\s*<folder>)[^<]*(</folder>)',
        rf'\g<1>{config_dir_rel}\g<2>',
        text, count=1
    )

    # per-ROI target volumes
    text = apply_volumes(text, type_volumes)

    out_path = CONFIG_DIR / f"PhysiCell_settings_{roi_stem}.xml"
    out_path.write_text(text)
    (DEPLOY_DIR / out_path.name).write_text(text)
    generated.append(roi_stem)
    n_types = len(type_volumes)
    print(f"  {roi_stem}: domain x[{x_min},{x_max}] y[{y_min},{y_max}]  "
          f"{n_types} type volumes updated  → {out_path.name}")

print(f"\nGenerated {len(generated)} config files.")
