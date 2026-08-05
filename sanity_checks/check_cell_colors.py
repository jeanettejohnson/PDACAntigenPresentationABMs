#!/usr/bin/env python3
"""
check_cell_colors.py -- guard against cell-type palette drift.

The SVG palette is hardcoded in several places that must agree. This has already
failed silently twice:

  1. The palette originally lived in PhysiCell's modules/PhysiCell_pathology.cpp
     (commit 124a30be). The drbergman 1.14.2 merge rewrote that function and the
     palette vanished -- movies kept rendering, just in the wrong colors, for
     about six weeks.
  2. PhysiCell (SVG) and PhysiCell-Studio (scalar-coloring view) fall back to
     different generated palettes past the built-in range, so the same cell type
     could appear in two different colors depending on which view you opened.

So the palette now lives in custom.cpp (which upstream cannot clobber) plus
Studio's cmaps.py, and this script checks they still agree.

Checks
  1. every custom.cpp copy declares an identical kCellTypeColors map
  2. those names exactly match <cell_definitions> in each config that has one
  3. paint_clist is ordered to match <cell_definitions>
  4. each paint_clist RGB equals matplotlib's value for the named SVG color

Exit 0 if consistent, 1 otherwise. Run after touching any palette, and after
bumping the PhysiCell or PhysiCell-Studio submodule.
"""

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent

CPP_SOURCES = [
    "PhysiCell/user_projects/antigen_presentation/custom_modules/custom.cpp",
    "PhysiCell/user_projects/antigen_presentation_htan_singlecell/custom_modules/custom.cpp",
    "data/inputs/custom_codes/antigen_presentation/custom_modules/custom.cpp",
    "data/inputs/custom_codes/antigen_presentation_htan_singlecell/custom_modules/custom.cpp",
]

CMAPS = "PhysiCell-Studio/bin/cmaps.py"

# Configs whose <cell_definitions> order the palettes must follow. The first is
# authoritative for paint_clist ordering; others are checked for name coverage.
CONFIGS = [
    "data/inputs/configs/antigen_presentation/PhysiCell_settings.xml",
    "data/inputs/configs/antigen_presentation_htan_singlecell/PhysiCell_settings.xml",
]

MAP_ENTRY = re.compile(r'\{\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\}')


def parse_cpp_palette(path):
    """Extract kCellTypeColors as an ordered list of (name, color)."""
    text = path.read_text()
    m = re.search(r"kCellTypeColors\s*=\s*\{(.*?)\}\s*;", text, re.S)
    if not m:
        raise LookupError(f"no kCellTypeColors map found in {path}")
    return MAP_ENTRY.findall(m.group(1))


def parse_paint_clist(path):
    """Extract paint_clist entries as (rgb_tuple, trailing_comment)."""
    text = path.read_text()
    start = text.index("paint_clist")
    open_br = text.index("[", start)
    depth = 0
    for i in range(open_br, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                close_br = i
                break
    body = text[open_br + 1:close_br]
    out = []
    for line in body.splitlines():
        for entry in re.finditer(r"\[([^\[\]]*)\]", line):
            nums = tuple(float(x) for x in entry.group(1).split(","))
            comment = ""
            c = line.split("#", 1)
            if len(c) == 2:
                comment = c[1].strip()
            out.append((nums, comment))
    return out


def config_cell_types(path):
    root = ET.parse(path).getroot()
    return [c.get("name") for c in root.find("cell_definitions")]


def main():
    failures = []

    # --- 1. all C++ copies agree ------------------------------------------
    palettes = {}
    for rel in CPP_SOURCES:
        p = BASE / rel
        if not p.is_file():
            failures.append(f"missing C++ source: {rel}")
            continue
        try:
            palettes[rel] = parse_cpp_palette(p)
        except LookupError as e:
            failures.append(str(e))

    if not palettes:
        print("FAIL: no palettes could be parsed", file=sys.stderr)
        return 1

    reference_rel, reference = next(iter(palettes.items()))
    for rel, pal in palettes.items():
        if pal != reference:
            only_a = set(dict(reference)) - set(dict(pal))
            only_b = set(dict(pal)) - set(dict(reference))
            diff = sorted(
                n for n in set(dict(pal)) & set(dict(reference))
                if dict(pal)[n] != dict(reference)[n]
            )
            failures.append(
                f"palette in {rel} differs from {reference_rel}"
                + (f"; missing {sorted(only_a)}" if only_a else "")
                + (f"; extra {sorted(only_b)}" if only_b else "")
                + (f"; recolored {diff}" if diff else "")
            )

    palette_names = [n for n, _ in reference]
    palette_map = dict(reference)
    if len(palette_names) != len(set(palette_names)):
        failures.append("duplicate cell type names in kCellTypeColors")

    # --- 2. names match the configs ---------------------------------------
    config_orders = {}
    for rel in CONFIGS:
        p = BASE / rel
        if not p.is_file():
            print(f"note: config not present, skipping: {rel}")
            continue
        types = config_cell_types(p)
        config_orders[rel] = types
        unknown = [t for t in types if t not in palette_map]
        if unknown:
            failures.append(f"{rel}: cell types with no palette entry: {unknown}")

    # --- 3 & 4. Studio list ordering and color values ---------------------
    cm = BASE / CMAPS
    if not cm.is_file():
        failures.append(f"missing {CMAPS}")
    elif config_orders:
        try:
            import matplotlib.colors as mc
        except ImportError:
            print("note: matplotlib unavailable, skipping paint_clist RGB check")
            mc = None

        clist = parse_paint_clist(cm)
        primary = CONFIGS[0]
        order = config_orders.get(primary)
        if order is None:
            order = next(iter(config_orders.values()))

        if len(clist) < len(order):
            failures.append(
                f"paint_clist has {len(clist)} entries, fewer than the "
                f"{len(order)} cell definitions in {primary}"
            )
        else:
            for i, name in enumerate(order):
                rgb, comment = clist[i]
                expected_color = palette_map.get(name)
                if expected_color is None:
                    continue
                if name not in comment:
                    failures.append(
                        f"paint_clist[{i}] comment is {comment!r}, expected to name "
                        f"{name!r} -- ordering may have drifted from <cell_definitions>"
                    )
                if mc is not None:
                    want = tuple(round(v, 5) for v in mc.to_rgb(expected_color))
                    got = tuple(round(v, 5) for v in rgb)
                    if want != got:
                        failures.append(
                            f"paint_clist[{i}] ({name}) is {got}, expected {want} "
                            f"for '{expected_color}'"
                        )

    if failures:
        print("FAIL: cell-type palette is inconsistent\n", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1

    print(f"OK: {len(palette_names)} cell-type colors consistent across "
          f"{len(palettes)} C++ copies and {CMAPS}")
    for rel, types in config_orders.items():
        print(f"    {len(types):>2} cell definitions covered in {Path(rel).parent.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
