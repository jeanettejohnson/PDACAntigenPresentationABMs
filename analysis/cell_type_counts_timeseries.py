"""
Cell-type counts at every timepoint for all HT* simulation runs.

For each run, iterates over every output XML file and records the number
of cells of each type at that timepoint.

Writes (in analysis/cell_type_counts/):
  - <run>_cell_type_counts.csv  : tidy long table for that run
  - all_cell_type_counts.csv    : all runs combined

Output columns: run, time, cell_type, count
"""

import os
import glob
import time
import traceback
import pandas as pd
import pcdl

OUTPUTS_DIR  = os.path.join(os.path.dirname(os.path.dirname(__file__)), "PhysiCell", "outputs")
OUT_DIR      = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures", "cell_type_counts")


def get_ht_folders():
    return sorted(glob.glob(os.path.join(OUTPUTS_DIR, "HT*")))


def get_xml_files(run_path: str) -> list[str]:
    """Return sorted list of output_*.xml files (excludes initial.xml / final.xml)."""
    xmls = sorted(glob.glob(os.path.join(run_path, "output*.xml")))
    return xmls


def counts_for_run(run_path: str) -> pd.DataFrame:
    run_name = os.path.basename(run_path)
    xmls = get_xml_files(run_path)
    if not xmls:
        print(f"  SKIP — no output*.xml found in {run_name}")
        return pd.DataFrame()

    rows = []
    n = len(xmls)
    t0 = time.time()

    for i, xml in enumerate(xmls):
        mcds = pcdl.TimeStep(xml, microenv=False, graph=False, verbose=False)
        t    = mcds.get_time()
        df   = mcds.get_cell_df()[["cell_type"]]
        vc   = df["cell_type"].value_counts()
        for cell_type, count in vc.items():
            rows.append({
                "run":       run_name,
                "time":      t,
                "cell_type": cell_type,
                "count":     int(count),
            })

        if (i + 1) % 50 == 0 or (i + 1) == n:
            elapsed = time.time() - t0
            rate    = (i + 1) / elapsed
            eta     = (n - i - 1) / rate if rate > 0 else 0
            print(f"    {i+1:3d}/{n}  t={t:.0f} min  "
                  f"[{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining]")

    return pd.DataFrame(rows)


def main():
    total_start = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    folders = get_ht_folders()
    total   = len(folders)
    print(f"Found {total} HT* simulation folders.")
    print(f"Output directory: {OUT_DIR}\n")

    all_frames = []
    succeeded  = []
    failed     = []

    for idx, folder in enumerate(folders, start=1):
        run_name = os.path.basename(folder)
        print(f"[{idx}/{total}] {run_name}")

        try:
            df = counts_for_run(folder)
        except Exception as e:
            print(f"  ERROR: {e}")
            traceback.print_exc()
            failed.append(run_name)
            continue

        if df.empty:
            failed.append(run_name)
            continue

        out_path = os.path.join(OUT_DIR, f"{run_name}_cell_type_counts.csv")
        df.to_csv(out_path, index=False)
        print(f"  Saved {len(df):,} rows → {out_path}")

        all_frames.append(df)
        succeeded.append(run_name)

        elapsed = time.time() - total_start
        avg     = elapsed / idx
        eta     = avg * (total - idx)
        print(f"  Overall: {idx}/{total} runs | {elapsed:.0f}s elapsed | "
              f"~{eta:.0f}s (~{eta/60:.1f} min) remaining\n")

    # ── combined file ──────────────────────────────────────────────────────
    if all_frames:
        combined_path = os.path.join(OUT_DIR, "all_cell_type_counts.csv")
        pd.concat(all_frames, ignore_index=True).to_csv(combined_path, index=False)
        print(f"Saved combined: {combined_path}")

    elapsed_total = time.time() - total_start
    print(f"\n── Summary ──────────────────────────────────")
    print(f"  Succeeded: {len(succeeded)}/{total}")
    if failed:
        print(f"  Failed:    {failed}")
    print(f"  Total time: {elapsed_total:.1f}s ({elapsed_total/60:.1f} min)")
    print(f"  Outputs in: {OUT_DIR}/")


if __name__ == "__main__":
    main()
