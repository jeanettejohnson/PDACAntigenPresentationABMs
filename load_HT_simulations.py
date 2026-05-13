"""
Load and analyze all HT* PhysiCell simulation output folders.

Produces per-run outputs in analysis/:
  - <run>_cell_counts.csv   : cell counts by type at every timestep
  - <run>_final_cells.csv   : full cell DataFrame at the final timestep
  - <run>_substrates.csv    : mean substrate concentrations at every timestep

And two aggregate files:
  - all_cell_counts.csv     : all runs combined (includes 'run' column)
  - all_substrates.csv      : all runs combined
"""

import os
import glob
import time
import pandas as pd
import pcdl

OUTPUTS_DIR = os.path.join(os.path.dirname(__file__), "PhysiCell", "outputs")
ANALYSIS_DIR = os.path.join(os.path.dirname(__file__), "analysis")


def get_ht_folders():
    pattern = os.path.join(OUTPUTS_DIR, "HT*")
    return sorted(glob.glob(pattern))


def process_run(run_path: str, run_idx: int, total_runs: int) -> dict[str, pd.DataFrame]:
    run_name = os.path.basename(run_path)
    run_start = time.time()
    print(f"\n[{run_idx}/{total_runs}] ── {run_name} ──────────────────────────")
    print(f"  Loading TimeSeries from: {run_path}")

    mcdsts = pcdl.TimeSeries(
        output_path=run_path,
        microenv=True,
        graph=False,
        verbose=False,
    )

    mcds_list = mcdsts.get_mcds_list()
    n_steps = len(mcds_list)
    print(f"  Loaded {n_steps} timesteps")
    print(f"  Time range: {mcds_list[0].get_time():.1f} – {mcds_list[-1].get_time():.1f} min")
    print(f"  Substrates: {mcds_list[0].get_substrate_list()}")
    print(f"  Cell types present: {sorted(mcds_list[-1].get_cell_df()['cell_type'].unique())}")

    # ── cell counts over time ──────────────────────────────────────────────
    print(f"  Extracting cell counts over {n_steps} timesteps...")
    count_rows = []
    for i, mcds in enumerate(mcds_list):
        t = mcds.get_time()
        df = mcds.get_cell_df()
        type_counts = df["cell_type"].value_counts().to_dict()
        type_counts["time"] = t
        type_counts["run"] = run_name
        count_rows.append(type_counts)
        if (i + 1) % 50 == 0 or (i + 1) == n_steps:
            print(f"    ... {i+1}/{n_steps} timesteps processed (t={t:.1f} min, {len(df)} cells)")

    df_counts = pd.DataFrame(count_rows).fillna(0)
    cols = ["run", "time"] + [c for c in df_counts.columns if c not in ("run", "time")]
    df_counts = df_counts[cols]
    print(f"  Cell counts table: {df_counts.shape[0]} rows x {df_counts.shape[1]} cols")

    # ── substrate means over time ──────────────────────────────────────────
    print(f"  Extracting substrate concentrations over {n_steps} timesteps...")
    sub_rows = []
    for i, mcds in enumerate(mcds_list):
        t = mcds.get_time()
        df_conc = mcds.get_conc_df()
        sub_names = mcds.get_substrate_list()
        row = {"run": run_name, "time": t}
        for sub in sub_names:
            if sub in df_conc.columns:
                row[f"{sub}_mean"] = df_conc[sub].mean()
        sub_rows.append(row)
        if (i + 1) % 50 == 0 or (i + 1) == n_steps:
            print(f"    ... {i+1}/{n_steps} timesteps processed (t={t:.1f} min)")

    df_subs = pd.DataFrame(sub_rows)
    print(f"  Substrates table: {df_subs.shape[0]} rows x {df_subs.shape[1]} cols")

    # ── final timepoint full cell table ───────────────────────────────────
    print(f"  Extracting full cell table at final timestep (t={mcds_list[-1].get_time():.1f} min)...")
    final_mcds = mcds_list[-1]
    df_final = final_mcds.get_cell_df()
    df_final.insert(0, "run", run_name)
    print(f"  Final cell table: {df_final.shape[0]} cells x {df_final.shape[1]} columns")

    elapsed = time.time() - run_start
    print(f"  Completed {run_name} in {elapsed:.1f}s")

    return {"counts": df_counts, "substrates": df_subs, "final_cells": df_final}


def main():
    total_start = time.time()
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    folders = get_ht_folders()
    total = len(folders)
    print(f"Found {total} HT* simulation folders.")
    print(f"Output directory: {ANALYSIS_DIR}\n")

    all_counts = []
    all_subs = []
    succeeded = []
    failed = []

    for idx, folder in enumerate(folders, start=1):
        run_name = os.path.basename(folder)
        try:
            results = process_run(folder, run_idx=idx, total_runs=total)
        except Exception as e:
            print(f"  ERROR loading {run_name}: {e}")
            import traceback; traceback.print_exc()
            failed.append(run_name)
            continue

        counts_path = os.path.join(ANALYSIS_DIR, f"{run_name}_cell_counts.csv")
        subs_path = os.path.join(ANALYSIS_DIR, f"{run_name}_substrates.csv")
        final_path = os.path.join(ANALYSIS_DIR, f"{run_name}_final_cells.csv")

        results["counts"].to_csv(counts_path, index=False)
        print(f"  Saved: {counts_path}")
        results["substrates"].to_csv(subs_path, index=False)
        print(f"  Saved: {subs_path}")
        results["final_cells"].to_csv(final_path, index=False)
        print(f"  Saved: {final_path}")

        all_counts.append(results["counts"])
        all_subs.append(results["substrates"])
        succeeded.append(run_name)

        elapsed_total = time.time() - total_start
        avg_per_run = elapsed_total / idx
        remaining = avg_per_run * (total - idx)
        print(f"  Progress: {idx}/{total} runs done | elapsed {elapsed_total:.0f}s | est. remaining {remaining:.0f}s")

    # ── aggregate outputs ──────────────────────────────────────────────────
    print("\n── Writing aggregate files ──────────────────────────────────────")
    if all_counts:
        agg_counts_path = os.path.join(ANALYSIS_DIR, "all_cell_counts.csv")
        pd.concat(all_counts, ignore_index=True).to_csv(agg_counts_path, index=False)
        print(f"  Saved: {agg_counts_path}")
    if all_subs:
        agg_subs_path = os.path.join(ANALYSIS_DIR, "all_substrates.csv")
        pd.concat(all_subs, ignore_index=True).to_csv(agg_subs_path, index=False)
        print(f"  Saved: {agg_subs_path}")

    total_elapsed = time.time() - total_start
    print(f"\n── Summary ──────────────────────────────────────────────────────")
    print(f"  Succeeded: {len(succeeded)}/{total}")
    if failed:
        print(f"  Failed:    {failed}")
    print(f"  Total time: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
    print(f"  Outputs in: {ANALYSIS_DIR}/")


if __name__ == "__main__":
    main()
