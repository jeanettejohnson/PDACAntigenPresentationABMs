"""
run_window_sweep.py
-------------------
Trains the four supplementary-experiment S1 models (window size / rasterisation
resolution) on JHH369 only. Long-running -- run it from a terminal, not a
notebook. `S1.window_size_sweep.ipynb` reads the checkpoints afterwards.

    python tissuemosaic/run_window_sweep.py --dry-run     # show the plan
    python tissuemosaic/run_window_sweep.py               # train all four
    python tissuemosaic/run_window_sweep.py --only B_fov192_px2.0
    python tissuemosaic/run_window_sweep.py --skip-trained

Budget: ~13.2 h total on one RTX 5060 Ti (1.6 + 3.6 + 1.6 + 6.4 h), 6400 steps
each. That is a COMPARATIVE SCREEN, not converged models -- every config gets an
identical schedule so the ranking is fair, but no single model is finished.

Each config trains in its own subprocess. That is deliberate: CUDA context and
pytorch-lightning global state do not reset cleanly in-process, and one config's
failure should not take the rest of the sweep with it.

CRITICAL: AnndataFolderDM caches train_dataset.pt / test_dataset.pt next to the
h5ad, and those rasters are built for ONE pixel_size. Reusing them across configs
would silently train on the wrong resolution, with no error. This script deletes
them before every run.
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imc_tm


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--only", nargs="*", default=None, help="train only these config names")
    p.add_argument("--skip-trained", action="store_true", help="skip configs that already have ckpt_last.pt")
    p.add_argument("--dry-run", action="store_true", help="write configs and print the plan, train nothing")
    p.add_argument("--run-dir", default=None)
    return p.parse_args(argv)


def prepare(run_dir):
    """One symlink folder of JHH369 h5ad, plus a config directory per model."""
    dest, linked = imc_tm.sweep_subset_folder()
    print("data folder : {}  ({} ROIs: {})".format(dest, len(linked), ", ".join(sorted(linked))))

    written = []
    for spec in imc_tm.WINDOW_SWEEP:
        d = imc_tm.sweep_run_dir(spec["name"], run_dir)
        d.mkdir(parents=True, exist_ok=True)
        cfg = imc_tm.sweep_config_dict(spec, data_folder=dest)
        path = d / "config.yaml"
        with open(path, "w") as f:
            yaml.safe_dump(cfg, f, sort_keys=False)
        written.append((spec, path, cfg))
    return dest, written


def clear_caches(data_folder):
    """See the CRITICAL note above -- these are pixel_size-specific."""
    removed = [p.name for p in Path(data_folder).glob("*.pt")]
    for p in Path(data_folder).glob("*.pt"):
        p.unlink()
    return removed


def train_one(spec, cfg_path, run_dir, data_folder):
    name = spec["name"]
    d = imc_tm.sweep_run_dir(name, run_dir)
    removed = clear_caches(data_folder)
    if removed:
        print("  cleared pixel_size-specific caches: {}".format(", ".join(removed)))

    entry = Path(__file__).resolve().parent / "run_train_imc.py"
    log = d / "train.log"
    print("  training -> {}".format(d))
    print("  log      -> {}".format(log))
    t0 = time.time()
    with open(log, "w") as fh:
        rc = subprocess.run([sys.executable, str(entry), "--config", str(cfg_path)],
                            cwd=str(d), stdout=fh, stderr=subprocess.STDOUT)
    dt = (time.time() - t0) / 3600.0
    ok = rc.returncode == 0 and imc_tm.sweep_checkpoint(name, run_dir).is_file()
    print("  {} in {:.2f} h (exit {})".format("DONE" if ok else "FAILED", dt, rc.returncode))
    if not ok:
        print("  --- last 15 log lines ---")
        print("\n".join(log.read_text().splitlines()[-15:]))
    return ok, dt


def main(argv):
    args = parse_args(argv)
    print(imc_tm.compat_summary())
    print()

    data_folder, written = prepare(args.run_dir)

    print("\n{:<16s} {:>7s} {:>6s} {:>7s} {:>7s} {:>8s} {:>10s}".format(
        "config", "FOV um", "px", "global", "local", "epochs", "est. hours"))
    # Hours per config, derived from MEASURED runs rather than a compute model:
    # config A (gs=64) took 95 min and B (gs=96) ~126 min at 4000 crops/epoch.
    # n_crops_for_tissue_train=512 halves that to 2048 crops/epoch. D (gs=128) is
    # scaled from B by (128/96)^2. A pure gs^2 model under-predicts the small configs
    # because the single-threaded crop/rasterise path carries fixed overhead.
    # measured: gs64 = 1.64 h; scaling fitted at global_size^0.55 (not ^2 -- the
    # GPU idles behind a single-threaded cropping path, so bigger crops are cheap)
    est = {32: 1.1, 48: 1.4, 64: 1.64, 96: 2.05, 128: 2.61}
    total = 0.0
    plan = []
    for spec, path, cfg in written:
        if args.only and spec["name"] not in args.only:
            continue
        if args.skip_trained and imc_tm.sweep_checkpoint(spec["name"], args.run_dir).is_file():
            print("{:<16s} already trained -- skipping".format(spec["name"]))
            continue
        h = est.get(cfg["global_size"], float("nan"))
        total += h
        plan.append((spec, path, cfg))
        print("{:<16s} {:>7.0f} {:>6.1f} {:>7d} {:>7d} {:>8d} {:>10.1f}".format(
            spec["name"], cfg["pixel_size"] * cfg["global_size"], cfg["pixel_size"],
            cfg["global_size"], cfg["local_size"], cfg["max_epochs"], h))
    print("{:<16s} {:>49s} {:>10.1f}".format("TOTAL", "", total))

    if args.dry_run:
        print("\n--dry-run: configs written, nothing trained.")
        return 0
    if not plan:
        print("\nnothing to do.")
        return 0

    print("\nstarting sweep of {} config(s)\n".format(len(plan)))
    results = []
    for i, (spec, path, cfg) in enumerate(plan, 1):
        print("[{}/{}] {}".format(i, len(plan), spec["name"]))
        results.append((spec["name"], *train_one(spec, path, args.run_dir, data_folder)))
        print()

    print("=" * 56)
    for name, ok, dt in results:
        print("  {:<16s} {:<7s} {:.2f} h".format(name, "ok" if ok else "FAILED", dt))
    print("  total {:.2f} h".format(sum(r[2] for r in results)))
    failed = [r[0] for r in results if not r[1]]
    if failed:
        print("\nFAILED: {}".format(", ".join(failed)))
        return 1
    print("\nAll trained. Next: S1.window_size_sweep.ipynb")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
