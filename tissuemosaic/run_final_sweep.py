"""
run_final_sweep.py
------------------
The confirmed final sweep: FOV 96 vs 128 um, 8 ROIs, 1000 epochs.

    python tissuemosaic/run_final_sweep.py --dry-run     # plan only
    python tissuemosaic/run_final_sweep.py --parallel    # both at once, ~17 h
    python tissuemosaic/run_final_sweep.py               # sequential, ~32 h

Differs from run_window_sweep.py in three ways:

  1. Both configs run CONCURRENTLY under --parallel. Training is single-thread CPU
     bound (one core pinned, GPU ~31%, 14 of 16 cores idle) and peak VRAM is 3.6-3.9
     GB each against 16 GB, so two runs barely contend. ~32 h sequential -> ~17 h.

  2. Each config gets its OWN h5ad folder. AnndataFolderDM caches its raster to
     <data_folder>/train_dataset.pt, and both configs share pixel_size 2.0, so a
     shared folder would have them racing to write the same file.

  3. No --skip-trained equivalent by default: these are long runs and a half-finished
     one should be resumed or removed deliberately, not silently skipped.

Still one subprocess per config, for the same reason as the screen: CUDA context and
Lightning global state do not reset cleanly in-process between runs.
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imc_tm  # noqa: E402  -- applies tm_compat on import

import yaml  # noqa: E402


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--only", nargs="*", default=None, help="train only these config names")
    p.add_argument("--parallel", action="store_true", help="run both configs concurrently")
    p.add_argument("--dry-run", action="store_true", help="write configs, print the plan, train nothing")
    p.add_argument("--run-dir", default=None)
    p.add_argument("--rebuild-data", action="store_true", help="repopulate the per-config h5ad folders")
    return p.parse_args(argv)


# measured on this machine during the S1 screen, at 32 steps/epoch
IT_PER_S = {48: 1.18, 64: 1.07}


def plan(specs, run_dir):
    """Print the plan and return the estimated hours per config."""
    n_roi = len(imc_tm.FINAL_ROIS)
    ncpt = imc_tm.FINAL_OVERRIDES["n_crops_for_tissue_train"]
    bs = imc_tm.FINAL_OVERRIDES["batch_size_per_gpu"]
    epochs = imc_tm.FINAL_OVERRIDES["max_epochs"]
    steps_per_epoch = n_roi * ncpt // bs

    print("{} ROIs x {} crops / batch {} = {} steps/epoch, {} epochs = {:,} steps".format(
        n_roi, ncpt, bs, steps_per_epoch, epochs, steps_per_epoch * epochs))
    print()
    hdr = "{:<14s} {:>6s} {:>4s} {:>8s} {:>10s} {:>9s}"
    print(hdr.format("config", "FOV", "gs", "it/s", "est.hours", "ckpts"))
    hours = {}
    for spec in specs:
        rate = IT_PER_S.get(spec["global_size"], 1.1)
        h = steps_per_epoch * epochs / rate / 3600.0
        hours[spec["name"]] = h
        print(hdr.format(spec["name"], "{}um".format(spec["fov_um"]), str(spec["global_size"]),
                         "{:.2f}".format(rate), "{:.1f}".format(h),
                         str(epochs // imc_tm.FINAL_OVERRIDES["checkpoint_interval_epochs"])))
    print()
    print("run dir: {}".format(run_dir))
    return hours


def write_config(spec, run_dir):
    cfg = imc_tm.final_config_dict(spec)
    d = run_dir / spec["name"]
    d.mkdir(parents=True, exist_ok=True)
    cfg_path = d / "config.yaml"
    with open(cfg_path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return cfg_path, d, Path(cfg["data_folder"])


def clear_caches(data_folder):
    """The raster cache is built for ONE pixel_size; a stale one trains on the wrong
    resolution silently. Always clear before a run."""
    gone = []
    for name in ("train_dataset.pt", "test_dataset.pt"):
        f = Path(data_folder) / name
        if f.is_file():
            f.unlink()
            gone.append(name)
    return gone


def launch(spec, cfg_path, out_dir):
    """Start one training subprocess, unwaited. Returns (Popen, log path).

    The entry point MUST be the local run_train_imc.py shim, not upstream's
    main_1_train_ssl.py. tm_compat is applied by importing imc_tm, and a subprocess
    inherits none of this process's imports -- launching upstream directly means
    weights_only stays True and TissueMosaic dies reading back its own
    train_dataset.pt cache (SparseImage is not an allowed global).
    """
    entry = Path(__file__).resolve().parent / "run_train_imc.py"
    log = out_dir / "train.log"
    fh = open(log, "w")
    proc = subprocess.Popen([sys.executable, "-u", str(entry), "--config", str(cfg_path)],
                            cwd=str(out_dir), stdout=fh, stderr=subprocess.STDOUT)
    proc._log_fh = fh          # keep the handle alive for the process lifetime
    return proc, log


def main(argv):
    args = parse_args(argv)
    run_dir = Path(args.run_dir) if args.run_dir else imc_tm.FINAL_RUN_DIR
    specs = imc_tm.FINAL_SWEEP
    if args.only:
        specs = [s for s in specs if s["name"] in args.only]
        if not specs:
            raise SystemExit("no config matched --only {}".format(args.only))

    print(imc_tm.compat_summary())
    print()
    for d, status in imc_tm.build_final_h5ad_dirs(overwrite=args.rebuild_data):
        print("  data: {:<28s} {}".format(d.name, status))
    print()
    plan(specs, run_dir)

    prepared = []
    for spec in specs:
        cfg_path, out_dir, data_folder = write_config(spec, run_dir)
        gone = clear_caches(data_folder)
        if gone:
            print("  cleared stale cache in {}: {}".format(data_folder.name, ", ".join(gone)))
        prepared.append((spec, cfg_path, out_dir))

    if args.dry_run:
        print("\ndry run -- configs written, nothing trained")
        return 0

    t0 = time.time()
    if args.parallel:
        print("\nlaunching {} config(s) concurrently\n".format(len(prepared)))
        running = []
        for spec, cfg_path, out_dir in prepared:
            proc, log = launch(spec, cfg_path, out_dir)
            print("  [{}] pid {} -> {}".format(spec["name"], proc.pid, log))
            running.append((spec, proc))
        rc_all = 0
        for spec, proc in running:
            rc = proc.wait()
            rc_all |= rc
            print("  [{}] exit {} after {:.2f} h".format(spec["name"], rc, (time.time() - t0) / 3600.0))
        return rc_all

    rc_all = 0
    for i, (spec, cfg_path, out_dir) in enumerate(prepared, 1):
        print("\n[{}/{}] {}".format(i, len(prepared), spec["name"]))
        t1 = time.time()
        proc, log = launch(spec, cfg_path, out_dir)
        print("  log -> {}".format(log))
        rc = proc.wait()
        rc_all |= rc
        print("  exit {} after {:.2f} h".format(rc, (time.time() - t1) / 3600.0))
    return rc_all


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
