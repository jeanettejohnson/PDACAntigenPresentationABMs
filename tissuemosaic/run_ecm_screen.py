"""
run_ecm_screen.py
-----------------
The 11-channel ECM screen: FOV 96 vs 128 um, 8 ROIs, 200 epochs.

    python tissuemosaic/run_ecm_screen.py --dry-run     # plan only
    python tissuemosaic/run_ecm_screen.py --parallel    # both at once, ~5.6 h
    python tissuemosaic/run_ecm_screen.py --control     # the missing no-ECM control

Same structure as run_final_sweep.py -- one subprocess per config through the local
run_train_imc.py shim, one h5ad folder per config so the two do not race on
<data_folder>/train_dataset.pt -- with one addition: --control.

WHY --control EXISTS. The screen as run has no matched no-ECM arm, so its results are
compared against final-sweep checkpoints whose LR protocol is (0, 100, 900, 1000) rather
than (0, 20, 180, 200). At epoch 199 the final sweep is at PEAK LR while the ECM runs
have fully annealed, and that difference favours ECM. --control trains 10-channel models
on the ECM screen's exact schedule so the comparison isolates the channel.

It does NOT remove the second difference between the arms. ECM values are all < 1, which
flips CropperSparseTensor to its one-hot branch (dataset.py:365) where n_elements counts
nonzero entries rather than summing cell mass; matrix alone then satisfies
n_element_min_for_crop = 20 at 128 um. The control keeps threshold 10 on cell mass, which
is the final sweep's operating point, so the ECM arm still saw a larger, emptier tile set.
Only a run that deliberately matches the tile sets would close that gap.
"""

import argparse
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
    p.add_argument("--long", action="store_true",
                   help="the 1000-epoch ECM test (warm 100/100) into tissuemosaic/runs/ecm_final")
    p.add_argument("--control", action="store_true",
                   help="10-channel runs on the ECM schedule, into tissuemosaic/runs/ecm_control")
    p.add_argument("--run-dir", default=None)
    p.add_argument("--rebuild-data", action="store_true", help="repopulate the per-config h5ad folders")
    return p.parse_args(argv)


CONTROL_RUN_DIR = imc_tm.RUNS_DIR / "ecm_control"


def control_specs_and_config(spec):
    """A 10-channel config on the ECM screen's schedule, reusing the final-sweep h5ads."""
    cfg = imc_tm.final_config_dict(spec)
    for k in ("max_epochs", "warm_up_epochs", "warm_down_epochs",
              "checkpoint_interval_epochs"):
        cfg[k] = imc_tm.ECM_OVERRIDES[k]
    return cfg


def build_jobs(args):
    jobs = []
    if args.long:
        run_dir = Path(args.run_dir) if args.run_dir else imc_tm.ECM_FINAL_RUN_DIR
        for spec in imc_tm.ecm_final_specs():
            jobs.append(dict(name=spec["name"], spec=spec,
                             cfg=imc_tm.ecm_final_config_dict(spec),
                             data=imc_tm.ecm_h5ad_dir(spec["name"])))
    elif args.control:
        run_dir = Path(args.run_dir) if args.run_dir else CONTROL_RUN_DIR
        for spec in imc_tm.FINAL_SWEEP:
            jobs.append(dict(name=spec["name"].replace("final_", "control_"), spec=spec,
                             cfg=control_specs_and_config(spec),
                             data=imc_tm.final_h5ad_dir(spec["name"])))
    else:
        run_dir = Path(args.run_dir) if args.run_dir else imc_tm.ECM_RUN_DIR
        for spec in imc_tm.ECM_SWEEP:
            jobs.append(dict(name=spec["name"], spec=spec,
                             cfg=imc_tm.ecm_config_dict(spec),
                             data=imc_tm.ecm_h5ad_dir(spec["name"])))
    if args.only:
        jobs = [j for j in jobs if j["name"] in args.only]
        if not jobs:
            raise SystemExit("no config matched --only {}".format(args.only))
    return jobs, run_dir


def prepare(job, run_dir, clear_cache=True):
    d = run_dir / job["name"]
    d.mkdir(parents=True, exist_ok=True)
    cfg_path = d / "config.yaml"
    with open(cfg_path, "w") as fh:
        yaml.safe_dump(job["cfg"], fh, sort_keys=False)
    # The raster cache is built for ONE pixel_size and channel count; a stale one trains
    # on the wrong image silently, so it is always cleared before a real launch. NOT on
    # --dry-run: deleting another run's cache is a side effect, and these folders are
    # shared with the final sweep, which would then rebuild on resume for no reason.
    gone = [n for n in ("train_dataset.pt", "test_dataset.pt")
            if (Path(job["cfg"]["data_folder"]) / n).is_file()]
    if clear_cache:
        for n in gone:
            (Path(job["cfg"]["data_folder"]) / n).unlink()
    return cfg_path, d, gone


def launch(cfg_path, out_dir):
    """Start one training subprocess, unwaited.

    The entry point MUST be the local run_train_imc.py shim, not upstream's
    main_1_train_ssl.py: tm_compat is applied by importing imc_tm, a subprocess inherits
    none of this process's imports, and without it weights_only stays True and
    TissueMosaic dies reading back its own train_dataset.pt cache.
    """
    import subprocess
    entry = Path(__file__).resolve().parent / "run_train_imc.py"
    log = out_dir / "train.log"
    fh = open(log, "w")
    proc = subprocess.Popen([sys.executable, "-u", str(entry), "--config", str(cfg_path)],
                            cwd=str(out_dir), stdout=fh, stderr=subprocess.STDOUT)
    proc._log_fh = fh
    return proc, log


def main(argv):
    args = parse_args(argv)
    jobs, run_dir = build_jobs(args)

    print(imc_tm.compat_summary())
    print()
    if args.control:
        print("CONTROL: 10 channels on the ECM schedule ({} epochs, warm {}/{})".format(
            imc_tm.ECM_OVERRIDES["max_epochs"], imc_tm.ECM_OVERRIDES["warm_up_epochs"],
            imc_tm.ECM_OVERRIDES["warm_down_epochs"]))
    elif args.long:
        o = imc_tm.ECM_FINAL_OVERRIDES
        print("LONG: 11 channels, {} epochs, warm {}/{}, n_element_min_for_crop {}".format(
            o["max_epochs"], o["warm_up_epochs"], o["warm_down_epochs"],
            o["n_element_min_for_crop"]))
        for job in jobs:
            d = Path(job["cfg"]["data_folder"])
            have = sorted(f.stem for f in d.glob("*.h5ad")) if d.is_dir() else []
            if args.rebuild_data or set(have) != set(imc_tm.FINAL_ROIS):
                print("building 11-channel h5ads in {}".format(d.name))
                imc_tm.build_ecm_anndata(imc_tm.FINAL_ROIS, d)
            else:
                print("  data: {:<28s} {} ROIs present".format(d.name, len(have)))
    else:
        for job in jobs:
            d = Path(job["cfg"]["data_folder"])
            have = sorted(f.stem for f in d.glob("*.h5ad")) if d.is_dir() else []
            if args.rebuild_data or set(have) != set(imc_tm.FINAL_ROIS):
                print("building 11-channel h5ads in {}".format(d.name))
                imc_tm.build_ecm_anndata(imc_tm.FINAL_ROIS, d)
            else:
                print("  data: {:<28s} {} ROIs present".format(d.name, len(have)))
    print()

    o = imc_tm.ECM_FINAL_OVERRIDES if args.long else imc_tm.ECM_OVERRIDES
    n_roi = len(imc_tm.FINAL_ROIS)
    steps = n_roi * o["n_crops_for_tissue_train"] // o["batch_size_per_gpu"]
    print("{} ROIs -> {} steps/epoch x {} epochs = {:,} steps per config".format(
        n_roi, steps, o["max_epochs"], steps * o["max_epochs"]))
    print("run dir: {}".format(run_dir))
    print()

    prepared = []
    for job in jobs:
        cfg_path, out_dir, gone = prepare(job, run_dir, clear_cache=not args.dry_run)
        print("  {:<15s} fov {:>3d} um  ch {:>2d}  data {}{}".format(
            job["name"], job["spec"]["fov_um"], job["cfg"].get("image_in_ch", 10),
            Path(job["cfg"]["data_folder"]).name,
            ("  (would clear {})" if args.dry_run else "  (cleared {})").format(", ".join(gone)) if gone else ""))
        prepared.append((job, cfg_path, out_dir))

    if args.dry_run:
        print("\ndry run -- configs written, nothing trained")
        return 0

    t0 = time.time()
    if args.parallel:
        print("\nlaunching {} config(s) concurrently\n".format(len(prepared)))
        running = []
        for job, cfg_path, out_dir in prepared:
            proc, log = launch(cfg_path, out_dir)
            print("  [{}] pid {} -> {}".format(job["name"], proc.pid, log))
            running.append((job, proc))
        rc_all = 0
        for job, proc in running:
            rc = proc.wait()
            rc_all |= rc
            print("  [{}] exit {} after {:.2f} h".format(job["name"], rc, (time.time() - t0) / 3600.0))
        return rc_all

    rc_all = 0
    for i, (job, cfg_path, out_dir) in enumerate(prepared, 1):
        print("\n[{}/{}] {}".format(i, len(prepared), job["name"]))
        t1 = time.time()
        proc, log = launch(cfg_path, out_dir)
        print("  log -> {}".format(log))
        rc = proc.wait()
        rc_all |= rc
        print("  exit {} after {:.2f} h".format(rc, (time.time() - t1) / 3600.0))
    return rc_all


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
