"""
run_featurize_screen.py
-----------------------
Featurize every config of a screen with main_2_featurize and record what produced it.

    python tissuemosaic/run_featurize_screen.py --screen ecm   --epoch 199
    python tissuemosaic/run_featurize_screen.py --screen final --epoch 574

Two things this adds over calling run_featurize_imc.py by hand.

  1. It writes `_provenance.json` next to the h5ad, naming the checkpoint and epoch. The
     featurized files record no training epoch anywhere, so without this a qualitative
     notebook cannot tell whether its two configs came from the same point in training --
     and comparing different epochs would confound FOV with training length.

  2. It clears the output folder first. main_2_featurize writes per-ROI files and does not
     remove stale ones, so a re-run at a different epoch leaves a folder mixing two
     models, which merge_featurized will happily concatenate.

Channel count is NOT passed here: main_2_featurize rebuilds its datamodule from
`model._hparams` (main_2_featurize.py:157), so an 11-channel ECM checkpoint brings its own
`categories_to_channels`, `pixel_size` and `global_size` with it.
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imc_tm  # noqa: E402  -- applies tm_compat on import

# name -> (specs, checkpoint resolver, source h5ad dir, featurised output dir)
SCREENS = {
    "ecm": (lambda: imc_tm.ECM_SWEEP,
            lambda n, e: imc_tm.ecm_checkpoint(n, e),
            lambda n: imc_tm.ecm_h5ad_dir(n), lambda n: imc_tm.feat_dir(n)),
    "final": (lambda: imc_tm.FINAL_SWEEP,
              lambda n, e: imc_tm.final_checkpoint(n, e),
              lambda n: imc_tm.final_h5ad_dir(n), lambda n: imc_tm.feat_dir(n)),
    # the 1000-epoch ECM test. Same rasters as "ecm", different checkpoints, and a
    # SEPARATE output dir so featurising it cannot delete the 200-epoch screen's data.
    "ecm-long": (lambda: imc_tm.ecm_final_specs(),
                 lambda n, e: imc_tm.ecm_final_checkpoint(n, e),
                 lambda n: imc_tm.ecm_h5ad_dir(n), lambda n: imc_tm.ecm_final_feat_dir(n)),
}


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--screen", choices=sorted(SCREENS), required=True)
    p.add_argument("--epoch", type=int, default=None,
                   help="periodic checkpoint epoch; default is the latest available")
    p.add_argument("--only", nargs="*", default=None)
    p.add_argument("--feature_key", default="dino")
    p.add_argument("--frac_overlap", type=float, default=0.0,
                   help="0.0 keeps windows disjoint; overlapping windows share cells and "
                        "inflate any spatial-coherence metric computed on them")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    specs_fn, ckpt_fn, src_fn, out_fn = SCREENS[args.screen]
    specs = specs_fn()
    if args.only:
        specs = [s for s in specs if s["name"] in args.only]
        if not specs:
            raise SystemExit("no config matched --only {}".format(args.only))

    print(imc_tm.compat_summary())
    print()

    entry = Path(__file__).resolve().parent / "run_featurize_imc.py"
    rc_all = 0
    for spec in specs:
        name = spec["name"]
        ckpt = ckpt_fn(name, args.epoch)
        src, dst = src_fn(name), out_fn(name)
        epoch = args.epoch if args.epoch is not None else -1
        print("{:<13s} fov {:>3d} um".format(name, spec["fov_um"]))
        print("  ckpt  {}".format(ckpt))
        print("  in    {}".format(src))
        print("  out   {}".format(dst))
        if args.dry_run:
            print("  (dry run)\n")
            continue

        if dst.exists():
            shutil.rmtree(dst)      # never merge two epochs into one folder
        dst.mkdir(parents=True)

        t0 = time.time()
        cmd = [sys.executable, "-u", str(entry),
               "--anndata_in", str(src), "--anndata_out", str(dst),
               "--ckpt_in", str(ckpt), "--feature_key", args.feature_key,
               "--patch_strategy", "tiling", "--frac_overlap", str(args.frac_overlap)]
        rc = subprocess.call(cmd)
        rc_all |= rc
        n = len(list(dst.glob("*.h5ad")))
        if rc != 0:
            print("  FAILED rc={} after {:.1f} min\n".format(rc, (time.time() - t0) / 60))
            continue

        (dst / "_provenance.json").write_text(json.dumps({
            "config": name,
            "epoch": int(epoch),
            "checkpoint": str(ckpt.relative_to(imc_tm.REPO)),
            "checkpoint_exists": ckpt.is_file(),
            "n_h5ad": n,
            "featurized_at": datetime.now().replace(microsecond=0).isoformat(),
            "patch_strategy": "tiling",
            "frac_overlap": args.frac_overlap,
            "feature_key": args.feature_key,
        }, indent=1) + "\n")
        print("  wrote {} h5ad + _provenance.json in {:.1f} min\n".format(
            n, (time.time() - t0) / 60))
    return rc_all


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
