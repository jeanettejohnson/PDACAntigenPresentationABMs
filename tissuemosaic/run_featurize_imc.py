"""
run_featurize_imc.py
--------------------
Applies the compatibility shim, then hands off to TissueMosaic's featurization
entry point, run/main_2_featurize.py.

    python tissuemosaic/run_featurize_imc.py \
        --anndata_in  tissuemosaic/imc_anndata \
        --anndata_out tissuemosaic/imc_anndata_featurized \
        --ckpt_in     tissuemosaic/runs/dino_imc/ckpt_last.pt \
        --feature_key dino --n_patches 400 --ncv_k 25 100

Unlike main_1 this has no yaml config -- everything is a CLI flag, and
--anndata_in / --anndata_out / --ckpt_in / --feature_key are required. Writing
to a different --anndata_out than --anndata_in is safer: the featurized h5ad
carry the full sparse-image state and are much larger.
"""

import runpy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imc_tm  # applies tm_compat on import

ENTRY = imc_tm.tissuemosaic_run_script("main_2_featurize.py")
print(imc_tm.compat_summary())

# main_2 writes into --anndata_out but never creates it
for i, a in enumerate(sys.argv):
    if a == "--anndata_out" and i + 1 < len(sys.argv):
        Path(sys.argv[i + 1]).mkdir(parents=True, exist_ok=True)

print("running {} with args {}".format(ENTRY, sys.argv[1:]))
sys.argv = [str(ENTRY)] + sys.argv[1:]
runpy.run_path(str(ENTRY), run_name="__main__")
