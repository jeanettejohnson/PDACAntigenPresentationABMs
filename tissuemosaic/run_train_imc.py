"""
run_train_imc.py
----------------
Applies the compatibility shim, then hands off to TissueMosaic's own training
entry point, run/main_1_train_ssl.py.

    mkdir -p tissuemosaic/runs/dino_imc && cd tissuemosaic/runs/dino_imc
    python ../../tissuemosaic/run_train_imc.py --config ../../tissuemosaic/config_dino_ssl_imc.yaml

Arguments are forwarded verbatim. Note main_1_train_ssl.py resolves parameters
in the order  yaml > CLI > defaults, so a value present in the yaml cannot be
overridden from the command line -- copy the yaml to vary one.

Checkpoints and the offline Neptune log land in the CURRENT WORKING DIRECTORY,
so run this from a dedicated output directory.
"""

import runpy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imc_tm  # applies tm_compat on import

ENTRY = imc_tm.tissuemosaic_run_script("main_1_train_ssl.py")
print(imc_tm.compat_summary())
print("running {} with args {}".format(ENTRY, sys.argv[1:]))

sys.argv = [str(ENTRY)] + sys.argv[1:]
runpy.run_path(str(ENTRY), run_name="__main__")
