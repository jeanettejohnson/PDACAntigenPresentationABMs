"""
make_imc_anndata.py
-------------------
CLI wrapper around imc_tm.build_anndata(). All the logic lives in imc_tm so the
notebooks and this script cannot drift apart; this file exists so the conversion
stays runnable headless, without a Jupyter kernel:

    python tissuemosaic/make_imc_anndata.py

Converts PhysiCell/user_projects/antigen_presentation/config/ics/JHH_IMC/*.csv
into tissuemosaic/imc_anndata/<ROI>.h5ad. DELETES any existing .h5ad and the
train/test_dataset.pt caches first, so the output never mixes generations.

See imc_tm for the channel definition, the coordinate frame, and the schema.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imc_tm


def main(argv=None):
    print(imc_tm.compat_summary())
    check = imc_tm.check_channels()
    print("channel/yaml consistency: {}".format("OK" if check.get("ok") else check))
    imc_tm.build_anndata(verbose=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
