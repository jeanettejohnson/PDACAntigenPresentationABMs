"""
tm_compat.py
------------
One compatibility shim, applied on import.

The tissue-mosaic-260819 environment resolves three of the four API breakages
between TissueMosaic and its modern dependency stack natively:

  numpy 1.26.4          restores np.Inf, which pytorch-lightning still calls
  pytorch-lightning 1.9.4  uses neptune.init_run() and exists("sys/id")

What it cannot fix is torch >= 2.6 flipping torch.load's `weights_only` default
to True, because the GPU floors torch at >= 2.7 (an RTX 5060 Ti is sm_120
Blackwell; the oldest conda-forge pytorch with a CUDA 12.8+ build is 2.7.1).
TissueMosaic predates that change and calls torch.load in three places.

This module fixes that WITHOUT disabling the security check globally:

  1. `add_safe_globals` for numpy's scalar/dtype types. Purely additive, so
     Lightning checkpoints load with weights_only=True still ENFORCED.
  2. `weights_only=False` only for the two files TissueMosaic pickles itself,
     in-process, mid-run: train_dataset.pt and test_dataset.pt. Those hold
     SparseImage objects, whose allowlist tail runs through `slice`, pandas and
     anndata internals and does not converge -- so relaxation is unavoidable
     there, and is scoped to exactly those names.

Anything else you torch.load is still refused under weights_only=True; a pickle
carrying an arbitrary callable is rejected as before.

Scoping is by BASENAME, so any file named train_dataset.pt is relaxed wherever
it lives. That is the intended trade: tightening it to a configured data_folder
buys little and couples this module to the datamodule config.

Import it before the first torch.load. There is no longer any constraint
relative to importing pytorch_lightning -- the numpy and neptune shims that
needed that ordering are gone.

Exposes `allowlisted` (count) and `scoped_files` so callers can report what
changed.
"""

import os

import numpy
import numpy.dtypes
import torch

# Files TissueMosaic writes itself in AnndataFolderDM.prepare_data and reads
# straight back in setup(). See datamodule.py:736 and :748.
SCOPED_FILES = frozenset({"train_dataset.pt", "test_dataset.pt"})


def _numpy_safe_globals():
    """The numpy types a Lightning checkpoint's hyper_parameters block contains.

    The module path moved between numpy versions -- 1.26 exposes a `numpy._core`
    stub that has no `multiarray`, so probing with hasattr(numpy, '_core') picks
    the wrong branch. Try both and take whichever actually resolves.
    """
    scalar = None
    for path in ("numpy.core.multiarray", "numpy._core.multiarray"):
        try:
            scalar = getattr(__import__(path, fromlist=["scalar"]), "scalar")
            break
        except Exception:
            continue
    if scalar is None:
        raise RuntimeError("could not resolve multiarray.scalar for numpy {}".format(numpy.__version__))

    dtypes = [getattr(numpy.dtypes, n) for n in dir(numpy.dtypes) if n.endswith("DType")]
    return [scalar, numpy.dtype] + dtypes


def _apply():
    if getattr(torch.load, "_tm_compat_patched", False):
        return 0

    allow = _numpy_safe_globals()
    torch.serialization.add_safe_globals(allow)

    _original_load = torch.load

    def load(f, *args, **kwargs):
        if ("weights_only" not in kwargs
                and isinstance(f, (str, os.PathLike))
                and os.path.basename(str(f)) in SCOPED_FILES):
            kwargs["weights_only"] = False
        return _original_load(f, *args, **kwargs)

    load._tm_compat_patched = True
    load.__doc__ = (
        "tm_compat shim: weights_only relaxed only for {}.\n\n".format(sorted(SCOPED_FILES))
        + (_original_load.__doc__ or "")
    )
    torch.load = load
    return len(allow)


allowlisted = _apply()
scoped_files = sorted(SCOPED_FILES)

__all__ = ["allowlisted", "scoped_files", "SCOPED_FILES"]
