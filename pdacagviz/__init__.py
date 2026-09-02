"""Figure styling and charts for the PDAC antigen-presentation model.

Pandas in, matplotlib out. No anndata, no scanpy.

    import pdacagviz

    pdacagviz.settings.mode = "poster"
    fig, axes = pdacagviz.grid(2, 3, width="half")
    with pdacagviz.settings.using(mode="article"):
        ...

Plot functions live in :mod:`pdacagviz.plots`, one module per data question,
with the mark selected by ``kind=``. Parameter names follow pandas and seaborn
-- ``data``, ``x``, ``y``, ``hue``, ``kind``, ``palette``, ``ax`` -- so the
call sites read the way the rest of the ecosystem does.

Supported floors are as low as python 3.10 allows: matplotlib 3.4,
seaborn 0.11.1, pandas 1.3, numpy 1.21. Differences across that span are
handled in :mod:`pdacagviz._compat` rather than at the call sites.
"""

from ._compat import MPL_VERSION, SNS_VERSION
from ._defaults import MPL_ALIASES, isiterable, set_default
from .modes import MODE_META, MODES
from .plots import (
    bar,
    composition,
    distribution,
    heatmap,
    percentages,
    relationship,
    timecourse,
)
from .palettes import (
    ANTIGEN_CLASS,
    ATLAS,
    GLASBEY_DARK,
    PALETTES,
    PATIENT,
    SIMULATION,
    categorical,
    colors_for,
    normalize,
)
from .axes import (
    bar_labels,
    despine,
    legend_outside,
    marker_size,
    pval_stars,
    thousands,
)
from .export import save, stack_versions
from .settings import Settings, configure, settings
from .sizing import figsize, grid, width_of

__all__ = [
    # state
    "settings",
    "Settings",
    "configure",
    "MODES",
    "MODE_META",
    # plots
    "bar",
    "composition",
    "distribution",
    "relationship",
    "timecourse",
    "heatmap",
    "percentages",
    # sizing
    "figsize",
    "grid",
    "width_of",
    # colour
    "ATLAS",
    "SIMULATION",
    "ANTIGEN_CLASS",
    "PATIENT",
    "PALETTES",
    "GLASBEY_DARK",
    "categorical",
    "colors_for",
    "normalize",
    # axes helpers
    "despine",
    "thousands",
    "legend_outside",
    "bar_labels",
    "marker_size",
    "pval_stars",
    # output
    "save",
    "stack_versions",
    # internals worth exposing
    "set_default",
    "isiterable",
    "MPL_ALIASES",
    "MPL_VERSION",
    "SNS_VERSION",
]

__version__ = "0.1.0.dev0"
