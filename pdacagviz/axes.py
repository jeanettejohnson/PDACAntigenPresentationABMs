"""The lines every archived script wrote by hand.

Nothing here is invented: each function is something that appeared in at least
two of the eleven scripts being retired, or in scanpy_extensions' ``_helper``.
The bar for adding to this module is that a call site already exists twice --
otherwise it belongs in the script that needs it.
"""

import numpy as np
from matplotlib.ticker import FuncFormatter

__all__ = [
    "despine",
    "thousands",
    "legend_outside",
    "bar_labels",
    "marker_size",
    "pval_stars",
    "PVAL_THRESHOLDS",
]

#: Cumulative: a p-value at or below each threshold earns another mark.
PVAL_THRESHOLDS = (1e-3, 1e-2, 5e-2)

#: U+204E LOW ASTERISK. Sits on the baseline rather than riding high like a
#: footnote marker, which is what keeps a row of them level with the tick text.
_STAR = "⁎"
_NOT_SIGNIFICANT = "ns"

#: Bar-label offset, as a fraction of the axis data range.
_LABEL_OFFSET = 5e-3

#: Headroom past the longest bar so its label is not clipped.
_LABEL_HEADROOM = 1.12


def despine(ax, top=True, right=True, left=False, bottom=False):
    """Hide the spines named True. Top and right by default.

    Uses list indexing on ``ax.spines``, which is matplotlib 3.4 and up -- the
    same floor the package declares.
    """
    drop = [
        name
        for name, wanted in (("top", top), ("right", right), ("left", left), ("bottom", bottom))
        if wanted
    ]
    if drop:
        ax.spines[drop].set_visible(False)
    return ax


def thousands(axis, decimals=0):
    """Format an axis with thousands separators: 12345 renders as 12,345.

    Pass an axis, not an Axes -- ``thousands(ax.xaxis)``. Being explicit about
    which one avoids the guessing that a bare Axes would require, and reads
    correctly for the horizontal bar charts where only x needs it.
    """
    fmt = f"{{x:,.{decimals}f}}"
    axis.set_major_formatter(FuncFormatter(lambda x, _pos: fmt.format(x=x)))
    return axis


def legend_outside(ax, title=None, ncols=1, **kwargs):
    """Place the legend outside the axes, top-left aligned to the right edge.

    Keeping it out of the data area matters more at article size than on
    screen: at 1.75 inches wide there is no interior space a legend can occupy
    without covering something.
    """
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None

    kwargs.setdefault("loc", "upper left")
    kwargs.setdefault("bbox_to_anchor", (1.02, 1.01))
    kwargs.setdefault("frameon", False)
    kwargs.setdefault("ncols" if _HAS_NCOLS else "ncol", ncols)
    return ax.legend(handles, labels, title=title, **kwargs)


def _detect_ncols_kwarg():
    """``ncols`` replaced ``ncol`` on legends in matplotlib 3.6."""
    import matplotlib as mpl
    from packaging.version import Version

    return Version(mpl.__version__) >= Version("3.6")


_HAS_NCOLS = _detect_ncols_kwarg()


def bar_labels(ax, bars, values, fmt="{:,}", horizontal=True, headroom=True, **kwargs):
    """Write each bar's value at its far end.

    ``fmt`` is applied to every value, so pass ``"{:.1%}"`` for proportions.
    With ``headroom``, the value axis is extended so the longest label is not
    clipped -- the archived scripts each set that limit by hand, and one of
    them got it wrong.
    """
    values = list(values)
    if len(values) != len(bars):
        raise ValueError(f"{len(bars)} bars but {len(values)} values")
    if not values:
        return ax

    span = max(abs(v) for v in values) or 1.0
    offset = span * _LABEL_OFFSET
    kwargs.setdefault("va", "center" if horizontal else "bottom")
    kwargs.setdefault("ha", "left" if horizontal else "center")

    for bar, value in zip(bars, values):
        if horizontal:
            x = bar.get_width() + offset
            y = bar.get_y() + bar.get_height() / 2
        else:
            x = bar.get_x() + bar.get_width() / 2
            y = bar.get_height() + offset
        ax.text(x, y, fmt.format(value), **kwargs)

    if headroom:
        limit = span * _LABEL_HEADROOM
        if horizontal:
            ax.set_xlim(min(0, min(values)), limit)
        else:
            ax.set_ylim(min(0, min(values)), limit)
    return ax


def marker_size(n, figsize=None, scale=1.0):
    """Marker area that stays readable as the point count changes.

    Ported from ``scanpy_extensions._helper``. Size falls as the square root
    of the count and rises with the drawing area, with a floor so that a
    handful of points on a small figure do not become invisible. The archived
    scatter scripts used a fixed ``s=60``, which is right for one figure size
    and wrong for every other.
    """
    import matplotlib.pyplot as plt

    if n <= 0:
        raise ValueError(f"marker_size needs at least one point, got {n}")

    figsize = figsize if figsize is not None else plt.rcParams["figure.figsize"]
    area = figsize[0] * figsize[1]
    fontsize = plt.rcParams["font.size"]

    floor = fontsize / area
    scaled = (fontsize * 10.0 * area) / np.sqrt(n)
    return max(floor, scaled) * scale


def pval_stars(pval, as_stars=True):
    """Render a p-value as significance marks, or as scientific notation.

    Cumulative against :data:`PVAL_THRESHOLDS`, so 0.04 gives one mark, 0.004
    two, and 0.0004 three. Anything above the largest threshold is ``"ns"``,
    spelled out rather than left blank so a reader can tell a
    non-significant result from a missing one.
    """
    if not as_stars:
        return f"{pval:.2e}"
    marks = _STAR * sum(pval <= t for t in PVAL_THRESHOLDS)
    return marks or _NOT_SIGNIFICANT
