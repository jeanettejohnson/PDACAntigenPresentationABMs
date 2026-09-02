"""Figure dimensions, derived from content rather than typed in.

Two ideas carry this module. First, width comes from a named token -- ``full``,
``half``, ``third`` -- that resolves differently per mode, so the same call
gives a 7-inch column figure in article mode and a 14-inch one on a poster.
Second, height follows from what is being drawn: a bar chart of forty cell
types needs forty rows of vertical space, and that is a calculation rather
than a guess.

The scripts this replaces wrote ``figsize=(9, max(6, n * 0.32 + 1.5))``. The
formula was right; what was missing was that the row pitch has to scale with
the type size, or a poster figure gets poster type crammed into article
spacing.
"""

import matplotlib.pyplot as plt

from .settings import settings

__all__ = ["width_of", "figsize", "grid"]

#: Height:width for a figure with no categorical rows to count.
_DEFAULT_RATIO = 0.75

#: Gaps between panels, as a fraction of one panel. Used to grow the figure so
#: that panels keep their requested size instead of shrinking to fit.
_WSPACE = 5e-2
_HSPACE = 2.5e-2


def width_of(width="full"):
    """Resolve a width token, or pass a number through unchanged.

    ``width`` may be ``"full"``, ``"half"``, ``"third"``, or an explicit
    number of inches for the cases a token does not cover.
    """
    if isinstance(width, (int, float)):
        return float(width)

    widths = settings.meta["widths"]
    if width not in widths:
        known = ", ".join(sorted(widths))
        raise ValueError(f"unknown width {width!r}; expected one of: {known}, or a number")
    return widths[width]


def figsize(width="full", rows=None, ratio=None):
    """Return ``(width, height)`` in inches for the current mode.

    Give ``rows`` for anything with a categorical axis -- bars, a heatmap --
    and the height is computed from the mode's row pitch, so the same call
    stays legible whether there are five categories or fifty. Give ``ratio``
    for continuous plots, where height is a proportion of width. Give neither
    and ``ratio`` defaults to 0.75.

    Passing both is an error rather than a silent precedence rule: they are
    two different intents, and picking one for the caller would hide a
    mistake.
    """
    if rows is not None and ratio is not None:
        raise ValueError(
            "pass rows= for a categorical axis or ratio= for a continuous one, not both"
        )

    w = width_of(width)
    if rows is not None:
        if rows < 0:
            raise ValueError(f"rows must not be negative, got {rows}")
        meta = settings.meta
        h = max(meta["min_height"], rows * meta["row_pitch"] + meta["min_height"])
    else:
        h = w * (_DEFAULT_RATIO if ratio is None else ratio)
    return (w, h)


def grid(nrows=1, ncols=1, width="full", rows=None, ratio=None, panel_size=None, **kwargs):
    """Create a panel grid whose panels are the requested size.

    The figure grows with the panel count rather than dividing a fixed canvas
    between them, so a two-by-three grid has panels the same size as a single
    plot -- which is what keeps type legible when panels are added.

    ``panel_size`` overrides the per-panel dimensions directly. Otherwise they
    come from ``width``, ``rows`` and ``ratio`` exactly as in :func:`figsize`.

    Returns ``(fig, axes)`` with ``axes`` always 2-D, so indexing does not
    change shape when a grid happens to have one row.
    """
    if nrows < 1 or ncols < 1:
        raise ValueError(f"grid needs at least one row and column, got {nrows}x{ncols}")

    pw, ph = panel_size if panel_size is not None else figsize(width, rows=rows, ratio=ratio)
    fig_w = ((1 + _WSPACE) * ncols - _WSPACE) * pw
    fig_h = ((1 + _HSPACE) * nrows - _HSPACE) * ph

    kwargs.setdefault("squeeze", False)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_w, fig_h), **kwargs)
    return fig, axes
