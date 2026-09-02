"""heatmap() -- a value across two categorical axes.

Ported from ``_feat_aggr.aggr``, which is the most anndata-entangled of the
sources: ``AggrFigure`` carried its own figure creation, legend, main and
annotation methods. Only the drawing logic survives here; the rest was
scaffolding for threading an AnnData through scanpy.

``kind="dot"`` is scanpy's dotplot: colour still encodes the value, and marker
area encodes a second quantity -- how many cells contributed, typically -- so
a strong mean resting on three cells does not read like one resting on three
thousand.
"""

import numpy as np
import pandas as pd

from .. import axes as _axes
from .._defaults import set_default
from ._base import EMPTY, finish, make_axes, pop_fig_kwargs, require_columns

__all__ = ["heatmap"]

KINDS = ("matrix", "dot")


def heatmap(
    data,
    x=None,
    y=None,
    value=None,
    size=None,
    kind="matrix",
    cmap=None,
    vcenter=None,
    colorbar=True,
    ax=None,
    plot_kwargs=EMPTY,
    **fig_kwargs,
):
    """``value`` across the grid formed by ``x`` and ``y``.

    Parameters
    ----------
    x, y
        The two categorical axes.
    value
        What colour encodes. Averaged where a cell has several rows.
    size
        For ``kind="dot"``, what marker area encodes. Defaults to the number
        of rows behind each cell, which is the count a dotplot usually shows.
    vcenter
        Anchor the colour scale at this value, for a diverging quantity where
        zero should read as neutral rather than as mid-range.

    Returns the Axes.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}; expected one of: {', '.join(KINDS)}")
    require_columns(data, x=x, y=y, value=value, size=size)
    if x is None or y is None or value is None:
        raise ValueError("heatmap() needs x=, y= and value=")

    grid = data.pivot_table(index=y, columns=x, values=value, aggfunc="mean")
    rows, cols = list(grid.index), list(grid.columns)

    fk = pop_fig_kwargs(fig_kwargs)
    fk.setdefault("rows", len(rows))
    fk.setdefault("x_rotation", 90)
    fk.setdefault("despine", False)
    fk.setdefault("legend", False)
    fig, cur_ax = make_axes(ax, fk)

    norm = None
    if vcenter is not None:
        import matplotlib.colors as mcolors

        span = float(np.nanmax(np.abs(grid.to_numpy() - vcenter)))
        norm = mcolors.TwoSlopeNorm(
            vmin=vcenter - span, vcenter=vcenter, vmax=vcenter + span
        )

    params = dict(plot_kwargs)
    set_default("cmap", cmap or "viridis", params)
    if norm is not None:
        set_default("norm", norm, params)

    if kind == "matrix":
        mappable = cur_ax.pcolormesh(
            np.arange(len(cols) + 1), np.arange(len(rows) + 1), grid.to_numpy(), **params
        )
        offset = 0.5
    else:
        if size is None:
            counts = data.pivot_table(index=y, columns=x, values=value, aggfunc="size")
        else:
            counts = data.pivot_table(index=y, columns=x, values=size, aggfunc="mean")
        counts = counts.reindex(index=rows, columns=cols).fillna(0.0).to_numpy()
        largest = counts.max() or 1.0
        area = _axes.marker_size(max(len(rows) * len(cols), 1))
        xs, ys = np.meshgrid(np.arange(len(cols)), np.arange(len(rows)))
        params.pop("norm", None) if norm is None else None
        mappable = cur_ax.scatter(
            xs.ravel(), ys.ravel(),
            c=grid.to_numpy().ravel(),
            s=(counts.ravel() / largest) * area,
            **params,
        )
        offset = 0.0

    cur_ax.set_xticks(np.arange(len(cols)) + offset)
    cur_ax.set_xticklabels([str(c) for c in cols])
    cur_ax.set_yticks(np.arange(len(rows)) + offset)
    cur_ax.set_yticklabels([str(r) for r in rows])
    cur_ax.set_xlabel(x)
    cur_ax.set_ylabel(y)
    if kind == "dot":
        cur_ax.set_xlim(-0.5, len(cols) - 0.5)
        cur_ax.set_ylim(-0.5, len(rows) - 0.5)

    if colorbar:
        bar = fig.colorbar(mappable, ax=cur_ax, fraction=4.6e-2, pad=2e-2)
        bar.set_label(value)
        bar.outline.set_linewidth(0.0)

    return finish(cur_ax, fk)
