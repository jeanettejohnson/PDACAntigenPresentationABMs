"""bar() -- one precomputed value per category.

Distinct from ``distribution(kind="bar")``, which aggregates many observations
into a mean with an error bar. Here each row is already the number to draw,
which is what the archived agent-count figures needed.
"""

import numpy as np
import pandas as pd

from .. import axes as _axes
from .._defaults import MPL_ALIASES, set_default
from ._base import EMPTY, facet, finish, pop_fig_kwargs, require_columns, resolve_palette

__all__ = ["bar"]


def bar(
    data,
    x=None,
    y=None,
    hue=None,
    orient="h",
    sort=True,
    labels=True,
    label_format="{:,}",
    panel=None,
    palette=None,
    strict=None,
    ax=None,
    plot_kwargs=EMPTY,
    **fig_kwargs,
):
    """Bar chart of one value per category.

    Parameters
    ----------
    x, y
        Category and value. With ``orient="h"`` the categories run down the
        y axis, so pass the category as ``y`` -- or as ``x`` and let it be
        swapped, since both readings are natural for a horizontal bar.
    hue
        Column supplying colours. Defaults to the category column, which is
        what makes cell-type colours consistent without being asked for.
    sort
        Order bars by value. Long category lists are unreadable unsorted.
    labels
        Write each value at the end of its bar.

    Returns the Axes, or a 2-D array of them when ``panel`` is given.
    """
    if orient not in ("h", "v"):
        raise ValueError(f"orient must be 'h' or 'v', got {orient!r}")
    require_columns(data, x=x, y=y, hue=hue, panel=panel)

    category, value = (y, x) if orient == "h" else (x, y)
    if category is None or value is None:
        raise ValueError("bar() needs both a category column and a value column")
    if not pd.api.types.is_numeric_dtype(data[value]):
        raise TypeError(
            f"the value column {value!r} is not numeric. For counts of rows, "
            "aggregate first or use composition()."
        )

    colour_key = hue if hue is not None else category
    colours = resolve_palette(data, colour_key, palette=palette, strict=strict)

    fk = pop_fig_kwargs(fig_kwargs)
    fk.setdefault("rows", int(data[category].nunique()))
    fk.setdefault("thousands", "x" if orient == "h" else "y")
    fk.setdefault("legend", hue is not None)

    def draw(subset, cur_ax, title):
        frame = subset.sort_values(value, ascending=True) if sort else subset
        params = dict(plot_kwargs)
        set_default(MPL_ALIASES["edgecolor"], "white", params)
        set_default(MPL_ALIASES["linewidth"], 0.4, params)
        colour_list = (
            [colours[k] for k in frame[colour_key]] if colours else None
        )
        plot = cur_ax.barh if orient == "h" else cur_ax.bar
        bars = plot(
            frame[category].astype(str), frame[value], color=colour_list, **params
        )
        if labels:
            _axes.bar_labels(
                cur_ax, bars, frame[value], fmt=label_format, horizontal=(orient == "h")
            )
        cur_ax.set_xlabel(value if orient == "h" else category)
        cur_ax.set_ylabel(category if orient == "h" else value)
        if title:
            cur_ax.set_title(title)

    result = facet(data, panel, draw, ax=ax, fig_kwargs=fk)
    for cur_ax in np.atleast_1d(result).ravel():
        finish(cur_ax, fk, hue=hue)
    return result
