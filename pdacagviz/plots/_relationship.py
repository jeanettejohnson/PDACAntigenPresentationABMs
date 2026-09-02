"""relationship() -- how two continuous variables move together.

Ported from ``_feat_vis.rel``. Marker size follows the point count rather than
being fixed, which is the one thing the archived scatter scripts got wrong in
a way that showed: ``s=60`` is right at one figure size and wrong at every
other.
"""

import numpy as np
import seaborn as sns

from .. import axes as _axes
from .._compat import translate_seaborn
from .._defaults import MPL_ALIASES, set_default
from ._base import EMPTY, facet, finish, pop_fig_kwargs, require_columns, resolve_palette

__all__ = ["relationship"]

KINDS = {"scatter": sns.scatterplot, "hist": sns.histplot, "kde": sns.kdeplot}


def relationship(
    data,
    x=None,
    y=None,
    hue=None,
    kind="scatter",
    panel=None,
    palette=None,
    strict=None,
    ax=None,
    plot_kwargs=EMPTY,
    **fig_kwargs,
):
    """Relationship between ``x`` and ``y``.

    ``kind="scatter"`` draws the points, ``"hist"`` bins them, ``"kde"``
    smooths them. Binning and smoothing are the answer when the points overlap
    so heavily that a scatter is a solid blob.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}; expected one of: {', '.join(KINDS)}")
    require_columns(data, x=x, y=y, hue=hue, panel=panel)
    if x is None or y is None:
        raise ValueError("relationship() needs both x= and y=")

    colours = resolve_palette(data, hue, palette=palette, strict=strict)

    fk = pop_fig_kwargs(fig_kwargs)
    fk.setdefault("ratio", 1.0)
    fk.setdefault("legend", hue is not None)

    def draw(subset, cur_ax, title):
        params = dict(plot_kwargs)
        if kind == "scatter":
            set_default(
                MPL_ALIASES["size"],
                _axes.marker_size(max(len(subset), 1), scale=0.5),
                params,
            )
            set_default(MPL_ALIASES["linewidth"], 0.0, params)
            set_default(MPL_ALIASES["edgecolor"], "none", params)
            set_default("legend", False, params)
            if hue is None:
                set_default(["color", "c"], "black", params)
        else:
            set_default("fill", True, params)
            if kind == "hist":
                set_default("bins", 25, params)
            else:
                set_default("bw_method", "silverman", params)
                set_default("cut", 0.5, params)
            if hue is None:
                set_default("cmap", "viridis", params)
        params, _dropped = translate_seaborn(params)

        KINDS[kind](
            data=subset, x=x, y=y, hue=hue, ax=cur_ax,
            palette=colours if (colours and hue is not None) else None,
            **params,
        )
        if title:
            cur_ax.set_title(title)

    result = facet(data, panel, draw, ax=ax, fig_kwargs=fk)
    for cur_ax in np.atleast_1d(result).ravel():
        finish(cur_ax, fk, hue=hue)
    return result
