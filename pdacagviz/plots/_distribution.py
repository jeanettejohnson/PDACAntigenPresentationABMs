"""distribution() -- how a value is spread across categories.

Ported from ``_feat_vis.dis``. The seaborn keyword differences that module
branched on inline are handled once in :mod:`pdacagviz._compat` instead.
"""

import numpy as np
import seaborn as sns

from .._compat import translate_seaborn
from .._defaults import MPL_ALIASES, set_default
from ._base import EMPTY, facet, finish, pop_fig_kwargs, require_columns, resolve_palette

__all__ = ["distribution"]

KINDS = {"violin": sns.violinplot, "box": sns.boxplot, "bar": sns.barplot}


def distribution(
    data,
    x=None,
    y=None,
    hue=None,
    kind="violin",
    panel=None,
    palette=None,
    strict=None,
    ax=None,
    plot_kwargs=EMPTY,
    **fig_kwargs,
):
    """Distribution of ``y`` across the categories in ``x``.

    ``kind`` selects the mark: ``"violin"`` shows the full shape, ``"box"``
    the quartiles, ``"bar"`` a mean with an error bar. Swapping between them
    is a one-word edit precisely because the surrounding call does not change.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}; expected one of: {', '.join(KINDS)}")
    require_columns(data, x=x, y=y, hue=hue, panel=panel)

    colour_key = hue if hue is not None else x
    colours = resolve_palette(data, colour_key, palette=palette, strict=strict)

    fk = pop_fig_kwargs(fig_kwargs)
    fk.setdefault("x_rotation", 90)
    fk.setdefault("legend", hue is not None)

    def draw(subset, cur_ax, title):
        params = dict(plot_kwargs)
        set_default(MPL_ALIASES["linewidth"], 0.5, params)
        if kind == "bar":
            set_default("errorbar", None, params)
        elif kind == "violin":
            set_default("density_norm", "width", params)
            set_default("bw_method", "silverman", params)
            set_default("inner", None, params)
            set_default("cut", 1.5, params)
        # seaborn 0.13 wants a hue to accompany a palette; passing the category
        # column satisfies that without changing what is drawn.
        set_default("hue", colour_key, params)
        set_default("legend", False, params)
        params, _dropped = translate_seaborn(params)

        KINDS[kind](
            data=subset, x=x, y=y, ax=cur_ax,
            palette=colours if colours else None, **params
        )
        if cur_ax.get_legend() is not None:
            cur_ax.get_legend().remove()
        if title:
            cur_ax.set_title(title)

    result = facet(data, panel, draw, ax=ax, fig_kwargs=fk)
    for cur_ax in np.atleast_1d(result).ravel():
        finish(cur_ax, fk, hue=hue)
    return result
