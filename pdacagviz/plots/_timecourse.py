"""timecourse() -- a value tracked over time.

New; scanpy_extensions has no equivalent, because a trajectory over simulation
time is not a question single-cell data asks. It is most of what the ABM
output is.
"""

import numpy as np
import seaborn as sns

from .._compat import translate_seaborn
from .._defaults import MPL_ALIASES, set_default
from ._base import EMPTY, facet, finish, pop_fig_kwargs, require_columns, resolve_palette

__all__ = ["timecourse"]


def timecourse(
    data,
    x=None,
    y=None,
    hue=None,
    band="ci",
    panel=None,
    palette=None,
    strict=None,
    ax=None,
    plot_kwargs=EMPTY,
    **fig_kwargs,
):
    """``y`` over ``x``, one line per level of ``hue``.

    ``band`` controls the spread drawn around each line when a time point has
    several observations -- replicate runs, say. ``"ci"`` for a bootstrap
    interval, ``"sd"`` for standard deviation, ``None`` for the line alone.
    With one observation per point there is nothing to draw either way, and
    ``None`` avoids paying for the bootstrap.
    """
    require_columns(data, x=x, y=y, hue=hue, panel=panel)
    if x is None or y is None:
        raise ValueError("timecourse() needs both x= (time) and y= (the value)")

    colours = resolve_palette(data, hue, palette=palette, strict=strict)

    fk = pop_fig_kwargs(fig_kwargs)
    fk.setdefault("ratio", 0.62)
    fk.setdefault("legend", hue is not None)

    def draw(subset, cur_ax, title):
        params = dict(plot_kwargs)
        set_default(MPL_ALIASES["linewidth"], None, params)  # follow rcParams
        set_default("errorbar", band, params)
        set_default("legend", False, params)
        if hue is None:
            set_default(["color", "c"], "black", params)
        params = {k: v for k, v in params.items() if v is not None or k == "errorbar"}
        params, _dropped = translate_seaborn(params)

        sns.lineplot(
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
