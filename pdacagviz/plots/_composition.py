"""composition() -- parts of a whole, per group.

Ported from ``_comp_vis.comp_bar`` and ``div_comp_bar``, which shared a
``get_percents`` front end and differed only in layout. That shared
computation is why ``diverging`` is a ``kind=`` here rather than a second
function.
"""

import numpy as np
import pandas as pd

from .._defaults import MPL_ALIASES, set_default
from ..settings import settings
from ._base import EMPTY, facet, finish, make_axes, pop_fig_kwargs, require_columns, resolve_palette

__all__ = ["composition", "percentages"]

KINDS = ("stacked", "grouped", "diverging")


def percentages(data, group, category, value=None, norm=True):
    """Cross-tabulate into a group x category table.

    With ``value=None`` the rows are counted; otherwise ``value`` is summed.
    ``norm`` turns each group's row into percentages, which is almost always
    what a composition figure wants -- absolute counts hide the comparison
    when group sizes differ.
    """
    if value is None:
        table = pd.crosstab(data[group], data[category])
    else:
        table = data.pivot_table(
            index=group, columns=category, values=value, aggfunc="sum", fill_value=0
        )
    if norm:
        totals = table.sum(axis=1).replace(0, np.nan)
        table = table.div(totals, axis=0) * 100.0
    return table.fillna(0.0)


def composition(
    data,
    x=None,
    hue=None,
    y=None,
    kind="stacked",
    norm=True,
    panel=None,
    palette=None,
    strict=None,
    ax=None,
    plot_kwargs=EMPTY,
    **fig_kwargs,
):
    """Composition of ``hue`` within each level of ``x``.

    Parameters
    ----------
    data
        Long-form frame: one row per observation, or per pre-aggregated cell.
    x
        Grouping column -- one bar, or one bar pair, per level.
    hue
        The composition being broken out, and what the colours encode.
    y
        Optional column to sum. Omit to count rows.
    kind
        ``"stacked"`` for one bar per group, ``"grouped"`` for side by side,
        ``"diverging"`` for two groups mirrored about zero. Diverging requires
        ``x`` to have exactly two levels, since that is what it means.
    norm
        Percentages rather than raw totals.

    Returns the Axes, or a 2-D array of them when ``panel`` is given.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}; expected one of: {', '.join(KINDS)}")
    require_columns(data, x=x, hue=hue, y=y, panel=panel)
    if x is None or hue is None:
        raise ValueError("composition() needs both x= (the groups) and hue= (the parts)")

    fk = pop_fig_kwargs(fig_kwargs)
    fk.setdefault("ylabel" if kind != "diverging" else "xlabel",
                  "Percent of group" if norm else "Count")

    levels = list(pd.unique(data[hue].dropna()))
    colours = resolve_palette(data, hue, palette=palette, strict=strict)

    if kind == "diverging":
        groups = list(pd.unique(data[x].dropna()))
        if len(groups) != 2:
            raise ValueError(
                f"kind='diverging' needs exactly two levels in {x!r}, found {len(groups)}: "
                f"{groups}. Use kind='grouped' for more."
            )

    def draw(subset, cur_ax, title):
        table = percentages(subset, x, hue, value=y, norm=norm)
        table = table.reindex(columns=[c for c in levels if c in table.columns])
        params = dict(plot_kwargs)
        set_default(MPL_ALIASES["edgecolor"], "white", params)
        set_default(MPL_ALIASES["linewidth"], 0.4, params)

        if kind == "diverging":
            _draw_diverging(table, cur_ax, colours, params)
        else:
            table.plot(
                kind="bar",
                stacked=(kind == "stacked"),
                color=[colours[c] for c in table.columns] if colours else None,
                ax=cur_ax,
                legend=False,
                **params,
            )
            cur_ax.set_xlabel(x)
        if title:
            cur_ax.set_title(title)

    result = facet(data, panel, draw, ax=ax, fig_kwargs=fk)

    for cur_ax in np.atleast_1d(result).ravel():
        if cur_ax.figure is not None:
            _legend_proxies(cur_ax, levels, colours)
            finish(cur_ax, fk, hue=hue)
    return result


def _draw_diverging(table, ax, colours, params):
    """Two groups mirrored about zero, one row per category."""
    groups = list(table.index)
    categories = list(table.columns)
    positions = np.arange(len(categories))

    for sign, group in zip((-1.0, 1.0), groups):
        widths = table.loc[group].to_numpy() * sign
        ax.barh(
            positions,
            widths,
            color=[colours[c] for c in categories] if colours else None,
            **params,
        )

    ax.axvline(0.0, color="black", linewidth=0.5, zorder=3)
    ax.set_yticks(positions)
    ax.set_yticklabels(categories)
    ax.xaxis.set_major_formatter(
        __import__("matplotlib").ticker.FuncFormatter(lambda v, _p: f"{abs(v):g}")
    )
    limit = float(np.abs(table.to_numpy()).max()) * 1.08 or 1.0
    ax.set_xlim(-limit, limit)
    ax.text(-limit, len(categories) - 0.4, str(groups[0]), ha="left", va="bottom")
    ax.text(limit, len(categories) - 0.4, str(groups[1]), ha="right", va="bottom")


def _legend_proxies(ax, levels, colours):
    """Give the axes handles to build a legend from.

    pandas' own bar legend is suppressed above so that the entries follow the
    palette order rather than column order, which is what keeps two figures of
    the same data legible side by side.
    """
    if not colours or ax.get_legend_handles_labels()[0]:
        return
    from matplotlib.patches import Patch

    ax.legend(
        handles=[Patch(facecolor=colours[c], label=str(c)) for c in levels if c in colours],
        labels=[str(c) for c in levels if c in colours],
    )
