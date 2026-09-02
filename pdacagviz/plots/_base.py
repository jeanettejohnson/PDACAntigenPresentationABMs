"""Shared machinery behind every plot function.

The signature every function takes::

    f(data, x=, y=, hue=, kind=, panel=None, palette=None, ax=None,
      plot_kwargs={}, **fig_kwargs)

``plot_kwargs`` reaches the underlying seaborn or matplotlib call untouched.
``**fig_kwargs`` configures the figure and axes. That two-tier split is the
best idea carried over from scanpy_extensions: it lets a function be fully
styled by default and still fully overridable, without either concern leaking
into the other.
"""

from types import MappingProxyType

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import axes as _axes
from ..palettes import categorical, colors_for, normalize
from ..settings import settings
from ..sizing import figsize as _figsize
from ..sizing import grid as _grid

__all__ = [
    "EMPTY",
    "resolve_palette",
    "make_axes",
    "finish",
    "facet",
    "require_columns",
    "pop_fig_kwargs",
]

#: Shared immutable default for the ``*_kwargs`` parameters. A fresh ``{}`` per
#: call signature would be a mutable default; this cannot be mutated at all.
EMPTY = MappingProxyType({})

#: Everything ``**fig_kwargs`` accepts. Anything else is a typo, and is
#: reported as one rather than silently ignored -- a misspelled styling option
#: that does nothing is the kind of thing found only after a figure is printed.
FIG_KWARGS = frozenset(
    {
        "width", "rows", "ratio", "figsize", "panel_size",
        "title", "xlabel", "ylabel",
        "xlim", "ylim", "x_rotation", "y_rotation",
        "despine", "legend", "legend_title", "legend_ncols",
        "thousands", "grid",
    }
)


def pop_fig_kwargs(kwargs):
    """Split figure options out of ``**kwargs``, rejecting unknown names."""
    unknown = set(kwargs) - FIG_KWARGS
    if unknown:
        known = ", ".join(sorted(FIG_KWARGS))
        raise TypeError(
            f"unknown option(s) {', '.join(map(repr, sorted(unknown)))}. "
            f"Figure options are: {known}. "
            "Anything meant for the underlying plotting call goes in plot_kwargs."
        )
    return dict(kwargs)


def require_columns(data, **named):
    """Check that named columns exist, naming all the misses at once."""
    if not isinstance(data, pd.DataFrame):
        raise TypeError(f"data must be a DataFrame, got {type(data).__name__}")
    missing = {role: col for role, col in named.items() if col is not None and col not in data}
    if missing:
        detail = ", ".join(f"{role}={col!r}" for role, col in sorted(missing.items()))
        raise KeyError(
            f"column(s) not in data: {detail}. Available: {', '.join(map(str, data.columns))}"
        )


def resolve_palette(data, key, palette=None, strict=None):
    """Colours for the categories in ``data[key]``, as a dict.

    Returns ``None`` when there is no categorical key to colour by, which the
    callers treat as "let the mark use a single colour".

    Strict lookup applies only to columns that are plausibly cell types. A
    column none of whose values appear in the palette is a different kind of
    column -- ``condition``, ``sample``, ``treatment`` -- and gets matplotlib's
    default cycle instead. Raising there would force every figure coloured by
    treatment arm to pass ``strict=False``, which would in turn disarm the
    check for the cell-type figures that need it.

    A column where *some* values match is the case worth failing on: that is a
    cell-type column with entries missing from the palette, which is exactly
    how four agent types came to be drawn in silent grey.
    """
    if key is None:
        return None

    levels = list(pd.unique(data[key].dropna()))
    if not levels:
        return None

    table = _palette_table(palette)
    known = sum(1 for level in levels if level in table or normalize(level) in table)

    if known == 0:
        # Glasbey rather than matplotlib's ten-colour cycle: the columns that land
        # here are conditions, arms and patients, and this data has 21 of the last
        # one. A ten-colour cycle silently repeats itself at the eleventh group,
        # which reads as two categories being the same thing.
        return dict(zip(levels, categorical(len(levels))))

    return colors_for(levels, palette=palette, strict=strict, as_dict=True)


def _palette_table(palette):
    """The mapping a palette argument refers to, without validating it here."""
    from ..palettes import PALETTES

    if palette is None:
        palette = settings.palette
    if isinstance(palette, str):
        return PALETTES.get(palette, {})
    return dict(palette)


def make_axes(ax=None, fig_kwargs=EMPTY):
    """Return ``(fig, ax)``, creating a figure only when one was not given.

    Size comes from the mode's tokens unless ``figsize`` is passed outright,
    so a call that says nothing about size still gets one appropriate to
    article or poster.
    """
    if ax is not None:
        return ax.get_figure(), ax

    explicit = fig_kwargs.get("figsize")
    if explicit is not None:
        size = explicit
    else:
        size = _figsize(
            fig_kwargs.get("width", "full"),
            rows=fig_kwargs.get("rows"),
            ratio=fig_kwargs.get("ratio"),
        )
    fig, grid_axes = _grid(1, 1, panel_size=size)
    return fig, grid_axes[0, 0]


def finish(ax, fig_kwargs=EMPTY, hue=None):
    """Apply the figure options that every plot function shares.

    Called last so that anything set here wins over the underlying plotting
    call's own defaults -- seaborn in particular likes to name axes after the
    dataframe columns, which is rarely what a figure caption wants.
    """
    fk = fig_kwargs

    if fk.get("title") is not None:
        ax.set_title(fk["title"])
    if fk.get("xlabel") is not None:
        ax.set_xlabel(fk["xlabel"])
    if fk.get("ylabel") is not None:
        ax.set_ylabel(fk["ylabel"])
    if fk.get("xlim") is not None:
        ax.set_xlim(fk["xlim"])
    if fk.get("ylim") is not None:
        ax.set_ylim(fk["ylim"])

    if fk.get("x_rotation") is not None:
        for label in ax.get_xticklabels():
            label.set_rotation(fk["x_rotation"])
            label.set_ha("right" if fk["x_rotation"] % 360 not in (0, 180) else "center")
    if fk.get("y_rotation") is not None:
        for label in ax.get_yticklabels():
            label.set_rotation(fk["y_rotation"])

    which = fk.get("thousands")
    if which:
        for name in which if isinstance(which, (list, tuple)) else [which]:
            _axes.thousands(getattr(ax, f"{name}axis"))

    if fk.get("grid"):
        ax.grid(True, axis=fk["grid"] if isinstance(fk["grid"], str) else "both")

    if fk.get("despine", True):
        _axes.despine(ax)

    legend = fk.get("legend", hue is not None)
    if legend:
        _axes.legend_outside(
            ax, title=fk.get("legend_title", hue), ncols=fk.get("legend_ncols", 1)
        )
    elif ax.get_legend() is not None:
        ax.get_legend().remove()

    return ax


def facet(data, panel, draw, ax=None, fig_kwargs=EMPTY):
    """Draw one panel per level of ``panel``, or a single panel when it is None.

    ``draw`` is called as ``draw(subset, ax, title)``. Returns a single Axes
    for the unfacetted case and a 2-D array otherwise, matching what callers
    expect from ``ax=`` having been honoured.
    """
    if panel is None:
        fig, single = make_axes(ax, fig_kwargs)
        draw(data, single, fig_kwargs.get("title"))
        return single

    if ax is not None:
        raise ValueError("pass either ax= for a single panel or panel= to facet, not both")

    levels = list(pd.unique(data[panel].dropna()))
    ncols = fig_kwargs.get("legend_ncols") if False else min(len(levels), 4)
    nrows = int(np.ceil(len(levels) / ncols))

    explicit = fig_kwargs.get("figsize")
    size = explicit if explicit is not None else _figsize(
        fig_kwargs.get("width", "half"),
        rows=fig_kwargs.get("rows"),
        ratio=fig_kwargs.get("ratio"),
    )
    fig, grid_axes = _grid(nrows, ncols, panel_size=size)

    for index, level in enumerate(levels):
        cur = grid_axes[index // ncols][index % ncols]
        draw(data[data[panel] == level], cur, str(level))
    for index in range(len(levels), nrows * ncols):
        grid_axes[index // ncols][index % ncols].remove()

    return grid_axes
