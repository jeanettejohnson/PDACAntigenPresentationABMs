"""The article and poster rcParam tables.

A mode fixes type scale, line weights, tick geometry, and save resolution --
everything about how a figure is drawn at a given physical size. It does not
fix the size itself: ``figure.figsize`` is deliberately absent from these
tables, because figure dimensions here follow from content (how many
categorical rows, how many panels) and belong to :mod:`pdacagviz.sizing`.

Tables may name rcParam keys from any matplotlib version;
:func:`pdacagviz._compat.filter_rcparams` drops what the installed one does not
recognise. ``figure.labelsize`` and ``figure.labelweight`` are here for exactly
that reason -- they arrived in 3.6 and are skipped below it.
"""

__all__ = ["MODES", "MODE_META", "resolve"]

_GRID_COLOR = "#ababab"

#: Shared by both modes. Anything a mode overrides is repeated in its own
#: table below rather than mutated here, so each mode reads as a whole.
_BASE = {
    # Lines and patches
    "lines.linewidth": 0.5,
    "lines.markeredgecolor": "none",
    "lines.markeredgewidth": 0.0,
    "patch.linewidth": 0.5,
    # Type
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.0,
    "mathtext.fontset": "dejavuserif",
    "mathtext.default": "regular",
    # Axes
    "axes.linewidth": 0.5,
    "axes.titlesize": 7.0,
    "axes.titleweight": "normal",
    "axes.titlepad": 4.0,
    "axes.labelsize": 7.0,
    "axes.labelpad": 2.0,
    "axes.labelweight": "normal",
    "axes.axisbelow": True,
    "axes.xmargin": 2.5e-2,
    "axes.ymargin": 2.5e-2,
    # Ticks
    "xtick.major.size": 2.0,
    "xtick.major.width": 0.5,
    "xtick.major.pad": 2.0,
    "xtick.labelsize": 6.0,
    "ytick.major.size": 2.0,
    "ytick.major.width": 0.5,
    "ytick.major.pad": 2.0,
    "ytick.labelsize": 6.0,
    # Grid
    "grid.color": _GRID_COLOR,
    "grid.linewidth": 0.5,
    # Legend
    "legend.edgecolor": _GRID_COLOR,
    "legend.markerscale": 0.5,
    "legend.fontsize": 6.0,
    "legend.title_fontsize": 6.0,
    "legend.borderpad": 1 / 3,
    "legend.handlelength": 1.0,
    "legend.handleheight": 0.5,
    "legend.handletextpad": 0.5,
    "legend.columnspacing": 1.5,
    # Figure
    "figure.titlesize": 9.0,
    "figure.titleweight": "bold",
    "figure.labelsize": 7.0,       # matplotlib >= 3.6, filtered below it
    "figure.labelweight": "normal",  # matplotlib >= 3.6, filtered below it
    # Layout. constrained_layout owns the outer margin, so savefig.bbox is
    # left alone -- setting it to "tight" would re-crop after constrained
    # layout resolved, making savefig.pad_inches the real margin and these
    # pads partly decorative.
    "figure.constrained_layout.use": True,
    "figure.constrained_layout.h_pad": 1 / 36,
    "figure.constrained_layout.w_pad": 1 / 36,
    "figure.constrained_layout.hspace": 1e-3,
    "figure.constrained_layout.wspace": 1e-3,
    # Marks
    "scatter.edgecolors": "none",
    # Output. fonttype 42 embeds real fonts, which is what keeps text
    # editable in Illustrator instead of arriving as outlines.
    "savefig.transparent": True,
    "ps.fonttype": 42,
    "ps.useafm": False,
    "pdf.fonttype": 42,
    "pdf.use14corefonts": False,
    "pdf.inheritcolor": False,
}

ARTICLE = dict(_BASE)
ARTICLE.update(
    {
        "figure.dpi": 150.0,
        "savefig.dpi": 300.0,
    }
)

POSTER = dict(_BASE)
POSTER.update(
    {
        "font.sans-serif": ["DejaVu Sans"],
        "font.stretch": "condensed",
        "font.size": 18.0,
        "axes.titlesize": 18.0,
        "axes.titleweight": "normal",
        "axes.titlepad": 7.5,
        "axes.labelsize": 18.0,
        "axes.labelpad": 4.5,
        "axes.linewidth": 2.0,
        "xtick.labelsize": 16.0,
        "xtick.major.size": 4.0,
        "xtick.major.width": 1.25,
        "ytick.labelsize": 16.0,
        "ytick.major.size": 4.0,
        "ytick.major.width": 1.25,
        "figure.titlesize": 24.0,
        "figure.titleweight": "bold",
        "figure.labelsize": 20.0,
        "patch.linewidth": 2.0,
        "grid.linewidth": 1.25,
        "lines.linewidth": 2.0,
        "legend.fontsize": 16.0,
        "legend.title_fontsize": 16.0,
        "figure.constrained_layout.h_pad": 1 / 72,
        "figure.constrained_layout.w_pad": 1 / 72,
        "figure.dpi": 72.0,
        "savefig.dpi": 300.0,
    }
)

#: Registry. Adding a third mode is an entry here plus a row in MODE_META.
MODES = {
    "article": ARTICLE,
    "poster": POSTER,
}

#: Per-mode values that are not rcParams: the physical widths a size token
#: resolves to, and the vertical pitch one categorical row occupies. Poster
#: pitch is roughly double article, matching the type scale.
MODE_META = {
    "article": {
        "widths": {"full": 7.0, "half": 3.4, "third": 2.25},
        "row_pitch": 0.20,
        "min_height": 1.2,
    },
    "poster": {
        "widths": {"full": 14.0, "half": 6.8, "third": 4.5},
        "row_pitch": 0.42,
        "min_height": 2.5,
    },
}


def resolve(mode):
    """Return ``(rcparams, meta)`` for a mode name, raising on an unknown one."""
    if mode not in MODES:
        known = ", ".join(sorted(MODES))
        raise ValueError(f"unknown mode {mode!r}; expected one of: {known}")
    return MODES[mode], MODE_META[mode]
