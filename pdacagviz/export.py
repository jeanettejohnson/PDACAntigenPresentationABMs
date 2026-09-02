"""Writing figures out.

One call replaces the ``os.path.join(ANALYSIS_DIR, ...)`` plus ``fig.savefig``
pair that each archived script assembled for itself, and fixes the thing they
mostly got wrong: only one of the eleven wrote a vector copy. PDF is the
format that cannot be recovered after the fact -- a raster figure cannot be
turned back into editable curves when a journal asks at submission -- so it is
written every time unless deliberately turned off.
"""

from pathlib import Path

import matplotlib as mpl
import pandas as pd
import seaborn as sns

from .settings import settings

__all__ = ["save", "stack_versions"]

#: Formats that carry embeddable metadata, and the key to put a note under.
#: PDF uses the standard document-info fields; PNG takes arbitrary text keys.
_METADATA_KEY = {"pdf": "Creator", "png": "Software", "svg": "Creator"}


def stack_versions():
    """One-line record of what drew the figure.

    Stamped into saved files so a panel in a manuscript draft can be traced
    back to the stack that produced it -- which matters precisely because the
    figure stack is now pinned and could later be re-pinned.
    """
    return (
        f"pdacagviz (matplotlib {mpl.__version__}, seaborn {sns.__version__}, "
        f"pandas {pd.__version__})"
    )


def save(fig, name, formats=None, directory=None, dpi=None, metadata=True, **kwargs):
    """Write ``fig`` under ``name`` in every configured format.

    Parameters
    ----------
    fig
        The figure. An Axes is accepted too, since chart functions return one
        and ``save(ax, ...)`` is what people will type.
    name
        Base filename, without extension. A name carrying one is accepted and
        the extension is stripped, so ``save(fig, "counts.pdf")`` does not
        produce ``counts.pdf.pdf``.
    formats
        Defaults to ``settings.formats``. Pass a single string for one format.
    directory
        Defaults to ``settings.figure_dir``. Created if absent.
    dpi
        Defaults to the current mode's save resolution.
    metadata
        Stamp :func:`stack_versions` into the file. Pass a dict to add fields
        of your own, or False to write nothing.

    Returns the list of paths written, in the order given.
    """
    figure = getattr(fig, "figure", fig)

    formats = settings.formats if formats is None else formats
    if isinstance(formats, str):
        formats = (formats,)
    formats = tuple(formats)
    if not formats:
        raise ValueError("save() needs at least one format")

    directory = Path(settings.figure_dir if directory is None else directory)
    directory.mkdir(parents=True, exist_ok=True)

    stem = Path(name).stem if Path(name).suffix.lstrip(".").lower() in formats else name
    dpi = settings.save_dpi if dpi is None else dpi

    written = []
    for fmt in formats:
        path = directory / f"{stem}.{fmt}"
        meta = None
        if metadata:
            meta = {}
            key = _METADATA_KEY.get(fmt)
            if key:
                meta[key] = stack_versions()
            if isinstance(metadata, dict):
                meta.update(metadata)
            meta = meta or None
        figure.savefig(path, format=fmt, dpi=dpi, metadata=meta, **kwargs)
        written.append(path)
    return written
