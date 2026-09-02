"""Plot functions, one module per data question.

Each function answers a question about the data and takes the mark as
``kind=``, rather than splitting into one function per mark. That is the
consolidation carried over from scanpy_extensions: when a reviewer asks to see
a distribution as a box instead of a violin, the edit is one argument, and the
surrounding call -- data, columns, palette, sizing -- does not move.

Every function shares one signature::

    f(data, x=, y=, hue=, kind=, panel=None, palette=None, ax=None,
      plot_kwargs={}, **fig_kwargs)

``plot_kwargs`` reaches the underlying seaborn or matplotlib call untouched;
``**fig_kwargs`` configures the figure. Learning one function teaches the rest.
"""

from ._bar import bar
from ._composition import composition, percentages
from ._distribution import distribution
from ._heatmap import heatmap
from ._relationship import relationship
from ._timecourse import timecourse

__all__ = [
    "bar",
    "composition",
    "percentages",
    "distribution",
    "relationship",
    "timecourse",
    "heatmap",
]
