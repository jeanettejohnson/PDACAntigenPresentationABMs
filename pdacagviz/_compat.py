"""Version differences, confined to one module.

The floors are deliberately low -- matplotlib 3.4, seaborn 0.11.1 -- so that
pdacagviz runs anywhere python 3.10 does. Supporting that span costs three
mechanisms, and they all live here rather than at the call sites:

1. :func:`filter_rcparams` drops rcParam keys the installed matplotlib does
   not know, so the mode tables in :mod:`pdacagviz.modes` can name keys from any
   version. This works in both directions: old matplotlib missing new keys, and
   future matplotlib having removed old ones.
2. :func:`translate_seaborn` rewrites modern seaborn keywords down to their
   older spellings, so chart code only ever writes the current names.
3. :func:`get_colormap` looks a colormap up by name across the 3.5 registry
   boundary.

Anything added here should cover a class of difference, not a single call.
"""

import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from packaging.version import Version

__all__ = [
    "MPL_VERSION",
    "SNS_VERSION",
    "filter_rcparams",
    "translate_seaborn",
    "get_colormap",
]

MPL_VERSION = Version(mpl.__version__)
SNS_VERSION = Version(sns.__version__)

_SNS_012 = SNS_VERSION >= Version("0.12")
_SNS_013 = SNS_VERSION >= Version("0.13")


def filter_rcparams(params, warn=False):
    """Return ``params`` without keys the installed matplotlib rejects.

    Assigning an unknown key to ``mpl.rcParams`` raises, so a mode table that
    names ``figure.labelsize`` (matplotlib 3.6) would break every older
    install. Filtering first is what lets one table serve every supported
    version.

    Set ``warn=True`` to hear about what was dropped; the default is silence,
    because on an older matplotlib the drops are expected rather than notable.
    """
    known, unknown = {}, []
    for key, value in params.items():
        if key in mpl.rcParams:
            known[key] = value
        else:
            unknown.append(key)

    if unknown and warn:
        warnings.warn(
            f"matplotlib {MPL_VERSION} does not know these rcParams, "
            f"so they were skipped: {', '.join(sorted(unknown))}",
            RuntimeWarning,
            stacklevel=2,
        )
    return known


#: Modern seaborn keyword -> (older keyword, the version that introduced the
#: modern one). ``None`` as the older keyword means drop it entirely, because
#: the option did not exist before and there is nothing to translate to.
_SEABORN_RENAMES = {
    "errorbar": ("ci", Version("0.12")),
    "bw_method": ("bw", Version("0.13")),
    "density_norm": ("scale", Version("0.13")),
    "native_scale": (None, Version("0.13")),
    "legend": (None, Version("0.13")),
}


def translate_seaborn(kwargs):
    """Rewrite modern seaborn keywords for the installed seaborn.

    Chart code always writes today's names. On an older seaborn this maps them
    back: ``errorbar`` to ``ci`` below 0.12, ``bw_method`` to ``bw`` and
    ``density_norm`` to ``scale`` below 0.13. Options with no older equivalent
    are dropped.

    Returns ``(translated_kwargs, dropped_names)``. ``dropped_names`` is what
    the caller has to handle another way -- most importantly ``legend``, which
    below 0.13 means removing the legend from the axes after plotting.

    One caveat this cannot paper over: ``errorbar`` and ``ci`` agree on
    ``None`` but not on their richer values, so anything beyond ``None`` needs
    checking against the floor rather than trusting this map.
    """
    out, dropped = {}, []
    for key, value in kwargs.items():
        if key not in _SEABORN_RENAMES:
            out[key] = value
            continue

        old_key, introduced_in = _SEABORN_RENAMES[key]
        if SNS_VERSION >= introduced_in:
            out[key] = value
        elif old_key is None:
            dropped.append(key)
        else:
            out[old_key] = value
    return out, dropped


def get_colormap(name):
    """Look up a colormap by name across matplotlib versions.

    ``mpl.colormaps`` is the modern registry and arrived in 3.5. Below that,
    ``plt.get_cmap`` is the way in -- note it is ``matplotlib.cm.get_cmap``
    that 3.9 removed, not the pyplot function, so this stays valid at both
    ends of the supported range.
    """
    try:
        return mpl.colormaps[name]
    except (AttributeError, KeyError):
        return plt.get_cmap(name)
