"""Applying defaults without overriding what the caller asked for.

Ported from ``scanpy_extensions._common``, minus its scanpy import. The
original module carried a ``RandomState`` alias that pulled in scanpy purely
to pick between two type names; nothing here needs it.
"""

from collections.abc import Iterable, Mapping
from typing import Any

__all__ = ["isiterable", "set_default", "MPL_ALIASES"]


#: Matplotlib accepts several spellings for the same keyword. Treating them as
#: one key is what keeps ``set_default`` from adding ``linewidth`` alongside a
#: caller's ``lw`` and letting the two fight inside the plotting call.
MPL_ALIASES: dict[str, tuple[str, ...]] = {
    "linewidth": ("linewidth", "linewidths", "lw"),
    "linecolor": ("linecolor", "linecolors", "lc"),
    "edgecolor": ("edgecolor", "edgecolors", "ec"),
    "facecolor": ("facecolor", "facecolors", "fc"),
    "alpha": ("alpha", "a"),
    "size": ("size", "sizes", "s"),
    "markersize": ("markersize", "ms"),
    "markeredgecolor": ("markeredgecolor", "mec"),
    "markeredgewidth": ("markeredgewidth", "mew"),
    "markerfacecolor": ("markerfacecolor", "mfc"),
}


def isiterable(x: Any) -> bool:
    """True for iterables that are not strings.

    Strings are iterable but are almost always meant as a single value here --
    one column name, one colour -- so they are excluded deliberately.
    """
    return not isinstance(x, str) and isinstance(x, Iterable)


def set_default(key: Any, value: Any, config: dict[str, Any]) -> None:
    """Set ``key`` in ``config`` only if the caller has not already set it.

    ``key`` may be a single name or several names for the same option, in
    which case any one of them counts as already set and the first is used
    when writing. Pass a tuple from :data:`MPL_ALIASES` to cover matplotlib's
    spellings.

    When both the existing entry and ``value`` are mappings, ``value`` fills in
    absent sub-keys rather than replacing the mapping. That matters for nested
    style dicts such as ``boxprops``, where a caller who sets one property
    should not lose the defaults for the rest.

    Mutates ``config`` in place and returns nothing.
    """
    keys = tuple(key) if isiterable(key) else (key,)
    if not keys:
        raise ValueError("set_default() needs at least one key name")

    for name in keys:
        if name in config:
            existing = config[name]
            if isinstance(value, Mapping) and isinstance(existing, Mapping):
                for sub_key, sub_value in value.items():
                    existing.setdefault(sub_key, sub_value)
            return

    config[keys[0]] = value
