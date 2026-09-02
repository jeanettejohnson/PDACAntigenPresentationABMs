"""Module state, held in one object rather than in the environment.

Nothing here reads a shell variable. A figure's appearance is a property of
the code that drew it, not of the terminal it was launched from -- which is
what lets a figure be reproduced six months later from the script alone.

Three ways in, all reaching the same object::

    pdacagviz.settings.mode = "poster"            # applies immediately
    pdacagviz.configure(mode="poster")            # several at once
    with pdacagviz.settings.using(mode="poster"): # scoped, restores on exit
        ...

The split against rcParams is deliberate. rcParams validates its keys and
accepts only matplotlib's own, so ``palette``, ``figure_dir``, ``formats`` and
``strict_palette`` have no home there; smuggling them in under a custom prefix
is the kind of thing that breaks on a matplotlib upgrade. Settings owns
project-level state, rcParams owns matplotlib-level state, and assigning
``mode`` writes through to rcParams as a side effect.
"""

import contextlib
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib as mpl

# tomllib is stdlib from python 3.11. The floor here is 3.10, so fall back to
# tomli, and if neither is present simply skip the project-file layer -- it is
# a convenience over the built-in defaults, never a requirement.
try:
    import tomllib as _toml
except ModuleNotFoundError:  # python 3.10
    try:
        import tomli as _toml
    except ModuleNotFoundError:
        _toml = None

from . import modes as _modes
from ._compat import filter_rcparams

__all__ = ["Settings", "settings", "configure"]

_SENTINEL = object()


def _find_repo_root(start=None):
    """Walk up from this file looking for the repository root.

    Recognised by a .git entry or a pyproject.toml. Falls back to the parent
    of the package directory, which is correct for the layout as shipped.
    """
    here = Path(start or __file__).resolve()
    for candidate in here.parents:
        if (candidate / ".git").exists() or (candidate / "pyproject.toml").is_file():
            return candidate
    return here.parent.parent


def _load_project_defaults(root):
    """Read ``[tool.pdacagviz]`` out of the repository's pyproject.toml.

    A checked-in file is what makes a project-wide choice diffable and
    greppable, which is the property an environment variable lacks. Absent or
    unreadable, this returns nothing and the built-in defaults stand.
    """
    if _toml is None:
        return {}
    path = root / "pyproject.toml"
    if not path.is_file():
        return {}
    try:
        with path.open("rb") as fh:
            return _toml.load(fh).get("tool", {}).get("pdacagviz", {})
    except (OSError, _toml.TOMLDecodeError):
        return {}


@dataclass
class Settings:
    """Project-level plotting state.

    Attributes
    ----------
    mode
        ``"article"`` or ``"poster"``. Assigning applies that mode's rcParams
        immediately.
    palette
        Default palette name for charts that colour by cell type.
    figure_dir
        Where :func:`pdacagviz.save` writes. Resolved from the repository root.
    formats
        Extensions written on every save. PDF is included by default because
        it is the output that cannot be regenerated after the fact when a
        journal asks at submission.
    dpi
        ``None`` means use the mode's own value.
    display_dpi
        Screen resolution for inline notebook rendering. Kept separate from
        ``dpi`` on purpose: article versus poster is about physical size, not
        about screen versus file, so a legible preview must not change what
        lands in the saved file.
    strict_palette
        Raise on a cell type with no colour rather than returning grey. The
        silent-grey fallback is how four missing entries went unnoticed in the
        scripts this module replaces.
    """

    mode: str = "article"
    palette: str = "atlas"
    figure_dir: Path = field(default_factory=lambda: _find_repo_root() / "figures")
    formats: tuple = ("png", "pdf")
    dpi: float | None = None
    display_dpi: float = 96.0
    strict_palette: bool = True

    def __post_init__(self):
        for key, value in _load_project_defaults(_find_repo_root()).items():
            if not hasattr(self, key):
                raise ValueError(
                    f"[tool.pdacagviz] sets unknown option {key!r} in pyproject.toml"
                )
            object.__setattr__(self, key, value)
        self.apply()

    def __setattr__(self, name, value):
        if name == "mode":
            _modes.resolve(value)  # reject an unknown mode before storing it
        known = hasattr(type(self), name) or name in self.__dataclass_fields__
        if not known:
            raise AttributeError(f"unknown pdacagviz setting {name!r}")
        object.__setattr__(self, name, value)
        if name == "mode":
            self.apply()

    # -- applying -----------------------------------------------------------

    def apply(self):
        """Push the current mode's rcParams into matplotlib.

        Called automatically when ``mode`` is assigned. Useful directly after
        something else has reset rcParams underneath -- ``plt.style.use`` and
        ``sns.set_theme`` both do.
        """
        rcparams, _ = _modes.resolve(self.mode)
        mpl.rcParams.update(filter_rcparams(rcparams))
        if self.dpi is not None:
            mpl.rcParams["figure.dpi"] = self.dpi
        return self

    @property
    def meta(self):
        """Non-rcParam values for the current mode: widths, row pitch."""
        return _modes.resolve(self.mode)[1]

    @property
    def save_dpi(self):
        """Resolution :func:`pdacagviz.save` writes at."""
        if self.dpi is not None:
            return self.dpi
        return mpl.rcParams["savefig.dpi"]

    # -- scoping ------------------------------------------------------------

    @contextlib.contextmanager
    def using(self, **overrides):
        """Apply settings for the duration of a block, then restore them.

        This matters more in a notebook than in a script. A notebook keeps one
        long-lived interpreter, so a bare ``settings.mode = "poster"`` in one
        cell silently changes every figure in the cells after it. Scoping is
        the difference between one poster figure and a notebook that quietly
        stops matching the manuscript::

            with pdacagviz.settings.using(mode="poster"):
                ax = pdacagviz.bar(df, x="agent", y="count")
                pdacagviz.save(ax.figure, "agent_counts_poster")
        """
        for key in overrides:
            if key not in self.__dataclass_fields__:
                raise AttributeError(f"unknown pdacagviz setting {key!r}")

        previous = {key: getattr(self, key) for key in overrides}
        saved_rcparams = dict(mpl.rcParams)
        try:
            for key, value in overrides.items():
                setattr(self, key, value)
            yield self
        finally:
            for key, value in previous.items():
                object.__setattr__(self, key, value)
            mpl.rcParams.update(saved_rcparams)


#: The module-wide instance. Import and mutate this rather than making another.
settings = Settings()


def configure(mode=_SENTINEL, **overrides):
    """Set several options at once and apply the result.

    ``mode`` may be given positionally, since it is the one people reach for::

        pdacagviz.configure("poster")
        pdacagviz.configure(mode="poster", palette="simulation")
    """
    if mode is not _SENTINEL:
        overrides["mode"] = mode
    for key, value in overrides.items():
        setattr(settings, key, value)
    return settings.apply()
