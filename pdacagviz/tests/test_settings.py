"""Settings and mode switching.

The chart smoke test arrives with the charts. These cover the foundation the
charts will rest on: that a mode reaches rcParams, that scoping restores what
it borrowed, and that the compat layer absorbs version differences instead of
letting them reach a call site.
"""

import matplotlib as mpl
import pytest

import pdacagviz
from pdacagviz import _compat, modes
from pdacagviz._defaults import MPL_ALIASES, set_default


@pytest.fixture(autouse=True)
def restore_settings():
    """Leave global state as it was found, whatever a test does to it."""
    before_mode = pdacagviz.settings.mode
    before_rc = dict(mpl.rcParams)
    yield
    pdacagviz.settings.mode = before_mode
    mpl.rcParams.update(before_rc)


class TestModes:
    def test_switching_mode_reaches_rcparams(self):
        pdacagviz.settings.mode = "article"
        article_fontsize = mpl.rcParams["font.size"]
        pdacagviz.settings.mode = "poster"
        assert mpl.rcParams["font.size"] > article_fontsize

    def test_unknown_mode_is_rejected_before_it_is_stored(self):
        with pytest.raises(ValueError, match="unknown mode"):
            pdacagviz.settings.mode = "billboard"
        assert pdacagviz.settings.mode in modes.MODES

    def test_figsize_is_not_a_mode_concern(self):
        # Size follows from content, so a mode must not pin figure.figsize.
        for table in modes.MODES.values():
            assert "figure.figsize" not in table

    def test_constrained_layout_owns_the_margin(self):
        # savefig.bbox="tight" would re-crop after constrained layout resolved.
        for table in modes.MODES.values():
            assert table["figure.constrained_layout.use"] is True
            assert "savefig.bbox" not in table

    def test_vector_output_embeds_real_fonts(self):
        for table in modes.MODES.values():
            assert table["pdf.fonttype"] == 42
            assert table["ps.fonttype"] == 42

    def test_every_mode_has_matching_metadata(self):
        assert set(modes.MODES) == set(modes.MODE_META)
        for meta in modes.MODE_META.values():
            assert set(meta["widths"]) == {"full", "half", "third"}


class TestScoping:
    def test_using_restores_the_previous_mode(self):
        pdacagviz.settings.mode = "article"
        with pdacagviz.settings.using(mode="poster"):
            assert pdacagviz.settings.mode == "poster"
        assert pdacagviz.settings.mode == "article"

    def test_using_restores_rcparams_too(self):
        pdacagviz.settings.mode = "article"
        before = mpl.rcParams["font.size"]
        with pdacagviz.settings.using(mode="poster"):
            assert mpl.rcParams["font.size"] != before
        assert mpl.rcParams["font.size"] == before

    def test_using_restores_after_an_exception(self):
        pdacagviz.settings.mode = "article"
        with pytest.raises(RuntimeError):
            with pdacagviz.settings.using(mode="poster"):
                raise RuntimeError("boom")
        assert pdacagviz.settings.mode == "article"

    def test_unknown_setting_is_refused(self):
        with pytest.raises(AttributeError):
            with pdacagviz.settings.using(pallete="atlas"):  # misspelled
                pass


class TestConfigure:
    def test_mode_may_be_positional(self):
        pdacagviz.configure("poster")
        assert pdacagviz.settings.mode == "poster"

    def test_several_options_at_once(self):
        pdacagviz.configure(mode="poster", palette="simulation")
        assert pdacagviz.settings.palette == "simulation"

    def test_no_environment_variable_is_consulted(self, monkeypatch):
        monkeypatch.setenv("PDACAGVIZ_MODE", "poster")
        pdacagviz.configure(mode="article")
        assert pdacagviz.settings.mode == "article"


class TestCompat:
    def test_unknown_rcparams_are_dropped_not_raised(self):
        filtered = _compat.filter_rcparams(
            {"font.size": 7.0, "figure.thisWillNeverExist": 1}
        )
        assert filtered == {"font.size": 7.0}

    def test_mode_tables_survive_the_installed_matplotlib(self):
        for name, table in modes.MODES.items():
            filtered = _compat.filter_rcparams(table)
            mpl.rcParams.update(filtered)  # must not raise
            assert filtered, f"{name} filtered down to nothing"

    def test_seaborn_keywords_pass_through_on_current_seaborn(self):
        out, dropped = _compat.translate_seaborn({"density_norm": "width"})
        if _compat.SNS_VERSION >= _compat.Version("0.13"):
            assert out == {"density_norm": "width"} and not dropped
        else:
            assert out == {"scale": "width"}

    def test_unmapped_keywords_are_left_alone(self):
        out, dropped = _compat.translate_seaborn({"x": "a", "hue": "b"})
        assert out == {"x": "a", "hue": "b"} and not dropped

    def test_colormap_lookup_works_across_the_registry_boundary(self):
        assert _compat.get_colormap("viridis") is not None


class TestSetDefault:
    def test_caller_wins_over_default(self):
        params = {"lw": 2.0}
        set_default(MPL_ALIASES["linewidth"], 0.5, params)
        assert params == {"lw": 2.0}

    def test_default_applies_when_absent(self):
        params = {}
        set_default(MPL_ALIASES["linewidth"], 0.5, params)
        assert params == {"linewidth": 0.5}

    def test_nested_mappings_merge_rather_than_replace(self):
        params = {"boxprops": {"color": "red"}}
        set_default("boxprops", {"color": "black", "linewidth": 1.0}, params)
        assert params["boxprops"] == {"color": "red", "linewidth": 1.0}
