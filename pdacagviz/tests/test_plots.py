"""Every plot function, every kind, both modes.

This is the smoke test the version safeguards exist for: it renders each
function once with deprecation warnings escalated to errors, so the next
matplotlib or seaborn rename arrives as a red test rather than as a figure
that quietly changed.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import pdacagviz


@pytest.fixture(autouse=True)
def restore():
    before, rc = pdacagviz.settings.mode, dict(mpl.rcParams)
    yield
    pdacagviz.settings.mode = before
    mpl.rcParams.update(rc)
    plt.close("all")


def observations(n=240, seed=0):
    """Long-form frame shaped like the ABM output: cells with a type and a size."""
    rng = np.random.default_rng(seed)
    types = ["CD8_Tcell", "CD4_Tcell", "Treg", "CAF", "PDAC", "Macrophage"]
    return pd.DataFrame(
        {
            "cell_type": rng.choice(types, n),
            "condition": rng.choice(["control", "treated"], n),
            "sample": rng.choice(["HT056", "HT060", "HT064"], n),
            "volume": rng.gamma(4.0, 120.0, n),
            "signal": rng.normal(0.0, 1.0, n),
            "time": rng.integers(0, 12, n),
        }
    )


def counts():
    """One precomputed value per category, as bar() expects."""
    return pd.DataFrame(
        {
            "agent": ["CD8_Tcell", "CD4_Tcell", "Treg", "CAF", "PDAC"],
            "count": [48210, 31005, 9877, 60432, 91250],
        }
    )


class TestBar:
    def test_draws_one_bar_per_category(self):
        df = counts()
        ax = pdacagviz.bar(df, x="count", y="agent")
        assert len(ax.patches) == len(df)

    def test_values_are_labelled(self):
        ax = pdacagviz.bar(counts(), x="count", y="agent")
        assert any("," in t.get_text() for t in ax.texts)

    def test_vertical_orientation(self):
        ax = pdacagviz.bar(counts(), x="agent", y="count", orient="v")
        assert len(ax.patches) == 5

    def test_a_non_numeric_value_column_is_refused(self):
        with pytest.raises(TypeError, match="not numeric"):
            pdacagviz.bar(counts(), x="agent", y="agent")

    def test_missing_column_names_what_is_available(self):
        with pytest.raises(KeyError, match="Available"):
            pdacagviz.bar(counts(), x="nope", y="agent")


class TestComposition:
    @pytest.mark.parametrize("kind", ["stacked", "grouped"])
    def test_kinds_draw(self, kind):
        ax = pdacagviz.composition(observations(), x="condition", hue="cell_type", kind=kind)
        assert ax.patches

    def test_percentages_sum_to_a_hundred(self):
        table = pdacagviz.percentages(observations(), "condition", "cell_type")
        assert np.allclose(table.sum(axis=1), 100.0)

    def test_counts_when_not_normalised(self):
        df = observations()
        table = pdacagviz.percentages(df, "condition", "cell_type", norm=False)
        assert table.to_numpy().sum() == len(df)

    def test_diverging_needs_exactly_two_groups(self):
        with pytest.raises(ValueError, match="exactly two"):
            pdacagviz.composition(
                observations(), x="sample", hue="cell_type", kind="diverging"
            )

    def test_diverging_mirrors_about_zero(self):
        ax = pdacagviz.composition(
            observations(), x="condition", hue="cell_type", kind="diverging"
        )
        widths = [p.get_width() for p in ax.patches]
        assert min(widths) < 0 < max(widths)

    def test_unknown_kind_is_refused(self):
        with pytest.raises(ValueError, match="unknown kind"):
            pdacagviz.composition(observations(), x="condition", hue="cell_type", kind="pie")


class TestDistribution:
    @pytest.mark.parametrize("kind", ["violin", "box", "bar"])
    def test_kinds_draw(self, kind):
        ax = pdacagviz.distribution(observations(), x="cell_type", y="volume", kind=kind)
        assert ax.collections or ax.patches or ax.lines

    def test_hue_splits_the_categories(self):
        ax = pdacagviz.distribution(
            observations(), x="cell_type", y="volume", hue="condition", kind="box"
        )
        assert ax is not None

    def test_unknown_kind_is_refused(self):
        with pytest.raises(ValueError, match="unknown kind"):
            pdacagviz.distribution(observations(), x="cell_type", y="volume", kind="swarm")


class TestRelationship:
    @pytest.mark.parametrize("kind", ["scatter", "hist", "kde"])
    def test_kinds_draw(self, kind):
        ax = pdacagviz.relationship(observations(), x="volume", y="signal", kind=kind)
        assert ax.collections or ax.patches or ax.images

    def test_hue_is_coloured_from_the_palette(self):
        ax = pdacagviz.relationship(
            observations(), x="volume", y="signal", hue="cell_type", kind="scatter"
        )
        assert ax.collections

    def test_both_axes_are_required(self):
        with pytest.raises(ValueError, match="both x= and y="):
            pdacagviz.relationship(observations(), x="volume")


class TestTimecourse:
    def test_draws_a_line(self):
        ax = pdacagviz.timecourse(observations(), x="time", y="volume")
        assert ax.lines

    def test_one_line_per_hue_level(self):
        df = observations()
        ax = pdacagviz.timecourse(df, x="time", y="volume", hue="cell_type")
        assert len(ax.lines) >= df["cell_type"].nunique()

    def test_band_can_be_turned_off(self):
        ax = pdacagviz.timecourse(observations(), x="time", y="volume", band=None)
        assert ax.lines


class TestHeatmap:
    @pytest.mark.parametrize("kind", ["matrix", "dot"])
    def test_kinds_draw(self, kind):
        ax = pdacagviz.heatmap(
            observations(), x="condition", y="cell_type", value="volume", kind=kind
        )
        assert ax.collections

    def test_ticks_label_every_category(self):
        df = observations()
        ax = pdacagviz.heatmap(df, x="condition", y="cell_type", value="volume")
        assert len(ax.get_yticklabels()) == df["cell_type"].nunique()

    def test_vcenter_gives_a_diverging_scale(self):
        ax = pdacagviz.heatmap(
            observations(), x="condition", y="cell_type", value="signal", vcenter=0.0
        )
        assert ax.collections


class TestSharedSignature:
    def test_ax_is_honoured(self):
        _, ax = plt.subplots()
        assert pdacagviz.bar(counts(), x="count", y="agent", ax=ax) is ax

    def test_panel_facets_into_a_grid(self):
        result = pdacagviz.distribution(
            observations(), x="cell_type", y="volume", panel="condition"
        )
        assert result.ndim == 2

    def test_ax_and_panel_together_is_an_error(self):
        _, ax = plt.subplots()
        with pytest.raises(ValueError, match="not both"):
            pdacagviz.distribution(
                observations(), x="cell_type", y="volume", panel="condition", ax=ax
            )

    def test_a_misspelled_figure_option_is_reported(self):
        # Silently ignoring it would surface only after the figure was printed.
        with pytest.raises(TypeError, match="unknown option"):
            pdacagviz.bar(counts(), x="count", y="agent", despline=True)

    def test_plot_kwargs_reach_the_underlying_call(self):
        ax = pdacagviz.bar(counts(), x="count", y="agent", plot_kwargs={"alpha": 0.5})
        assert ax.patches[0].get_alpha() == 0.5

    def test_mode_changes_the_figure_size(self):
        with pdacagviz.settings.using(mode="article"):
            small = pdacagviz.bar(counts(), x="count", y="agent").figure.get_size_inches()[0]
        with pdacagviz.settings.using(mode="poster"):
            large = pdacagviz.bar(counts(), x="count", y="agent").figure.get_size_inches()[0]
        assert large > small


class TestPaletteFallback:
    """Strict lookup must protect cell-type columns without blocking others."""

    def test_a_non_cell_type_hue_gets_the_default_cycle(self):
        # condition/sample/treatment are not cell types; requiring strict=False
        # for them would disarm the check where it actually matters.
        ax = pdacagviz.distribution(
            observations(), x="cell_type", y="volume", hue="condition", kind="box"
        )
        assert ax is not None

    def test_a_partial_cell_type_column_still_raises(self):
        df = observations()
        df.loc[df.index[:20], "cell_type"] = "not_a_real_cell_type"
        with pytest.raises(KeyError, match="no colour"):
            pdacagviz.bar(
                df.groupby("cell_type", as_index=False)["volume"].sum(),
                x="volume",
                y="cell_type",
            )

    def test_cell_type_colours_still_come_from_the_atlas(self):
        ax = pdacagviz.bar(counts(), x="count", y="agent")
        drawn = {p.get_facecolor() for p in ax.patches}
        assert len(drawn) == len(counts())
