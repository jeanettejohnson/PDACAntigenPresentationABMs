"""Axes helpers: the lines the archived scripts each wrote by hand."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest

import pdacagviz
from pdacagviz.axes import PVAL_THRESHOLDS


@pytest.fixture(autouse=True)
def restore():
    rc = dict(mpl.rcParams)
    yield
    mpl.rcParams.update(rc)
    plt.close("all")


class TestDespine:
    def test_top_and_right_go_by_default(self):
        _, ax = plt.subplots()
        pdacagviz.despine(ax)
        assert not ax.spines["top"].get_visible()
        assert not ax.spines["right"].get_visible()
        assert ax.spines["left"].get_visible()

    def test_nothing_removed_when_all_false(self):
        _, ax = plt.subplots()
        pdacagviz.despine(ax, top=False, right=False)
        assert all(ax.spines[s].get_visible() for s in ("top", "right", "left", "bottom"))


class TestThousands:
    def test_separator_is_inserted(self):
        _, ax = plt.subplots()
        ax.set_xlim(0, 20000)
        pdacagviz.thousands(ax.xaxis)
        assert ax.xaxis.get_major_formatter()(12345, 0) == "12,345"

    def test_decimals_are_honoured(self):
        _, ax = plt.subplots()
        pdacagviz.thousands(ax.xaxis, decimals=1)
        assert ax.xaxis.get_major_formatter()(12345.67, 0) == "12,345.7"


class TestLegendOutside:
    def test_legend_is_placed_beyond_the_axes(self):
        _, ax = plt.subplots()
        ax.plot([0, 1], [0, 1], label="a")
        legend = pdacagviz.legend_outside(ax)
        assert legend is not None
        assert legend.get_bbox_to_anchor() is not None
        assert not legend.get_frame_on()

    def test_no_handles_means_no_legend(self):
        _, ax = plt.subplots()
        assert pdacagviz.legend_outside(ax) is None

    def test_ncols_works_across_the_3_6_rename(self):
        _, ax = plt.subplots()
        for i in range(4):
            ax.plot([0, 1], [i, i], label=str(i))
        assert pdacagviz.legend_outside(ax, ncols=2) is not None


class TestBarLabels:
    def test_one_label_per_bar(self):
        _, ax = plt.subplots()
        bars = ax.barh(["a", "b"], [10, 20])
        before = len(ax.texts)
        pdacagviz.bar_labels(ax, bars, [10, 20])
        assert len(ax.texts) - before == 2

    def test_values_are_formatted(self):
        _, ax = plt.subplots()
        bars = ax.barh(["a"], [12345])
        pdacagviz.bar_labels(ax, bars, [12345])
        assert ax.texts[-1].get_text() == "12,345"

    def test_custom_format(self):
        _, ax = plt.subplots()
        bars = ax.barh(["a"], [0.42])
        pdacagviz.bar_labels(ax, bars, [0.42], fmt="{:.0%}")
        assert ax.texts[-1].get_text() == "42%"

    def test_headroom_leaves_space_past_the_longest_bar(self):
        _, ax = plt.subplots()
        bars = ax.barh(["a"], [100])
        pdacagviz.bar_labels(ax, bars, [100])
        assert ax.get_xlim()[1] > 100

    def test_length_mismatch_is_an_error(self):
        _, ax = plt.subplots()
        bars = ax.barh(["a", "b"], [1, 2])
        with pytest.raises(ValueError, match="values"):
            pdacagviz.bar_labels(ax, bars, [1])

    def test_empty_is_a_no_op(self):
        _, ax = plt.subplots()
        pdacagviz.bar_labels(ax, [], [])
        assert not ax.texts


class TestMarkerSize:
    def test_size_falls_as_points_are_added(self):
        assert pdacagviz.marker_size(100) > pdacagviz.marker_size(10000)

    def test_size_rises_with_drawing_area(self):
        small = pdacagviz.marker_size(1000, figsize=(2, 2))
        large = pdacagviz.marker_size(1000, figsize=(8, 8))
        assert large > small

    def test_a_floor_keeps_huge_counts_visible(self):
        assert pdacagviz.marker_size(10**9, figsize=(2, 2)) > 0

    def test_scale_multiplies(self):
        assert pdacagviz.marker_size(100, scale=2.0) == pytest.approx(
            pdacagviz.marker_size(100) * 2.0
        )

    def test_zero_points_is_an_error(self):
        with pytest.raises(ValueError, match="at least one"):
            pdacagviz.marker_size(0)


class TestPvalStars:
    def test_marks_accumulate_with_significance(self):
        assert len(pdacagviz.pval_stars(0.04)) == 1
        assert len(pdacagviz.pval_stars(0.004)) == 2
        assert len(pdacagviz.pval_stars(0.0004)) == 3

    def test_thresholds_are_inclusive(self):
        assert pdacagviz.pval_stars(PVAL_THRESHOLDS[-1]) != "ns"

    def test_non_significant_is_spelled_out(self):
        # "ns" rather than "" so a reader can tell it from a missing value.
        assert pdacagviz.pval_stars(0.5) == "ns"

    def test_scientific_notation_on_request(self):
        assert pdacagviz.pval_stars(0.00012, as_stars=False) == "1.20e-04"
