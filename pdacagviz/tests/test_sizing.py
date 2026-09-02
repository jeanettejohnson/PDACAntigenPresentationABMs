"""Sizing: width tokens, content-driven height, panel grids."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest

import pdacagviz


@pytest.fixture(autouse=True)
def restore():
    before, rc = pdacagviz.settings.mode, dict(mpl.rcParams)
    yield
    pdacagviz.settings.mode = before
    mpl.rcParams.update(rc)
    plt.close("all")


class TestWidth:
    def test_tokens_scale_with_mode(self):
        pdacagviz.settings.mode = "article"
        article = pdacagviz.width_of("full")
        pdacagviz.settings.mode = "poster"
        assert pdacagviz.width_of("full") > article

    def test_tokens_are_ordered(self):
        assert pdacagviz.width_of("full") > pdacagviz.width_of("half") > pdacagviz.width_of("third")

    def test_a_number_passes_through(self):
        assert pdacagviz.width_of(4.25) == 4.25

    def test_unknown_token_is_refused(self):
        with pytest.raises(ValueError, match="unknown width"):
            pdacagviz.width_of("gigantic")


class TestFigsize:
    def test_height_grows_with_rows(self):
        assert pdacagviz.figsize("full", rows=40)[1] > pdacagviz.figsize("full", rows=5)[1]

    def test_row_pitch_scales_with_mode(self):
        pdacagviz.settings.mode = "article"
        article = pdacagviz.figsize("full", rows=20)[1]
        pdacagviz.settings.mode = "poster"
        assert pdacagviz.figsize("full", rows=20)[1] > article

    def test_ratio_applies_to_width(self):
        w, h = pdacagviz.figsize("full", ratio=0.5)
        assert h == pytest.approx(w * 0.5)

    def test_rows_and_ratio_together_is_an_error(self):
        # Two different intents; choosing for the caller would hide a mistake.
        with pytest.raises(ValueError, match="not both"):
            pdacagviz.figsize("full", rows=10, ratio=0.5)

    def test_negative_rows_is_an_error(self):
        with pytest.raises(ValueError, match="negative"):
            pdacagviz.figsize("full", rows=-1)

    def test_zero_rows_still_has_height(self):
        assert pdacagviz.figsize("full", rows=0)[1] > 0


class TestGrid:
    def test_axes_are_always_two_dimensional(self):
        _, axes = pdacagviz.grid(1, 1)
        assert axes.ndim == 2

    def test_figure_grows_with_panel_count(self):
        one = pdacagviz.grid(1, 1, width="half")[0].get_size_inches()
        six = pdacagviz.grid(2, 3, width="half")[0].get_size_inches()
        assert six[0] > one[0] and six[1] > one[1]

    def test_panels_keep_their_size_as_columns_are_added(self):
        # The figure grows rather than dividing a fixed canvas, which is what
        # keeps type legible when a panel is added.
        w1 = pdacagviz.grid(1, 1, width="half")[0].get_size_inches()[0]
        w3 = pdacagviz.grid(1, 3, width="half")[0].get_size_inches()[0]
        assert w3 == pytest.approx(w1 * 3, rel=0.1)

    def test_panel_size_overrides_tokens(self):
        fig, _ = pdacagviz.grid(1, 1, panel_size=(3.0, 2.0))
        assert tuple(fig.get_size_inches()) == pytest.approx((3.0, 2.0))

    def test_empty_grid_is_refused(self):
        with pytest.raises(ValueError, match="at least one"):
            pdacagviz.grid(0, 3)
