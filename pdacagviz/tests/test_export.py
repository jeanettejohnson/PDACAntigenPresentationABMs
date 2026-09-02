"""Saving figures: formats, resolution, and the version stamp."""

import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest

import pdacagviz


@pytest.fixture(autouse=True)
def restore():
    before = (pdacagviz.settings.mode, pdacagviz.settings.formats)
    rc = dict(mpl.rcParams)
    yield
    pdacagviz.settings.mode, pdacagviz.settings.formats = before
    mpl.rcParams.update(rc)
    plt.close("all")


def _figure():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    return fig, ax


class TestSave:
    def test_vector_is_written_every_time(self):
        # Only one of the eleven archived scripts emitted a PDF; it is the
        # output that cannot be recovered later.
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            written = pdacagviz.save(fig, "demo", directory=d)
            assert any(p.suffix == ".pdf" for p in written)
            assert all(p.exists() and p.stat().st_size > 0 for p in written)

    def test_default_formats_come_from_settings(self):
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            written = pdacagviz.save(fig, "demo", directory=d)
            assert [p.suffix.lstrip(".") for p in written] == list(
                pdacagviz.settings.formats
            )

    def test_an_axes_is_accepted_not_only_a_figure(self):
        _, ax = _figure()
        with tempfile.TemporaryDirectory() as d:
            assert pdacagviz.save(ax, "demo", directory=d, formats="png")

    def test_a_name_carrying_an_extension_is_not_doubled(self):
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            written = pdacagviz.save(fig, "demo.png", directory=d, formats="png")
            assert written[0].name == "demo.png"

    def test_a_dot_in_the_name_is_preserved(self):
        # "counts.v2" must not lose its suffix to extension-stripping.
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            written = pdacagviz.save(fig, "counts.v2", directory=d, formats="png")
            assert written[0].name == "counts.v2.png"

    def test_missing_directory_is_created(self):
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            target = Path(d) / "deep" / "nested"
            written = pdacagviz.save(fig, "demo", directory=target, formats="png")
            assert written[0].exists()

    def test_single_format_as_a_string(self):
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            assert len(pdacagviz.save(fig, "demo", directory=d, formats="png")) == 1

    def test_no_formats_is_an_error(self):
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            with pytest.raises(ValueError, match="at least one format"):
                pdacagviz.save(fig, "demo", directory=d, formats=())


class TestVersionStamp:
    def test_stamp_names_the_drawing_stack(self):
        stamp = pdacagviz.stack_versions()
        assert "matplotlib" in stamp and "seaborn" in stamp and "pandas" in stamp

    def test_stamp_reaches_the_pdf(self):
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            path = pdacagviz.save(fig, "demo", directory=d, formats="pdf")[0]
            assert b"pdacagviz" in path.read_bytes()

    def test_metadata_can_be_turned_off(self):
        fig, _ = _figure()
        with tempfile.TemporaryDirectory() as d:
            path = pdacagviz.save(fig, "demo", directory=d, formats="pdf", metadata=False)[0]
            assert path.exists()
