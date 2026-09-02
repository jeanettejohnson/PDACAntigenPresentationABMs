"""Palettes: the two authorities, the alias map, and strict lookup."""

import matplotlib.colors as mcolors
import pytest

import pdacagviz
from pdacagviz import palettes


@pytest.fixture(autouse=True)
def restore():
    before = (pdacagviz.settings.palette, pdacagviz.settings.strict_palette)
    yield
    pdacagviz.settings.palette, pdacagviz.settings.strict_palette = before


class TestTables:
    def test_both_authorities_are_present(self):
        assert len(palettes.ATLAS) == 46
        assert len(palettes.SIMULATION) == 21

    def test_every_colour_is_valid_to_matplotlib(self):
        for name, table in palettes.PALETTES.items():
            for key, colour in table.items():
                assert mcolors.is_color_like(colour), f"{name}[{key}] = {colour}"

    def test_simulation_keeps_the_css_names(self):
        # Byte-identical to kCellTypeColors in custom.cpp; the drift guard in
        # sanity_checks/ should assert this against the C++ source.
        assert palettes.SIMULATION["CAF"] == "yellow"
        assert palettes.SIMULATION["CD8_Tcell"] == "maroon"
        assert palettes.SIMULATION["apCAF"] == "grey"

    def test_atlas_keeps_the_hex_family_shading(self):
        assert palettes.ATLAS["CAF"] == "#f4a460"
        assert palettes.ATLAS["CD8_Tcell"] == "#8B0000"

    def test_the_two_authorities_genuinely_disagree(self):
        # If this ever passes trivially, one palette has been overwritten by
        # the other and the movie-matched figures have silently changed.
        shared = set(palettes.ATLAS) & set(palettes.SIMULATION)
        assert any(palettes.ATLAS[k] != palettes.SIMULATION[k] for k in shared)

    def test_antigen_class_is_the_plotting_variant(self):
        # Not the analysis/ variant, which reused tumour-lineage colours.
        assert palettes.ANTIGEN_CLASS["class I+II"] == "#1f78b4"
        assert "#00CED1" not in palettes.ANTIGEN_CLASS.values()

    def test_the_four_missing_agent_types_are_covered(self):
        # These are the entries plot_agent_counts_per_sample.py lacked, and
        # drew as silent grey.
        for name in (
            "epithelial_mesenchymal",
            "epithelial_mesenchymal_class1",
            "epithelial_mesenchymal_class2",
            "epithelial_mesenchymal_class1_class2",
        ):
            assert name in palettes.ATLAS


class TestNormalize:
    def test_simulation_names_map_to_atlas(self):
        assert palettes.normalize("epithelial_tumor") == "epithelial"
        assert palettes.normalize("B cell") == "B"
        assert palettes.normalize("macrophage") == "Macrophage"

    def test_atlas_names_pass_through(self):
        assert palettes.normalize("CD8_Tcell") == "CD8_Tcell"

    def test_every_alias_target_exists_in_atlas(self):
        for target in palettes._ALIASES.values():
            assert target in palettes.ATLAS


class TestLookup:
    def test_returns_colours_in_input_order(self):
        assert pdacagviz.colors_for(["CAF", "Treg"], palette="atlas") == [
            palettes.ATLAS["CAF"],
            palettes.ATLAS["Treg"],
        ]

    def test_as_dict(self):
        got = pdacagviz.colors_for(["CAF"], palette="atlas", as_dict=True)
        assert got == {"CAF": palettes.ATLAS["CAF"]}

    def test_simulation_names_resolve_against_atlas_via_aliases(self):
        assert pdacagviz.colors_for(["epithelial_tumor"], palette="atlas") == [
            palettes.ATLAS["epithelial"]
        ]

    def test_strict_lookup_raises_rather_than_greying(self):
        with pytest.raises(KeyError, match="no colour"):
            pdacagviz.colors_for(["not_a_cell_type"], palette="atlas", strict=True)

    def test_strict_error_names_every_missing_value_at_once(self):
        with pytest.raises(KeyError) as excinfo:
            pdacagviz.colors_for(["nope_a", "nope_b"], palette="atlas", strict=True)
        assert "nope_a" in str(excinfo.value) and "nope_b" in str(excinfo.value)

    def test_non_strict_falls_back_to_grey(self):
        assert pdacagviz.colors_for(["not_a_cell_type"], palette="atlas", strict=False) == [
            palettes.MISSING_COLOR
        ]

    def test_settings_supply_the_default_palette(self):
        pdacagviz.settings.palette = "simulation"
        assert pdacagviz.colors_for(["CAF"]) == ["yellow"]

    def test_an_explicit_mapping_is_accepted(self):
        assert pdacagviz.colors_for(["x"], palette={"x": "#123456"}) == ["#123456"]

    def test_unknown_palette_name_is_refused(self):
        with pytest.raises(ValueError, match="unknown palette"):
            pdacagviz.colors_for(["CAF"], palette="pastel")


class TestGlasbey:
    """The default categorical palette, vendored from colorcet."""

    def test_full_palette_is_present(self):
        assert len(palettes.GLASBEY_DARK) == 256

    def test_every_entry_is_a_valid_colour(self):
        assert all(mcolors.is_color_like(c) for c in palettes.GLASBEY_DARK)

    def test_no_duplicates(self):
        assert len(set(palettes.GLASBEY_DARK)) == len(palettes.GLASBEY_DARK)

    def test_categorical_takes_a_prefix(self):
        assert pdacagviz.categorical(5) == list(palettes.GLASBEY_DARK[:5])

    def test_categorical_cycles_past_the_end(self):
        got = pdacagviz.categorical(300)
        assert len(got) == 300
        assert got[256] == palettes.GLASBEY_DARK[0]

    def test_categorical_of_zero_is_empty(self):
        assert pdacagviz.categorical(0) == []

    def test_twenty_one_categories_stay_distinct(self):
        # matplotlib's ten-colour cycle repeats at the eleventh group; this data
        # has 21 patients, so a repeat would read as two patients being one.
        assert len(set(pdacagviz.categorical(21))) == 21


class TestPatientPalette:
    def test_covers_every_patient_in_the_run(self):
        assert len(palettes.PATIENT) == 21

    def test_colours_come_from_glasbey(self):
        assert set(palettes.PATIENT.values()) <= set(palettes.GLASBEY_DARK)

    def test_patients_are_distinguishable(self):
        assert len(set(palettes.PATIENT.values())) == len(palettes.PATIENT)

    def test_lookup_no_longer_raises_for_this_cohort(self):
        # The gap that motivated the rebuild: 8 entries against 21 patients, 7
        # overlapping, which is a partial match and so a strict-lookup failure.
        cohort = list(palettes.PATIENT)
        assert len(pdacagviz.colors_for(cohort, palette="patient", strict=True)) == 21
