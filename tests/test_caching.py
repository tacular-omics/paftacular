"""Parser caches: repeated substrings share one immutable component, and the caches stay bounded."""

import pytest

import paftacular as pft
from paftacular import parser
from paftacular.constants import MAX_CACHE_SIZE


@pytest.fixture(autouse=True)
def _empty_caches():
    parser._clear_caches()
    yield
    parser._clear_caches()


@pytest.mark.parametrize(
    ("text", "field"),
    [
        ("y5-H2O", "neutral_losses"),
        ("y5+2i13C", "isotopes"),
        ("y5[M+Na]", "adducts"),
    ],
)
def test_repeated_modifier_text_shares_components(text, field):
    first = pft.parse(text)
    second = pft.parse(text)
    assert first == second
    assert getattr(first, field) is getattr(second, field)


def test_repeated_ion_text_shares_ion_component():
    assert pft.parse("r[TMT127N]").ion_type is pft.parse("r[TMT127N]^2").ion_type


def test_repeated_mass_error_shares_component():
    assert pft.parse("y5/1.2ppm").mass_error is pft.parse("b3/1.2ppm").mass_error


def test_constructed_components_are_not_interned():
    # Constructors build fresh objects. Only the parser caches, by substring.
    assert pft.NeutralLoss(-1, base_formula="H2O") == pft.NeutralLoss(-1, base_formula="H2O")
    assert pft.NeutralLoss(-1, base_formula="H2O") is not pft.NeutralLoss(-1, base_formula="H2O")


def test_ion_cache_is_bounded():
    for index in range(MAX_CACHE_SIZE + 25):
        pft.parse(f"?{index}")
    assert len(parser._ION_CACHE) == MAX_CACHE_SIZE
    assert "?0" not in parser._ION_CACHE
    assert f"?{MAX_CACHE_SIZE + 24}" in parser._ION_CACHE


def test_modifier_cache_is_bounded():
    for index in range(MAX_CACHE_SIZE + 5):
        pft.parse(f"y2-{index + 1}.5")
    assert len(parser._LOSS_CACHE) == MAX_CACHE_SIZE


def test_clear_caches_empties_every_cache():
    pft.parse("y5-H2O+i[M+Na]/1.2ppm")
    parser._clear_caches()
    assert len(parser._ION_CACHE) == 0
    for cache in (parser._LOSS_CACHE, parser._ISOTOPE_CACHE, parser._ADDUCT_CACHE, parser._MASS_ERROR_CACHE):
        assert len(cache) == 0
