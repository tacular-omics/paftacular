from collections import Counter

import pytest
from tacular import ELEMENT_LOOKUP

from paftacular.comps import (
    composition_to_formula_string,
    composition_to_proforma_formula_string,
)


def test_hill_order_basic():
    """Hill order: C, H, then alphabetical"""
    comp = Counter(
        {
            ELEMENT_LOOKUP["O"]: 6,
            ELEMENT_LOOKUP["H"]: 12,
            ELEMENT_LOOKUP["C"]: 6,
        }
    )
    assert composition_to_formula_string(comp) == "C6H12O6"


def test_hill_order_basic2():
    """Hill order: C, H, then alphabetical"""
    comp = Counter(
        {
            ELEMENT_LOOKUP["O"]: 6,
            ELEMENT_LOOKUP["H"]: 12,
            ELEMENT_LOOKUP["13C"]: 6,
        }
    )
    assert composition_to_formula_string(comp) == "[13C6]H12O6"


def test_hill_order_basic3():
    """Hill order: C, H, then alphabetical"""
    comp = Counter(
        {
            ELEMENT_LOOKUP["O"]: 6,
            ELEMENT_LOOKUP["H"]: 12,
            ELEMENT_LOOKUP["C"]: 6,
            ELEMENT_LOOKUP["12C"]: 6,
            ELEMENT_LOOKUP["13C"]: 6,
        }
    )
    assert composition_to_formula_string(comp) == "C6[12C6][13C6]H12O6"


def test_hill_order_no_carbon():
    """No carbon: H first, then alphabetical"""
    comp = Counter({ELEMENT_LOOKUP["O"]: 1, ELEMENT_LOOKUP["H"]: 2})
    assert composition_to_formula_string(comp) == "H2O"


def test_proforma_negative_counts():
    """ProForma shows negative counts"""
    comp = Counter({ELEMENT_LOOKUP["H"]: -2, ELEMENT_LOOKUP["O"]: -1})
    assert composition_to_proforma_formula_string(comp) == "H-2O-1"


def test_formula_negative_counts():
    """Standard formula uses absolute values"""
    comp = Counter({ELEMENT_LOOKUP["H"]: -2, ELEMENT_LOOKUP["O"]: -1})
    assert composition_to_formula_string(comp) == "H2O"


def test_mixed_signs_raises():
    """Can't convert mixed positive/negative to standard formula"""
    comp = Counter({ELEMENT_LOOKUP["H"]: 2, ELEMENT_LOOKUP["O"]: -1})
    with pytest.raises(ValueError):
        composition_to_formula_string(comp)
