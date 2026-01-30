"""Tests for utility functions in util.py"""

from collections import Counter

import pytest

from paftacular.util import parse_formula


class TestParseFormula:
    """Test parse_formula function"""

    def test_simple_formula(self):
        """Test parsing a simple formula"""
        result = parse_formula("H2O")
        assert result == Counter({"H": 2, "O": 1})

    def test_complex_formula(self):
        """Test parsing a more complex formula"""
        result = parse_formula("C6H12O6")
        assert result == Counter({"C": 6, "H": 12, "O": 6})

    def test_single_element(self):
        """Test parsing a single element"""
        result = parse_formula("C")
        assert result == Counter({"C": 1})

    def test_element_with_count(self):
        """Test parsing an element with count"""
        result = parse_formula("O2")
        assert result == Counter({"O": 2})

    def test_two_letter_element(self):
        """Test parsing two-letter element symbols"""
        result = parse_formula("CaCl2")
        assert result == Counter({"Ca": 1, "Cl": 2})

    def test_isotope_notation(self):
        """Test parsing isotope notation"""
        result = parse_formula("[13C]H4")
        assert result == Counter({"13C": 1, "H": 4})

    def test_isotope_with_count(self):
        """Test parsing isotope with count"""
        result = parse_formula("[13C2]H6")
        assert result == Counter({"13C": 2, "H": 6})

    def test_multiple_isotopes(self):
        """Test parsing multiple isotopes"""
        result = parse_formula("[13C2][15N]H5")
        assert result == Counter({"13C": 2, "15N": 1, "H": 5})

    def test_mixed_regular_and_isotopes(self):
        """Test parsing mix of regular elements and isotopes"""
        result = parse_formula("C6[13C2]H14O2")
        assert result == Counter({"C": 6, "13C": 2, "H": 14, "O": 2})

    def test_whitespace_handling(self):
        """Test that whitespace is handled correctly"""
        result = parse_formula("H2 O")
        assert result == Counter({"H": 2, "O": 1})

    def test_empty_formula_raises(self):
        """Test that empty formula raises ValueError"""
        with pytest.raises(ValueError):
            parse_formula("")

    def test_unclosed_bracket_raises(self):
        """Test that unclosed bracket raises ValueError"""
        with pytest.raises(ValueError):
            parse_formula("[13C")

    def test_invalid_isotope_format_raises(self):
        """Test that invalid isotope format raises ValueError"""
        with pytest.raises(ValueError):
            parse_formula("[ABC]")

    def test_unexpected_character_raises(self):
        """Test that unexpected character raises ValueError"""
        with pytest.raises(ValueError):
            parse_formula("H2$O")

    def test_no_elements_found_raises(self):
        """Test that formula with no valid elements raises ValueError"""
        with pytest.raises(ValueError):
            parse_formula("   ")

    def test_isotope_two_letter_element(self):
        """Test parsing isotope with two-letter element"""
        result = parse_formula("[14Ca]O")
        assert result == Counter({"14Ca": 1, "O": 1})

    def test_large_counts(self):
        """Test parsing formulas with large counts"""
        result = parse_formula("C100H200")
        assert result == Counter({"C": 100, "H": 200})

    def test_isotope_large_count(self):
        """Test isotope with large count"""
        result = parse_formula("[13C10]")
        assert result == Counter({"13C": 10})
