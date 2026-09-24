import re
from collections import Counter
from decimal import Decimal
from enum import Enum
from math import isfinite

from .errors import PaftacularError


def to_enum[E: Enum](cls: type[E], value: object, what: str) -> E:
    """Convert ``value`` to the enum ``cls``, raising :class:`PaftacularError` naming ``what``."""
    try:
        return cls(value)
    except ValueError:
        choices = ", ".join(str(member.value) for member in cls)  # type: ignore[attr-defined]
        raise PaftacularError(f"Invalid {what} {value!r}. Expected one of: {choices}") from None


def validate_number(value: float) -> None:
    """Require a number representable as a finite float, excluding booleans."""
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PaftacularError("Expected a finite number")
    try:
        finite = isfinite(value)
    except OverflowError:
        finite = False
    if not finite:
        raise PaftacularError("Expected a finite number")


def format_number(value: float, *, minimum_places: int = 0) -> str:
    """Render a finite number as a decimal that round trips through float."""
    validate_number(value)
    whole, _, fraction = format(Decimal(str(value)), "f").partition(".")
    fraction = fraction.rstrip("0").ljust(minimum_places, "0")
    return whole + ("." + fraction if fraction else "")


def validate_integer(value: int, name: str, *, minimum: int = 0) -> None:
    """Reject booleans, nonintegral types, and values below the minimum."""
    if type(value) is not int or value < minimum:
        raise PaftacularError(f"{name} must be an integer >= {minimum}, got {value!r}")


def parse_formula(formula: str) -> Counter[str]:
    """
    Parse a chemical formula into element counts, supporting isotopes.

    Args:
        formula: Chemical formula string (e.g., "H2O", "CO2", "[13C2]H6")

    Returns:
        Counter mapping element symbols to their counts
        - Regular elements: "H", "O", "Ca"
        - Isotopes: "13C", "2H" (without brackets)
    """
    if not formula:
        raise PaftacularError("Empty formula string")

    element_counts = Counter()
    i = 0

    while i < len(formula):
        # Skip whitespace
        if formula[i].isspace():
            i += 1
            continue

        # Handle isotope notation: [13C2] or [13C]
        if formula[i] == "[":
            close = formula.find("]", i)
            if close == -1:
                raise PaftacularError(f"Unclosed bracket at position {i}")

            # Extract content inside brackets: "13C2" or "13C"
            content = formula[i + 1 : close]

            # Parse: isotope_number + element + optional_count
            # Pattern: digits followed by element (capital + optional lowercase) + optional digits
            match = re.fullmatch(r"([0-9]+)([A-Z][a-z]?)([0-9]*)", content)
            if not match:
                raise PaftacularError(f"Invalid isotope format: [{content}]")

            isotope_num, element, count_str = match.groups()
            count = int(count_str) if count_str else 1

            # Use isotope notation WITHOUT brackets as key: 13C
            element_key = f"{isotope_num}{element}"
            element_counts[element_key] += count

            i = close + 1

        # Handle regular element: C2, Ca, H
        elif "A" <= formula[i] <= "Z":
            # Get element symbol (capital + optional lowercase)
            element = formula[i]
            i += 1
            if i < len(formula) and "a" <= formula[i] <= "z":
                element += formula[i]
                i += 1

            # Get optional count
            count_str = ""
            # ASCII digits only: str.isdigit() also accepts superscripts that int() rejects.
            while i < len(formula) and formula[i] in "0123456789":
                count_str += formula[i]
                i += 1

            count = int(count_str) if count_str else 1
            element_counts[element] += count

        else:
            raise PaftacularError(f"Unexpected character '{formula[i]}' at position {i}")

    if not element_counts:
        raise PaftacularError(f"No elements found in formula: '{formula}'")

    return element_counts
