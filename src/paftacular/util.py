import re
from collections import Counter


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
        raise ValueError("Empty formula string")

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
                raise ValueError(f"Unclosed bracket at position {i}")

            # Extract content inside brackets: "13C2" or "13C"
            content = formula[i + 1 : close]

            # Parse: isotope_number + element + optional_count
            # Pattern: digits followed by element (capital + optional lowercase) + optional digits
            match = re.match(r"^(\d+)([A-Z][a-z]?)(\d*)$", content)
            if not match:
                raise ValueError(f"Invalid isotope format: [{content}]")

            isotope_num, element, count_str = match.groups()
            count = int(count_str) if count_str else 1

            # Use isotope notation WITHOUT brackets as key: 13C
            element_key = f"{isotope_num}{element}"
            element_counts[element_key] += count

            i = close + 1

        # Handle regular element: C2, Ca, H
        elif formula[i].isupper():
            # Get element symbol (capital + optional lowercase)
            element = formula[i]
            i += 1
            if i < len(formula) and formula[i].islower():
                element += formula[i]
                i += 1

            # Get optional count
            count_str = ""
            while i < len(formula) and formula[i].isdigit():
                count_str += formula[i]
                i += 1

            count = int(count_str) if count_str else 1
            element_counts[element] += count

        else:
            raise ValueError(f"Unexpected character '{formula[i]}' at position {i}")

    if not element_counts:
        raise ValueError(f"No elements found in formula: '{formula}'")

    return element_counts
