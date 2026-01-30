"""Utility functions for formula and composition conversions"""

from collections import Counter

from tacular import ELEMENT_LOOKUP, ElementInfo

from ..util import parse_formula


def formula_to_composition(formula: str) -> Counter[ElementInfo]:
    """Convert chemical formula string to elemental composition"""
    elem_counts: Counter[str] = parse_formula(formula)
    return Counter({ELEMENT_LOOKUP[elem]: count for elem, count in elem_counts.items()})


def composition_to_proforma_formula_string(comp: Counter[ElementInfo], hill_order: bool = True) -> str:
    """Convert composition to ProForma-style formula string"""
    keys = list(sorted(comp.keys()))
    comps = []
    for key in keys:
        elem_info, cnt = key, comp[key]
        if cnt == 0:
            continue
        comps.append(elem_info.serialize(cnt))
    return "".join(comps)


def composition_to_formula_string(comp: Counter[ElementInfo]) -> str:
    """Convert composition to standard chemical formula string"""
    if not comp:
        return ""

    all_positive = all(count >= 0 for count in comp.values())
    all_negative = all(count <= 0 for count in comp.values())
    if not (all_positive or all_negative):
        raise ValueError("Composition must have all positive or all negative counts to convert to formula string")

    keys = list(sorted(comp.keys()))
    comps = []
    for key in keys:
        elem_info, cnt = key, comp[key]
        if cnt == 0:
            continue
        comps.append(elem_info.serialize(abs(cnt)))
    return "".join(comps)
