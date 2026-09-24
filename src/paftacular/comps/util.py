"""Utility functions for formula and composition conversions"""

from collections import Counter
from functools import lru_cache

from tacular import ELEMENT_LOOKUP, REFMOL_LOOKUP, UNIMOD_LOOKUP, ElementInfo, RefMolInfo

from ..errors import PafUnknownReferenceError
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


@lru_cache(maxsize=1024)
def lookup_reference(name: str) -> RefMolInfo:
    """Find a reference molecule by name.

    mzPAF 1.0.1 sections 4.4.7 and 4.5: the mzPAF reference molecule list takes priority,
    then a Unimod entry name (for example ``Hex`` or ``HexNAc(2)``).
    """
    try:
        return REFMOL_LOOKUP[name]
    except KeyError:
        pass
    try:
        unimod = UNIMOD_LOOKUP[name]
    except KeyError:
        unimod = None
    composition = Counter(unimod.composition) if unimod is not None and unimod.name == name and unimod.composition else None
    if not composition:
        raise PafUnknownReferenceError(name)
    if any(count < 0 for count in composition.values()):
        raise PafUnknownReferenceError(name, f"Unimod entry '{name}' is a composition change, not a molecule, so it cannot be a reference")
    return RefMolInfo(
        name=name,
        label_type="Unimod",
        molecule_type="modification",
        chemical_formula=composition_to_formula_string(composition),
        monoisotopic_mass=sum(element.get_mass(True) * count for element, count in composition.items()),
        average_mass=sum(element.get_mass(False) * count for element, count in composition.items()),
        dict_composition={str(element): count for element, count in composition.items()},
    )
