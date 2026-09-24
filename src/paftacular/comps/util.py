"""Utility functions for formula and composition conversions"""

from collections import Counter
from functools import lru_cache

from tacular import ELEMENT_LOOKUP, REFMOL_LOOKUP, UNIMOD_LOOKUP, ElementInfo, RefMolInfo

from ..errors import PaftacularError, PafUnknownReferenceError
from ..util import parse_formula


@lru_cache(maxsize=4096)
def _formula_items(formula: str) -> tuple[tuple[ElementInfo, int], ...]:
    elem_counts: Counter[str] = parse_formula(formula)
    try:
        return tuple((ELEMENT_LOOKUP[elem], count) for elem, count in elem_counts.items())
    except KeyError as error:
        raise PaftacularError(f"Unknown element or isotope in formula {formula!r}: {error}") from None


def formula_to_composition(formula: str) -> Counter[ElementInfo]:
    """Convert chemical formula string to elemental composition (a new Counter each call)."""
    return Counter(dict(_formula_items(formula)))


def composition_to_proforma_formula_string(comp: Counter[ElementInfo]) -> str:
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
        raise PaftacularError("Composition must have all positive or all negative counts to convert to formula string")

    keys = list(sorted(comp.keys()))
    comps = []
    for key in keys:
        elem_info, cnt = key, comp[key]
        if cnt == 0:
            continue
        comps.append(elem_info.serialize(abs(cnt)))
    return "".join(comps)


def _listed_or_summed(listed: float | None, composition: Counter[ElementInfo], *, monoisotopic: bool) -> float:
    """The listed database mass, or the composition's mass when the entry lists none."""
    if listed is not None:
        return listed
    return sum(element.get_mass(monoisotopic=monoisotopic) * count for element, count in composition.items())


@lru_cache(maxsize=1024)
def lookup_reference(name: str) -> RefMolInfo:
    """Find a reference molecule by name.

    mzPAF 1.0.1 sections 4.4.7 and 4.5: the mzPAF reference molecule list takes priority,
    then a Unimod entry name (for example ``Hex`` or ``HexNAc(2)``). A Unimod entry keeps its
    listed masses, like peptacular, and its composition only supplies the formula.
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
    if unimod is None or not composition:
        raise PafUnknownReferenceError(name)
    if any(count < 0 for count in composition.values()):
        raise PafUnknownReferenceError(name, f"Unimod entry '{name}' is a composition change, not a molecule, so it cannot be a reference")
    return RefMolInfo(
        name=name,
        label_type="Unimod",
        molecule_type="modification",
        formula=composition_to_formula_string(composition),
        monoisotopic_mass=_listed_or_summed(unimod.monoisotopic_mass, composition, monoisotopic=True),
        average_mass=_listed_or_summed(unimod.average_mass, composition, monoisotopic=False),
        dict_composition={str(element): count for element, count in composition.items()},
    )
