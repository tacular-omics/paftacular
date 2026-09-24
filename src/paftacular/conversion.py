"""Convert peptacular fragments to PafAnnotation."""

from __future__ import annotations

import math
from collections import Counter
from functools import cache, lru_cache
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import peptacular as pt
else:
    try:
        import peptacular as pt
    except ImportError:
        pt = None

from tacular import FRAGMENT_ION_LOOKUP, NEUTRAL_DELTA_LOOKUP, FragmentIonInfo
from tacular import IonType as TacularIonType

from .annotation import PafAnnotation
from .comps import (
    Adduct,
    ImmoniumIon,
    InternalFragment,
    IonType,
    IsotopeSpecification,
    MassError,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    UnknownIon,
)
from .constants import INTERNAL_MASS_DIFFS, AminoAcids, IonSeries
from .errors import PaftacularError
from .util import parse_formula, to_enum

# Peptacular ion types that mzPAF writes as another series, plus the hydrogen change that keeps
# the mass. mzPAF z is the z-dot radical, so the Biemann z, z+H and c-H carry a hydrogen delta.
_SERIES_FOR_ION_TYPE: dict[TacularIonType, tuple[IonSeries, int]] = {
    TacularIonType.W_VALINE: (IonSeries.W, 0),
    TacularIonType.D_VALINE: (IonSeries.D, 0),
    TacularIonType.WB_ISOLEUCINE: (IonSeries.WB, 0),
    TacularIonType.WB_THREONINE: (IonSeries.WB, 0),
    TacularIonType.WA_ISOLEUCINE: (IonSeries.WA, 0),
    TacularIonType.WA_THREONINE: (IonSeries.WA, 0),
    TacularIonType.DB_ISOLEUCINE: (IonSeries.DB, 0),
    TacularIonType.DB_THREONINE: (IonSeries.DB, 0),
    TacularIonType.DA_ISOLEUCINE: (IonSeries.DA, 0),
    TacularIonType.DA_THREONINE: (IonSeries.DA, 0),
    TacularIonType.Z: (IonSeries.Z, -1),
    TacularIonType.Z_RADICAL: (IonSeries.Z, 0),
    TacularIonType.Z_PLUS_H: (IonSeries.Z, 1),
    TacularIonType.C_MINUS_H: (IonSeries.C, -1),
}
_HYDROGEN_LOSS = {-1: (NeutralLoss(-1, base_formula="H"),), 0: (), 1: (NeutralLoss(1, base_formula="H"),)}


def _require_peptacular() -> None:
    if pt is None:
        raise ImportError("peptacular is required for this feature. Install it with: pip install paftacular[peptacular]")


@lru_cache(maxsize=256)
def _parent(sequence: str) -> pt.ProFormaAnnotation:
    return pt.parse(sequence)


@lru_cache(maxsize=65536)
def _subsequence(parent: str, start: int, end: int) -> str:
    """The ProForma text of residues ``start:end`` (zero-based, end exclusive) without charge."""
    return _parent(parent)[start:end].serialize(exclude_charge=True)


@lru_cache(maxsize=4096)
def _whole_sequence(parent: str) -> str:
    return _parent(parent).serialize(exclude_charge=True)


def _fragment_span(info: FragmentIonInfo, frag: pt.Fragment) -> tuple[int, int] | None:
    """Zero-based residue span of the fragment in its parent, or None for the whole parent."""
    length = frag.parent_sequence_length
    if length is None:
        raise PaftacularError("Fragment has a parent sequence but no parent sequence length")
    position = frag.position
    if isinstance(position, int):
        if info.is_forward:
            span = (0, position)
        elif info.is_backward:
            span = (length - position, length)
        else:
            raise PaftacularError(f"Integer position {position} is not valid for ion type {info.ion_type}")
    elif isinstance(position, tuple):
        span = (position[0] - 1, position[1])
    else:
        return None
    if span[0] < 0 or span[1] > length or span[0] >= span[1]:
        raise PaftacularError(f"Fragment position {position!r} is outside its parent sequence of length {length}")
    return span


def _fragment_sequence(info: FragmentIonInfo, frag: pt.Fragment) -> str | None:
    """The fragment's own ProForma sequence without charge, or None without a parent sequence."""
    parent = frag.parent_sequence
    if parent is None:
        return None
    span = _fragment_span(info, frag)
    if span is None:
        return _whole_sequence(parent)
    return _subsequence(parent, span[0], span[1])


@cache
def _named_deltas() -> dict[frozenset[tuple[str, int]], str]:
    """Unsigned composition -> canonical mzPAF name of each known neutral delta (NH3, not H3N)."""
    named: dict[frozenset[tuple[str, int]], str] = {}
    for info in NEUTRAL_DELTA_LOOKUP.values():
        counts = parse_formula(info.formula)
        named.setdefault(frozenset(counts.items()), info.formula)
    return named


def _hill_key(has_carbon: bool):
    def key(item: tuple[str, int]) -> tuple[int, str]:
        symbol = item[0].lstrip("0123456789")
        if has_carbon and symbol in ("C", "H"):
            return (0 if symbol == "C" else 1, item[0])
        return (2, symbol + item[0])

    return key


@lru_cache(maxsize=4096)
def _formula_loss(formula: pt.ChargedFormula) -> NeutralLoss:
    """One unit of a peptacular delta formula as a NeutralLoss (count +1 or -1)."""
    counts: Counter[str] = Counter()
    for element in formula.formula:
        symbol = element.element.value
        counts[f"{element.isotope}{symbol}" if element.isotope else symbol] += element.occurance
    signs = {count > 0 for count in counts.values() if count}
    if len(signs) != 1:
        raise PaftacularError(f"Cannot write delta {formula} in mzPAF: it must be a pure gain or a pure loss")
    sign = 1 if signs.pop() else -1
    magnitudes = {key: abs(count) for key, count in counts.items() if count}
    name = _named_deltas().get(frozenset(magnitudes.items()))
    if name is None:
        parts: list[str] = []
        for key, count in sorted(magnitudes.items(), key=_hill_key("C" in magnitudes)):
            if key[0].isdigit():
                parts.append(f"[{key}{count}]")
            else:
                parts.append(f"{key}{count if count > 1 else ''}")
        name = "".join(parts)
    return NeutralLoss(sign, base_formula=name)


def _losses(frag: pt.Fragment) -> list[NeutralLoss]:
    losses: list[NeutralLoss] = []
    for key, count in frag.deltas.items():
        if isinstance(key, int | float):
            total = float(key) * count
            if not math.isfinite(total):
                raise PaftacularError(f"Fragment delta {key!r} x {count} is not a finite mass")
            # Rounded like peptacular's own label: 6 decimals, a delta that rounds to zero is left out.
            total = round(total, 6)
            if total != 0:
                losses.append(NeutralLoss(1 if total > 0 else -1, base_mass=abs(total)))
            continue
        unit = _formula_loss(key)
        losses.append(unit if count == 1 else NeutralLoss(unit.count * count, base_formula=unit.base_formula))
    return losses


def _isotopes(frag: pt.Fragment) -> tuple[IsotopeSpecification, ...]:
    isotopes = frag.isotopes
    if not isotopes:
        return ()
    if frag.is_c13:
        # A plain isotope count is the generic 13C-12C shift (mzPAF "+ni").
        return tuple(IsotopeSpecification(count) for count in isotopes.values() if count)
    return tuple(IsotopeSpecification(count, element=str(element)) for element, count in isotopes.items() if count)


@lru_cache(maxsize=1024)
def _adduct(text: str) -> Adduct:
    return Adduct.parse(text)


def _adducts(frag: pt.Fragment) -> tuple[Adduct, ...]:
    if frag.is_protonated:
        return ()
    adducts: list[Adduct] = []
    for mod in frag.charge_adducts.mods:
        carrier: pt.GlobalChargeCarrier = mod.value
        adduct = _adduct(carrier.to_mz_paf()[1:])  # "M+Na" -> "+Na"
        if mod.count != 1:
            adduct = Adduct(adduct.count * mod.count, adduct.base_formula)
        adducts.append(adduct)
    # mzPAF 4.7: several carriers SHOULD be in alphabetical order ([M+2H+Na]).
    adducts.sort(key=lambda adduct: adduct.base_formula)
    return tuple(adducts)


def _immonium(frag: pt.Fragment) -> ImmoniumIon:
    sequence = frag.sequence if frag.parent_sequence is not None else None
    if sequence is None:
        raise PaftacularError("An immonium fragment needs a parent sequence")
    annot = pt.parse(sequence)
    if len(annot.sequence) != 1:
        raise PaftacularError(f"Immonium ion sequence must be a single amino acid, got {annot.sequence}")
    # mzPAF allows one modification on an immonium ion. A terminal modification of the
    # residue (for example an N-terminal acetyl) adds the same mass, so it is written there.
    tags: list[str] = []
    for has_mods, get_mods in (
        (annot.has_internal_mods_at_index(0), lambda: annot.get_internal_mods_at_index(0)),
        (annot.has_nterm_mods, lambda: annot.nterm_mods),
        (annot.has_cterm_mods, lambda: annot.cterm_mods),
    ):
        if has_mods:
            for mod in get_mods().mods:
                tags.extend([str(mod.value)] * mod.count)
    if len(tags) > 1:
        raise PaftacularError(f"mzPAF allows one modification on an immonium ion, got {', '.join(tags)}")
    modification = tags[0] if tags else None
    return ImmoniumIon(to_enum(AminoAcids, annot.sequence, "immonium amino acid"), modification=modification)


def to_mzpaf(
    frag: pt.Fragment,
    *,
    confidence: float | None = None,
    mass_error: float | None = None,
    mass_error_type: Literal["ppm", "da"] = "ppm",
    include_sequence: bool = True,
) -> PafAnnotation:
    """Convert a peptacular Fragment to a PafAnnotation.

    Mass deltas are folded (``mass * count``) and rounded to 6 decimals like
    ``Fragment.to_mzpaf()``. A delta that rounds to zero is left out. A negative charge state
    is kept (``serialize()`` writes ``^-n``). With ``include_sequence`` the fragment's own
    sequence is embedded (peptide and internal ions) or stored as the resolved sequence
    (precursor ions).
    """
    _require_peptacular()

    sequence: str | None = None
    ion: IonType
    ion_type = frag.ion_type
    if ion_type is None:
        ion = UnknownIon()
        fixed_losses: tuple[NeutralLoss, ...] = ()
        kind = _UNKNOWN
    else:
        kind, series, fixed_losses = _plan(ion_type)
        if kind == _TERMINAL:
            position = frag.position
            if include_sequence:
                sequence = _fragment_sequence(FRAGMENT_ION_LOOKUP[ion_type], frag)
            assert series is not None
            ion = _peptide_ion(series, position if isinstance(position, int) else -1, sequence)
        elif kind == _IMMONIUM:
            ion = _immonium(frag)
        elif kind == _INTERNAL:
            start, end = frag.position if isinstance(frag.position, tuple) and len(frag.position) == 2 else (-1, -1)
            if include_sequence:
                sequence = _fragment_sequence(FRAGMENT_ION_LOOKUP[ion_type], frag)
            ion = InternalFragment(start, end, sequence=sequence)
        else:
            ion = PrecursorIon()
            if include_sequence:
                sequence = _fragment_sequence(FRAGMENT_ION_LOOKUP[ion_type], frag)

    charge = frag.charge_state
    if charge == 0:
        raise PaftacularError("Cannot write an uncharged fragment in mzPAF")
    deltas = _losses(frag)
    if kind == _INTERNAL:
        losses = (*deltas, *fixed_losses)
    elif fixed_losses:
        losses = (*fixed_losses, *deltas)
    else:
        losses = tuple(deltas)

    return PafAnnotation(
        ion,
        neutral_losses=losses,
        isotopes=_isotopes(frag),
        adducts=_adducts(frag),
        charge=charge,
        mass_error=None if mass_error is None else MassError(mass_error, unit=mass_error_type),
        confidence=confidence,
        resolved_sequence=sequence if kind == _PRECURSOR else None,
    )


_UNKNOWN, _TERMINAL, _IMMONIUM, _INTERNAL, _PRECURSOR = range(5)


@lru_cache(maxsize=128)
def _plan(ion_type: str) -> tuple[int, IonSeries | None, tuple[NeutralLoss, ...]]:
    """How to convert one peptacular ion type: (kind, mzPAF series, fixed losses).

    Cached per ion type, so the tacular flag checks run once per type, not once per fragment.
    """
    info = FRAGMENT_ION_LOOKUP[ion_type]
    if info.is_forward or info.is_backward:
        mapped = _SERIES_FOR_ION_TYPE.get(info.ion_type)
        series, hydrogen = mapped if mapped is not None else (to_enum(IonSeries, info.ion_type.value, "ion series"), 0)
        return _TERMINAL, series, _HYDROGEN_LOSS[hydrogen]
    if info.is_internal and info.ion_type == TacularIonType.IMMONIUM:
        return _IMMONIUM, None, ()
    if info.is_internal:
        key = tuple(info.ion_type.value)
        if len(key) != 2 or key not in INTERNAL_MASS_DIFFS:
            raise PaftacularError(f"Internal ion type {info.ion_type} is not supported in mzPAF")
        return _INTERNAL, None, _internal_losses(key[0], key[1])
    if info.is_intact and info.ion_type == TacularIonType.PRECURSOR:
        return _PRECURSOR, None, ()
    raise PaftacularError(f"Cannot convert fragment with ion type {ion_type} to mzPAF")


@lru_cache(maxsize=65536)
def _peptide_ion(series: IonSeries, position: int, sequence: str | None) -> PeptideIon:
    # Components are immutable, so equal ions share one object.
    return PeptideIon(series, position, sequence=sequence)


@lru_cache(maxsize=16)
def _internal_losses(nterm: str, cterm: str) -> tuple[NeutralLoss, ...]:
    """The physical correction for an internal cleavage, relative to by, as neutral losses."""
    correction = InternalFragment(1, 2, nterm_ion_type=IonSeries(nterm), cterm_ion_type=IonSeries(cterm)).cleavage_correction
    if not correction:
        return ()
    from .parser import _NEUTRAL_LOSS_TOKEN

    return tuple(NeutralLoss.parse(token.group()) for token in _NEUTRAL_LOSS_TOKEN.finditer(correction))
