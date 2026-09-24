"""mzPAF text parsing: parse, parse_multi and iter_parse."""

import re
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass

from .annotation import PafAnnotation
from .comps import (
    Adduct,
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    IonType,
    IsotopeSpecification,
    MassError,
    NamedCompound,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)
from .constants import (
    ADDUCT_REGEX_PATTERN,
    FULL_PAF_PATTERN,
    ISOTOPE_REGEX_PATTERN,
    MAX_CACHE_SIZE,
    NEUTRAL_LOSS_REGEX_PATTERN,
    PARTIAL_PAF_PATTERN,
    AminoAcids,
    IonSeries,
)
from .errors import PafParseError, PaftacularError
from .syntax import annotation_spans
from .util import to_enum

_ISOTOPE_TOKEN = re.compile(ISOTOPE_REGEX_PATTERN)
_NEUTRAL_LOSS_TOKEN = re.compile(NEUTRAL_LOSS_REGEX_PATTERN)
_ADDUCT_TOKEN = re.compile(ADDUCT_REGEX_PATTERN)


class _BoundedCache[V]:
    """Substring -> component cache. Components are immutable, so equal text shares one object.

    Only successful builds are stored, so an invalid substring fails every time it is parsed.
    The oldest entry is dropped once the cache holds ``MAX_CACHE_SIZE`` entries.
    """

    __slots__ = ("_build", "_data")

    def __init__(self, build: Callable[[str], V]):
        self._data: dict[str, V] = {}
        self._build = build

    def get(self, key: str) -> V:
        try:
            return self._data[key]
        except KeyError:
            pass
        value = self._build(key)
        data = self._data
        if len(data) >= MAX_CACHE_SIZE:
            del data[next(iter(data))]
        data[key] = value
        return value

    def clear(self) -> None:
        self._data.clear()

    def __len__(self) -> int:
        return len(self._data)


def _build_ion(groups: dict[str, str | None]) -> IonType:
    """Build the ion component from the regex groups of one annotation."""
    if series := groups["series"]:
        return PeptideIon(to_enum(IonSeries, series, "ion series"), int(groups["ordinal"] or ""), sequence=groups["sequence_ordinal"])
    if start := groups["internal_start"]:
        return InternalFragment(int(start), int(groups["internal_end"] or ""), sequence=groups["sequence_internal"])
    if groups["precursor"]:
        return PrecursorIon()
    if amino_acid := groups["immonium"]:
        return ImmoniumIon(to_enum(AminoAcids, amino_acid, "immonium amino acid"), modification=groups["immonium_modification"])
    if name := groups["reference_label"]:
        return ReferenceIon(name)
    if formula := groups["formula"]:
        return ChemicalFormula(formula)
    if name := groups["named_compound"]:
        return NamedCompound(name)
    if smiles := groups["smiles"]:
        return SMILESCompound(smiles)
    if groups["unannotated"]:
        label = groups["unannotated_label"]
        return UnknownIon(label=None if label is None else int(label))
    # The annotation regex guarantees one of the branches above.
    raise PaftacularError("Unable to parse ion type")


def _neutral_losses(text: str) -> tuple[NeutralLoss, ...]:
    return tuple(NeutralLoss.parse(token.group()) for token in _NEUTRAL_LOSS_TOKEN.finditer(text))


def _isotopes(text: str) -> tuple[IsotopeSpecification, ...]:
    isotopes: list[IsotopeSpecification] = []
    for sign, count_text, element in _ISOTOPE_TOKEN.findall(text):
        count = int(count_text) if count_text else 1
        if sign == "-":
            count = -count
        if element == "A":
            isotopes.append(IsotopeSpecification(count, is_average=True))
        else:
            isotopes.append(IsotopeSpecification(count, element=element or None))
    return tuple(isotopes)


def _adducts(text: str) -> tuple[Adduct, ...]:
    # The annotation regex guarantees the leading "M" and at least one token.
    adducts: list[Adduct] = []
    for sign, count_text, formula in _ADDUCT_TOKEN.findall(text[1:]):
        count = int(count_text) if count_text else 1
        adducts.append(Adduct(-count if sign == "-" else count, formula))
    return tuple(adducts)


def _mass_error(text: str) -> MassError:
    if text.endswith("ppm"):
        return MassError(float(text[:-3]), unit="ppm")
    return MassError(float(text))


# The ion component is built from the match groups, so it is cached inline by its text.
_ION_CACHE: dict[str, IonType] = {}
_LOSS_CACHE: _BoundedCache[tuple[NeutralLoss, ...]] = _BoundedCache(_neutral_losses)
_ISOTOPE_CACHE: _BoundedCache[tuple[IsotopeSpecification, ...]] = _BoundedCache(_isotopes)
_ADDUCT_CACHE: _BoundedCache[tuple[Adduct, ...]] = _BoundedCache(_adducts)
_MASS_ERROR_CACHE: _BoundedCache[MassError] = _BoundedCache(_mass_error)


def _clear_caches() -> None:
    """Empty every parser cache (for tests and memory measurements)."""
    _ION_CACHE.clear()
    for cache in (_LOSS_CACHE, _ISOTOPE_CACHE, _ADDUCT_CACHE, _MASS_ERROR_CACHE):
        cache.clear()


def _build_annotation(match: re.Match[str]) -> PafAnnotation:
    """Build a PafAnnotation from a match against FULL_PAF_PATTERN."""
    groups = match.groupdict()
    ion = groups["ion"] or ""
    ion_type = _ION_CACHE.get(ion)
    if ion_type is None:
        ion_type = _build_ion(groups)
        if len(_ION_CACHE) >= MAX_CACHE_SIZE:
            del _ION_CACHE[next(iter(_ION_CACHE))]
        _ION_CACHE[ion] = ion_type

    losses = groups["neutral_losses"]
    isotopes = groups["isotope"]
    adducts = groups["adducts"]
    charge = groups["charge"]
    analyte_reference = groups["analyte_reference"]
    mass_error = groups["mass_error"]
    confidence = groups["confidence"]
    if mass_error is not None and groups["mass_error_unit"]:
        mass_error += "ppm"
    return PafAnnotation(
        ion_type,
        analyte_reference=None if analyte_reference is None else int(analyte_reference),
        is_auxiliary=groups["is_auxiliary"] is not None,
        neutral_losses=_LOSS_CACHE.get(losses) if losses else (),
        isotopes=_ISOTOPE_CACHE.get(isotopes) if isotopes else (),
        adducts=_ADDUCT_CACHE.get(adducts) if adducts else (),
        charge=1 if charge is None else int(charge),
        mass_error=_MASS_ERROR_CACHE.get(mass_error) if mass_error else None,
        confidence=None if confidence is None else _float(confidence),
    )


def _float(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        raise PaftacularError(f"Expected a number, got {text!r}") from None


def _spans(text: str) -> Iterable[tuple[int, int]]:
    # Fast path: without a comma the whole text is one annotation. Delimiter errors are
    # reported by annotation_spans() when the regex match fails.
    if "," not in text:
        return ((0, len(text)),) if text.strip() else ()
    return annotation_spans(text)


def parse_multi(annotation_str: str) -> list[PafAnnotation]:
    """Parse annotations separated by commas outside labels and sequences.

    Returns an empty list for empty or blank input. Raises :class:`PafParseError`.
    """
    if not isinstance(annotation_str, str):
        raise TypeError(f"Expected an mzPAF string, got {type(annotation_str).__name__}")
    annotations: list[PafAnnotation] = []
    for index, (start, end) in enumerate(_spans(annotation_str)):
        segment = annotation_str[start:end]
        stripped = segment.strip()
        offset = start + len(segment) - len(segment.lstrip())
        match = FULL_PAF_PATTERN.fullmatch(stripped)
        if match is None:
            # Let the delimiter scan report unclosed brackets or braces first.
            for _ in annotation_spans(annotation_str):
                pass
            partial = PARTIAL_PAF_PATTERN.match(stripped)
            position = offset + (partial.end() if partial else 0)
            raise PafParseError(annotation_str, position, index, "Unexpected or missing annotation content")
        try:
            annotations.append(_build_annotation(match))
        except PaftacularError as error:
            raise PafParseError(annotation_str, offset, index, str(error)) from error
    return annotations


def parse(annotation_str: str) -> PafAnnotation:
    """Parse exactly one mzPAF annotation.

    Raises :class:`PafParseError` when the text holds zero or several annotations. Use
    :func:`parse_multi` for comma-separated input.
    """
    annotations = parse_multi(annotation_str)
    if len(annotations) != 1:
        raise PafParseError(annotation_str, 0, 0, f"Expected one annotation, got {len(annotations)}. Use parse_multi() for comma-separated annotations.")
    return annotations[0]


@dataclass(frozen=True, slots=True)
class ParseResult:
    """One input record and either its annotations or a structured error."""

    index: int
    text: str
    annotations: tuple[PafAnnotation, ...] = ()
    error: PafParseError | None = None

    @property
    def ok(self) -> bool:
        return self.error is None


def iter_parse(records: Iterable[str]) -> Iterator[ParseResult]:
    """Parse records lazily, retaining one error per failed input record."""
    if isinstance(records, str):
        raise TypeError("Expected an iterable of records, not one string")
    for index, text in enumerate(records):
        try:
            annotations = tuple(parse_multi(text))
        except PafParseError as error:
            yield ParseResult(index, text, error=error)
        else:
            yield ParseResult(index, text, annotations)
