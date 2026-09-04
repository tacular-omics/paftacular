import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import ClassVar

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
    NEUTRAL_LOSS_REGEX_PATTERN,
    PARTIAL_PAF_PATTERN,
    AminoAcids,
    IonSeries,
)
from .errors import PafParseError
from .syntax import annotation_spans


class mzPAFParser:
    _instance: ClassVar["mzPAFParser | None"] = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = object.__new__(cls)
        return cls._instance

    def parse(self, annotation_str: str) -> PafAnnotation:
        """Parse a single annotation string"""
        return parse_single(annotation_str)

    def _build_annotation(self, match: re.Match[str]) -> PafAnnotation:
        """Build a PafAnnotation from a match against FULL_PAF_PATTERN/PARTIAL_PAF_PATTERN"""
        groups = match.groupdict()

        return PafAnnotation(
            ion_type=self._parse_ion_type(groups),
            analyte_reference=self._parse_int(groups, "analyte_reference"),
            is_auxiliary=bool(groups.get("is_auxiliary")),
            neutral_losses=self._parse_neutral_losses(groups.get("neutral_losses")),
            isotopes=self._parse_isotopes(groups.get("isotope")),
            adducts=self._parse_adducts(groups.get("adducts")),
            charge=1 if groups.get("charge") is None else self._require_int(groups, "charge"),
            mass_error=self._parse_mass_error(groups),
            confidence=self._parse_float(groups, "confidence"),
        )

    def _parse_ion_type(self, groups: dict[str, str | None]) -> IonType:
        """Parse the ion type from regex groups using dispatch pattern"""

        # Peptide ion series
        if groups.get("series"):
            series_str = self._require(groups, "series")
            ion_series = IonSeries(series_str)
            return PeptideIon(
                series=ion_series,
                position=self._require_int(groups, "ordinal"),
                sequence=groups.get("sequence_ordinal"),
            )

        # Internal fragment
        if groups.get("internal_start"):
            return InternalFragment(
                start_position=self._require_int(groups, "internal_start"),
                end_position=self._require_int(groups, "internal_end"),
                sequence=groups.get("sequence_internal"),
            )

        # Precursor ion
        if groups.get("precursor"):
            return PrecursorIon()

        # Immonium ion
        if groups.get("immonium"):
            aa_str = self._require(groups, "immonium")
            amino_acid = AminoAcids(aa_str)
            return ImmoniumIon(
                amino_acid=amino_acid,
                modification=groups.get("immonium_modification"),
            )

        # Reference ion
        if groups.get("reference_label"):
            return ReferenceIon(name=self._require(groups, "reference_label"))

        # Chemical formula
        if groups.get("formula"):
            return ChemicalFormula(formula=self._require(groups, "formula"))

        # Named compound
        if groups.get("named_compound"):
            return NamedCompound(name=self._require(groups, "named_compound"))

        # SMILES compound
        if groups.get("smiles"):
            return SMILESCompound(smiles=self._require(groups, "smiles"))

        # Unknown/unannotated ion
        if groups.get("unannotated"):
            return UnknownIon(label=self._parse_int(groups, "unannotated_label"))

        # Should never reach here due to regex, but provide helpful error
        non_null = {k: v for k, v in groups.items() if v is not None}
        raise ValueError(f"Unable to parse ion type. Available groups: {non_null}")

    def _parse_isotopes(self, isotope_str: str | None) -> tuple[IsotopeSpecification, ...]:
        """Parse isotope notation into tuple of specifications

        Examples:
            "+i" -> (IsotopeSpecification(count=1),)
            "-2i13C" -> (IsotopeSpecification(count=-2, element="13C"),)
            "+i-2i13C+iA" -> (IsotopeSpecification(count=1), IsotopeSpecification(count=-2, element="13C"), IsotopeSpecification(count=1, is_average=True))
        """
        if not isotope_str:
            return ()

        # Extract individual isotope strings like "+i", "-2i13C", "+iA"
        isotope_matches = re.findall(ISOTOPE_REGEX_PATTERN, isotope_str)
        if not isotope_matches:
            return ()

        # Parse each isotope component
        isotopes: list[IsotopeSpecification] = []
        for match_groups in isotope_matches:
            # The outer annotation pattern has already validated these tokens.
            sign_str, count_str, element_or_avg = match_groups

            count = int(count_str) if count_str else 1
            isotopes.append(
                IsotopeSpecification(
                    count=-count if sign_str == "-" else count,
                    element=None if element_or_avg == "A" else element_or_avg or None,
                    is_average=element_or_avg == "A",
                )
            )

        return tuple(isotopes)

    def _parse_neutral_losses(self, losses_str: str | None) -> tuple[NeutralLoss, ...]:
        """Parse neutral losses/gains into tuple of NeutralLoss objects"""
        if not losses_str:
            return ()

        # Pattern matches: [+-] followed by number, formula, or named group

        loss_strings = re.findall(NEUTRAL_LOSS_REGEX_PATTERN, losses_str)

        losses: list[NeutralLoss] = []
        for loss_str in loss_strings:
            losses.append(NeutralLoss.parse(loss_str))
        return tuple(losses)

    def _parse_adducts(self, adduct_str: str | None) -> tuple[Adduct, ...]:
        """Parse adduct notation into tuple of Adduct objects

        Examples:
            "M+H" -> (Adduct(count=1, base_formula="H"),)
            "M+2H+Na" -> (Adduct(count=2, base_formula="H"), Adduct(count=1, base_formula="Na"))
            "M+NH4" -> (Adduct(count=1, base_formula="NH4"),)
            "M-H+2Na" -> (Adduct(count=-1, base_formula="H"), Adduct(count=2, base_formula="Na"))
        """
        if not adduct_str:
            return ()

        # Verify M prefix
        if not adduct_str.startswith("M"):
            raise ValueError(f"Adduct string must start with 'M': '{adduct_str}'")

        content = adduct_str[1:]  # Remove 'M'
        if not content:
            raise ValueError(f"Adduct string must have components after 'M': '{adduct_str}'")

        # Extract individual adduct strings like "+H", "+2Na", "-NH4"
        adduct_strings = re.findall(ADDUCT_REGEX_PATTERN, content)
        if not adduct_strings:
            raise ValueError(f"No adduct components found in '{adduct_str}'")

        # Parse each adduct component
        adducts: list[Adduct] = []
        for match_groups in adduct_strings:
            # Construct directly from the validated token groups.
            sign_str, count_str, formula = match_groups
            count = int(count_str) if count_str else 1
            adducts.append(Adduct(count=-count if sign_str == "-" else count, base_formula=formula))

        return tuple(adducts)

    def _parse_mass_error(self, groups: dict[str, str | None]) -> MassError | None:
        """Parse mass error value and unit"""
        if not groups.get("mass_error"):
            return None

        mass_error_str = groups.get("mass_error")
        if mass_error_str is None:
            raise ValueError("Mass error value is missing")
        value = float(mass_error_str)
        unit = groups.get("mass_error_unit")

        if unit == "ppm":
            return MassError(value, "ppm")
        elif unit is None:
            return MassError(value, "da")
        else:
            raise ValueError(f"Unknown mass error unit: '{unit}'")

    # Helper methods for common parsing patterns
    def _require(self, groups: dict[str, str | None], key: str) -> str:
        """Get required string value from groups, raise if missing"""
        value = groups.get(key)
        if value is None or not isinstance(value, str):
            raise ValueError(f"Required field '{key}' is missing or invalid")
        return value

    def _require_int(self, groups: dict[str, str | None], key: str) -> int:
        """Get required integer value from groups, raise if missing"""
        value = groups.get(key)
        if value is None:
            raise ValueError(f"Required field '{key}' is missing")
        try:
            return int(value)
        except ValueError as e:
            raise ValueError(f"Field '{key}' must be an integer, got '{value}'") from e

    def _parse_int(self, groups: dict[str, str | None], key: str) -> int | None:
        """Parse optional integer from groups"""
        value = groups.get(key)
        if value is None:
            return None
        try:
            return int(value)
        except ValueError as e:
            raise ValueError(f"Field '{key}' must be an integer, got '{value}'") from e

    def _parse_float(self, groups: dict[str, str | None], key: str) -> float | None:
        """Parse optional float from groups"""
        value = groups.get(key)
        if value is None:
            return None
        try:
            return float(value)
        except ValueError as e:
            raise ValueError(f"Field '{key}' must be a number, got '{value}'") from e

    def parse_multi(self, annotation_str: str) -> list[PafAnnotation]:
        """Parse annotations separated by commas outside labels and sequences."""
        annotations: list[PafAnnotation] = []
        for index, (start, end) in enumerate(annotation_spans(annotation_str)):
            segment = annotation_str[start:end]
            offset = start + len(segment) - len(segment.lstrip())
            segment = segment.strip()
            match = FULL_PAF_PATTERN.fullmatch(segment)
            if match is None:
                partial = PARTIAL_PAF_PATTERN.match(segment)
                position = offset + (partial.end() if partial else 0)
                raise PafParseError(annotation_str, position, index, "Unexpected or missing annotation content")
            try:
                annotations.append(self._build_annotation(match))
            except ValueError as error:
                raise PafParseError(annotation_str, offset, index, str(error)) from error
        return annotations


MZ_PAF_PARSER = mzPAFParser()


def parse_multi(annotation_str: str) -> list[PafAnnotation]:
    """parse mzPAF annotation string into list of PafAnnotation"""
    return MZ_PAF_PARSER.parse_multi(annotation_str)


def parse(annotation_str: str) -> PafAnnotation | list[PafAnnotation]:
    """parse single mzPAF annotation string into PafAnnotation"""
    annots = parse_multi(annotation_str)
    if len(annots) == 1:
        return annots[0]
    return annots


def parse_single(annotation_str: str) -> PafAnnotation:
    """Parse exactly one annotation."""
    annots = parse_multi(annotation_str)
    if len(annots) != 1:
        raise PafParseError(annotation_str, 0, 0, f"Expected single annotation, got {len(annots)}")
    return annots[0]


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


def parse_batch(records: Iterable[str]) -> list[ParseResult]:
    """Collect batch results without discarding malformed input records."""
    return list(iter_parse(records))
