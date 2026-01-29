import re
from collections import Counter
from dataclasses import asdict, dataclass

from tacular import ELEMENT_LOOKUP, ElementInfo

from .comps import (
    Adduct,
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    IsotopeSpecification,
    MassError,
    NamedCompound,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
    composition_to_formula_string,
    composition_to_proforma_formula_string,
)
from .constants import AminoAcids, IonSeries

# Type aliases for cleaner code
IonType = (
    PeptideIon
    | InternalFragment
    | ImmoniumIon
    | ReferenceIon
    | NamedCompound
    | ChemicalFormula
    | SMILESCompound
    | UnknownIon
    | PrecursorIon
)


@dataclass(frozen=True, slots=True)
class PafAnnotation:
    """Complete fragment ion annotation following mzPAF specification"""

    # Core ion description
    ion_type: IonType

    # Optional components
    analyte_reference: int | None = None
    is_auxiliary: bool = False
    neutral_losses: tuple[NeutralLoss, ...] = ()  # Immutable tuple instead of mutable list
    isotopes: tuple[IsotopeSpecification, ...] = ()  # Renamed from 'isotope' for clarity
    adducts: tuple[Adduct, ...] = ()  # Immutable tuple
    charge: int = 1
    mass_error: MassError | None = None
    confidence: float | None = None

    def __post_init__(self):
        """Validate annotation constraints"""
        if self.charge < 1:
            raise ValueError(f"Charge must be >= 1, got {self.charge}")
        if self.confidence is not None and not (0.0 <= self.confidence <= 1.0):
            raise ValueError(f"Confidence must be between 0.0 and 1.0, got {self.confidence}")

    def mass(self, monoisotopic: bool = True, calculate_sequence: bool = False) -> float:
        """Calculate the mass of the annotated ion including modifications"""
        base_mass = self.ion_type.mass(monoisotopic=monoisotopic)

        # Apply neutral losses/gains
        for loss in self.neutral_losses:
            base_mass += loss.mass(monoisotopic=monoisotopic)

        # Apply adducts
        for adduct in self.adducts:
            base_mass += adduct.mass(monoisotopic=monoisotopic)

        # Adjust for charge state (if no adducts specified) default protonation/deprotonation
        if self.charge != 0 and len(self.adducts) == 0:
            base_mass += self.charge * 1.007276466812

        if calculate_sequence is True and self.sequence is not None:
            # Additional mass calculations based on sequence can be added here
            import peptacular as pt  # ty: ignore

            annot = pt.parse(self.sequence)

            if annot.has_charge:
                raise ValueError("Sequence in annotation should not have charge for mass calculation")

            sequence_mass = annot.mass(monoisotopic=monoisotopic, ion_type="n")
            base_mass += sequence_mass

        return base_mass

    def mz(self, monoisotopic: bool = True, calculate_sequence: bool = False) -> float:
        """Calculate the m/z of the annotated ion"""
        total_mass = self.mass(monoisotopic=monoisotopic, calculate_sequence=calculate_sequence)
        return total_mass / self.charge

    def composition(self, calculate_sequence: bool = False) -> Counter[ElementInfo]:
        """Calculate the elemental composition of the annotated ion including modifications"""
        comp: Counter[ElementInfo] = Counter()

        # Base ion composition
        comp.update(self.ion_type.composition)

        # Apply neutral losses/gains
        for loss in self.neutral_losses:
            comp.update(loss.composition)

        # Apply adducts
        for adduct in self.adducts:
            comp.update(adduct.composition)

        # Adjust for charge state (if no adducts specified) default protonation/deprotonation
        if self.charge != 0 and len(self.adducts) == 0:
            proton = ELEMENT_LOOKUP["H"]
            comp[proton] += self.charge

        if calculate_sequence is True and self.sequence is not None:
            # Additional composition calculations based on sequence can be added here
            import peptacular as pt  # ty: ignore

            annot = pt.parse(self.sequence)

            if annot.has_charge:
                raise ValueError("Sequence in annotation should not have charge for mass calculation")

            seq_comp = annot.comp()
            comp.update(seq_comp)

        return comp

    def dict_composition(self, calculate_sequence: bool = False) -> dict[str, int]:
        """Get the elemental composition as a dictionary of element symbols to counts"""
        comp_counter = self.composition(calculate_sequence=calculate_sequence)
        return {str(elem): count for elem, count in comp_counter.items()}

    @property
    def sequence(self) -> str | None:
        """Get the peptide sequence if applicable, else None"""
        if isinstance(self.ion_type, PeptideIon):
            return self.ion_type.sequence
        elif isinstance(self.ion_type, InternalFragment):
            return self.ion_type.sequence
        return None

    def formula(self, calculate_sequence: bool = False) -> str:
        """Get the chemical formula string of the annotated ion"""
        return composition_to_formula_string(self.composition(calculate_sequence=calculate_sequence))

    def proforma_formula(self, calculate_sequence: bool = False) -> str:
        """Get the ProForma-style chemical formula string of the annotated ion"""
        return composition_to_proforma_formula_string(self.composition(calculate_sequence=calculate_sequence))

    def serialize(self) -> str:
        """Serialize the annotation back to mzPAF string format"""
        parts: list[str] = []

        # Auxiliary marker
        if self.is_auxiliary:
            parts.append("&")

        # Analyte reference
        if self.analyte_reference is not None:
            parts.append(f"{self.analyte_reference}@")

        # Ion type
        parts.append(str(self.ion_type))

        # Neutral losses
        for loss in self.neutral_losses:
            parts.append(str(loss))

        # Isotopes
        for iso in self.isotopes:
            if iso.count != 0:
                parts.append(str(iso))

        # Adducts - reconstruct full adduct string
        if self.adducts:
            adduct_str = "M" + "".join(str(a) for a in self.adducts)
            parts.append(f"[{adduct_str}]")

        # Charge state (only if > 1)
        if self.charge > 1:
            parts.append(f"^{self.charge}")

        # Mass error
        if self.mass_error:
            parts.append(f"/{self.mass_error}")

        # Confidence
        if self.confidence is not None:
            parts.append(f"*{self.confidence:g}")

        return "".join(parts)

    @staticmethod
    def parse(annotation_str: str) -> "PafAnnotation":
        """Parse a single mzPAF annotation string into a FragmentAnnotation object"""
        parser = mzPAFParser()
        return parser.parse_single(annotation_str)

    def as_dict(self) -> dict:
        """Convert the annotation to a dictionary representation"""
        ion_dict = {}
        ion_dict["ion_type"] = type(self.ion_type).__name__
        ion_dict.update(asdict(self.ion_type))
        return {
            "ion": str(self.ion_type),
            "analyte_reference": self.analyte_reference,
            "is_auxiliary": self.is_auxiliary,
            "neutral_losses": [str(nl) for nl in self.neutral_losses],
            "isotopes": [str(iso) for iso in self.isotopes],
            "adducts": [str(ad) for ad in self.adducts],
            "charge": self.charge,
            "mass_error": str(self.mass_error) if self.mass_error else None,
            "confidence": self.confidence,
        }

    def __str__(self) -> str:
        return self.serialize()

    def __repr__(self) -> str:
        return f"PafAnnotation({self.as_dict()})"


class mzPAFParser:
    """Parser for mzPAF annotation strings following PSI specification"""

    # Regex components for better readability
    _AUXILIARY = r"(?P<is_auxiliary>&)?"
    _ANALYTE_REF = r"(?:(?P<analyte_reference>\d+)@)?"

    # Ion type patterns
    _PEPTIDE_SERIES = r"(?:(?P<series>(?:da|db|wa|wb)|[axbyczdwv]\.?)(?P<ordinal>\d+)(?:\{(?P<sequence_ordinal>.+)\})?)"
    _INTERNAL = r"(?P<series_internal>m(?P<internal_start>\d+):(?P<internal_end>\d+)(?:\{(?P<sequence_internal>.+)\})?)"
    _PRECURSOR = r"(?P<precursor>p)"
    _IMMONIUM = r"(?:I(?P<immonium>[A-Z])(?:\[(?P<immonium_modification>(?:[^\]]+))\])?)"
    _REFERENCE = r"(?P<reference>r(?:(?:\[(?P<reference_label>[^\]]+)\])))"
    _FORMULA = r"(?:f\{(?P<formula>[A-Za-z0-9\[\]]+)\})"
    _NAMED = r"(?:_\{(?P<named_compound>[^\{\}\s,/]+)\})"
    _SMILES = r"(?:s\{(?P<smiles>[^\}]+)\})"
    _UNKNOWN = r"(?:(?P<unannotated>\?)(?P<unannotated_label>\d+)?)"

    # Combine all ion types
    _ION_TYPES = f"(?:{_PEPTIDE_SERIES}|{_INTERNAL}|{_PRECURSOR}|{_IMMONIUM}|{_REFERENCE}|{_FORMULA}|{_NAMED}|{_SMILES}|{_UNKNOWN})"

    # Modifiers
    _NEUTRAL_LOSSES = r"(?P<neutral_losses>(?:[+-](?:\d+(?:\.\d+)?|\d*(?:(?:(?:\[[0-9]+[A-Z][A-Za-z0-9]*\])|(?:[A-Z][A-Za-z0-9]*))+)|(?:\d*\[(?:(?:[A-Za-z0-9:\.]+)(?:\[(?:[A-Za-z0-9\.:\-]+)\])?)\])))+)?"
    _ISOTOPE = r"(?P<isotope>(?:(?:[+-]\d*)i(?:(?:\d+)?(?:[A-Z][a-z]*)?|A)?)+)?"
    _ADDUCTS = r"(?:\[(?P<adducts>M(?:[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?"
    _CHARGE = r"(?:\^(?P<charge>[+-]?\d+))?"
    _MASS_ERROR = r"(?:/(?P<mass_error>-?\d+(?:\.\d+)?)(?P<mass_error_unit>ppm)?)?"
    _CONFIDENCE = r"(?:\*(?P<confidence>\d*(?:\.\d+)?))?"

    # Full pattern
    PATTERN = re.compile(
        f"^{_AUXILIARY}{_ANALYTE_REF}{_ION_TYPES}{_NEUTRAL_LOSSES}{_ISOTOPE}{_ADDUCTS}{_CHARGE}{_MASS_ERROR}{_CONFIDENCE}$"
    )

    def parse_single(self, annotation_str: str) -> PafAnnotation:
        """Parse a single annotation string"""
        match = self.PATTERN.match(annotation_str)
        if not match:
            raise ValueError(f"Invalid mzPAF annotation: '{annotation_str}'")

        groups = match.groupdict()

        return PafAnnotation(
            ion_type=self._parse_ion_type(groups),
            analyte_reference=self._parse_int(groups, "analyte_reference"),
            is_auxiliary=bool(groups.get("is_auxiliary")),
            neutral_losses=self._parse_neutral_losses(groups.get("neutral_losses")),
            isotopes=self._parse_isotopes(groups.get("isotope")),
            adducts=self._parse_adducts(groups.get("adducts")),
            charge=self._parse_int(groups, "charge") or 1,
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
        """Parse isotope notation into tuple of specifications"""
        if not isotope_str:
            return ()

        isotopes: list[IsotopeSpecification] = []
        # Pattern: [+-]?digits?i element?|A?
        pattern = r"([+-]?)(\d*)i((?:\d+)?(?:[A-Z][a-z]*)?|A)?"

        for match in re.finditer(pattern, isotope_str):
            sign_str, count_str, element_or_avg = match.groups()

            # Calculate signed count
            sign = -1 if sign_str == "-" else 1
            count = (int(count_str) if count_str else 1) * sign

            # Check if averaged isotopomer
            if element_or_avg == "A":
                isotopes.append(IsotopeSpecification(count=count, is_average=True))
            else:
                element = element_or_avg if element_or_avg else None
                isotopes.append(IsotopeSpecification(count=count, element=element))

        return tuple(isotopes)

    def _parse_neutral_losses(self, losses_str: str | None) -> tuple[NeutralLoss, ...]:
        """Parse neutral losses/gains into tuple of NeutralLoss objects"""
        if not losses_str:
            return ()

        # Pattern matches: [+-] followed by number, formula, or named group
        pattern = r"[+-](?:\d+(?:\.\d+)?(?!\[)|\d*(?:(?:\[[0-9]+[A-Z][A-Za-z0-9]*\])|(?:[A-Z][A-Za-z0-9]*))+|\d*\[(?:[A-Za-z0-9:\.]+)(?:\[[A-Za-z0-9\.:\-]+\])?\])"
        loss_strings = re.findall(pattern, losses_str)

        losses: list[NeutralLoss] = []
        for loss_str in loss_strings:
            sign = loss_str[0]
            sign_mult: int
            if sign == "+":
                sign_mult = 1
            elif sign == "-":
                sign_mult = -1
            else:
                raise ValueError(f"Invalid sign in neutral loss: '{loss_str}'")
            content = loss_str[1:]  # Remove sign

            # Try to parse as mass (decimal number)
            if re.match(r"^\d+(?:\.\d+)?$", content):
                count = 1 * sign_mult
                losses.append(NeutralLoss(count=count, base_mass=float(content)))

            # Parse as reference group [Name] or COUNT[Name]
            elif "[" in content:  # Changed from content.startswith('[')
                # Extract count and reference name
                match = re.match(r"^(\d*)\[([^\]]+)\]$", content)
                if match:
                    count_str, ref_name = match.groups()
                    count = int(count_str) if count_str else 1
                    losses.append(NeutralLoss(count=count * sign_mult, base_reference=ref_name))

            # Parse as formula (with optional count prefix)
            else:
                # Extract count and formula
                match = re.match(r"^(\d*)([A-Z].*)$", content)
                if match:
                    count_str, formula = match.groups()
                    count = int(count_str) if count_str else 1
                    losses.append(NeutralLoss(count=count * sign_mult, base_formula=formula))

        return tuple(losses)

    def _parse_adducts(self, adduct_str: str | None) -> tuple[Adduct, ...]:
        """Parse adduct notation into tuple of Adduct objects

        Examples:
            "M+H" -> (Adduct(sign=1, count=1, formula="H"),)
            "M+2H+Na" -> (Adduct(sign=1, count=2, formula="H"), Adduct(sign=1, count=1, formula="Na"))
            "M+NH4" -> (Adduct(sign=1, count=1, formula="NH4"),)
        """
        if not adduct_str:
            return ()

        # Verify M prefix
        if not adduct_str.startswith("M"):
            raise ValueError(f"Adduct string must start with 'M': '{adduct_str}'")

        content = adduct_str[1:]  # Remove 'M'
        if not content:
            raise ValueError(f"Adduct string must have components after 'M': '{adduct_str}'")

        # Pattern: [+-] followed by optional count and formula
        # Matches: +2H, +Na, -NH4, +H, etc.
        pattern = r"([+-])(\d*)([A-Z][A-Za-z0-9]*)"

        adducts: list[Adduct] = []
        for match in re.finditer(pattern, content):
            sign_str, count_str, formula = match.groups()
            sign: int
            if sign_str == "+":
                sign = 1
            elif sign_str == "-":
                sign = -1
            else:
                raise ValueError(f"Invalid sign in adduct: '{match.group(0)}'")
            count = int(count_str) if count_str else 1
            adducts.append(Adduct(count=count * sign, base_formula=formula))

        if not adducts:
            raise ValueError(f"No adduct components found in '{adduct_str}'")

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

    def parse(self, annotation_str: str) -> list[PafAnnotation]:
        """Parse potentially multiple comma-separated annotations"""
        if not annotation_str:
            return []

        annotations: list[PafAnnotation] = []
        for part in annotation_str.split(","):
            part = part.strip()
            if part:
                annotations.append(self.parse_single(part))

        return annotations


MZ_PAF_PARSER = mzPAFParser()


def parse(annotation_str: str) -> list[PafAnnotation]:
    """Convenience function to parse mzPAF annotation string into list of PafAnnotation"""
    return MZ_PAF_PARSER.parse(annotation_str)


def parse_single(annotation_str: str) -> PafAnnotation:
    """Convenience function to parse single mzPAF annotation string into PafAnnotation"""
    return parse(annotation_str)[0]
