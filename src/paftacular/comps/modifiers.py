"""Modifier components for mzPAF annotations"""

import re
from collections import Counter
from dataclasses import KW_ONLY, dataclass
from typing import Literal

from tacular import ELEMENT_LOOKUP, ElementInfo, RefMolInfo

from ..constants import _ATOM_TOKEN, ADDUCT_REGEX_PATTERN, ISOTOPE_REGEX_PATTERN
from ..errors import PaftacularError, PafUnknownReferenceError
from ..util import format_number, validate_number
from .base import CompositionProvider, MassProvider, ScalableComposition, Serializable
from .util import composition_to_proforma_formula_string, formula_to_composition, lookup_reference

_ISOTOPE_ELEMENT = re.compile(r"\d+[A-Z][a-z]?")
_MASS_CONTENT = re.compile(r"\d+(?:\.\d+)?")
_FORMULA_CONTENT = re.compile(rf"(\d*)({_ATOM_TOKEN}+)")
_REFERENCE_CONTENT = re.compile(r"(\d*)\[([^\]]+)\]")


@dataclass(frozen=True, slots=True)
class MassError(Serializable):
    """Represents mass error with value and unit"""

    value: float
    _: KW_ONLY
    unit: Literal["da", "ppm"] = "da"

    def __post_init__(self):
        validate_number(self.value)
        if self.unit not in ("da", "ppm"):
            raise PaftacularError(f"Unknown mass error unit: {self.unit}")

    def serialize(self) -> str:
        if self.unit == "ppm":
            return f"{format_number(self.value)}ppm"
        return format_number(self.value)

    @staticmethod
    def parse(s: str) -> "MassError":
        """Parse mass error string like '0.55ppm' or '0.06'"""
        s = s.strip()
        text, unit = (s[:-3], "ppm") if s.endswith("ppm") else (s, "da")
        try:
            value = float(text)
        except ValueError:
            raise PaftacularError(f"Invalid mass error: {s!r}") from None
        return MassError(value, unit=unit)


@dataclass(frozen=True, slots=True)
class IsotopeSpecification(Serializable, CompositionProvider, MassProvider):
    """Represents isotope information"""

    count: int = 0  # number of isotopes above/below monoisotope
    _: KW_ONLY
    element: str | None = None  # e.g., "13C", "15N"
    is_average: bool = False  # True for averaged isotopomers

    def __post_init__(self):
        if type(self.count) is not int or type(self.is_average) is not bool:
            raise PaftacularError("Isotope count must be an integer and is_average must be a boolean")
        if self.element is not None and (not isinstance(self.element, str) or not _ISOTOPE_ELEMENT.fullmatch(self.element)):
            raise PaftacularError("An isotope element requires a nucleon count and element symbol")
        if self.is_average and self.element is not None:
            raise PaftacularError("Average isotopes cannot also specify an element")

    @property
    def _prefix(self) -> str:
        """Get prefix for serialization"""
        sign = "+" if self.count > 0 else "-"
        count_str = "" if abs(self.count) == 1 else str(abs(self.count))
        return f"{sign}{count_str}"

    def serialize(self) -> str:
        if self.count == 0:
            return ""

        if self.is_average is True:
            return f"{self._prefix}iA"
        elif self.element is not None:
            return f"{self._prefix}i{self.element}"
        else:
            return f"{self._prefix}i"

    @staticmethod
    def parse(s: str) -> "IsotopeSpecification":
        """Parse isotope string like '+i', '-2i13C', '+iA'"""
        s = s.strip()
        match = re.fullmatch(ISOTOPE_REGEX_PATTERN, s)
        if not match:
            raise PaftacularError(f"Invalid isotope specification: '{s}'")

        sign_str, count_str, element_or_avg = match.groups()
        sign = -1 if sign_str == "-" else 1
        count = (int(count_str) if count_str else 1) * sign

        if element_or_avg == "A":
            return IsotopeSpecification(count, is_average=True)
        elif element_or_avg:
            return IsotopeSpecification(count, element=element_or_avg)
        else:
            return IsotopeSpecification(count)

    def get_mass(self, *, monoisotopic: bool = True) -> float:
        """Calculate mass contribution of isotope specification"""
        if monoisotopic is False:
            raise PaftacularError("Cannot calculate mass shift for average isotopomer specification")

        if self.count == 0:
            return 0.0

        if self.is_average:
            raise PaftacularError("Cannot calculate mass shift for average isotopomer specification")

        if self.element is None:
            # Generic isotope (no element specified): mzPAF section 4.6 defines this as the
            # difference between 13C and 12C, regardless of which atom actually carries it.
            c13 = ELEMENT_LOOKUP["13C"].get_mass(monoisotopic=True)
            c12 = ELEMENT_LOOKUP.get_monoisotopic("C").get_mass(monoisotopic=True)
            return (c13 - c12) * self.count

        comp = self.composition
        m = 0.0
        for elem, count in comp.items():
            m += elem.get_mass(monoisotopic=True) * count
        return m

    @property
    def composition(self) -> Counter[ElementInfo]:
        # lose mono and gain isotope
        if self.count == 0:
            return Counter()

        if self.is_average:
            raise PaftacularError("Cannot calculate composition for average isotopomer specification")

        if self.element is None:
            # Generic isotope (no element specified): mzPAF section 4.6 defines this as a 13C
            # substitution (13C in place of 12C) regardless of which atom actually carries it, so
            # gain one 13C and lose one 12C per count. Mirrors the element-specified path below and
            # keeps composition consistent with mass() (which uses the same 13C-12C shift).
            comp = Counter()
            comp[ELEMENT_LOOKUP["13C"]] = self.count
            comp[ELEMENT_LOOKUP.get_monoisotopic("C")] = -self.count
            return comp

        if self.element not in ELEMENT_LOOKUP:
            raise PaftacularError(f"Unknown element for isotope specification: {self.element}")

        elem_info: ElementInfo = ELEMENT_LOOKUP[self.element]
        # Get monoisotopic using the base element symbol (e.g., "C" from "13C")
        base_symbol = elem_info.symbol
        mono_info: ElementInfo = ELEMENT_LOOKUP.get_monoisotopic(base_symbol)
        comp: Counter[ElementInfo] = Counter()
        comp[elem_info] += self.count
        comp[mono_info] -= self.count
        return comp


@dataclass(frozen=True, slots=True)
class NeutralLoss(
    Serializable,
    ScalableComposition,
    MassProvider,
):
    """Represents a neutral loss or gain"""

    count: int
    _: KW_ONLY
    base_formula: str | None = None  # e.g., "H2O", "NH3"
    base_mass: float | None = None  # e.g., 17.03 for direct mass specification
    base_reference: str | None = None  # e.g., "Phospho", "iTRAQ115" (without brackets)

    def __post_init__(self):
        if type(self.count) is not int or self.count == 0:
            raise PaftacularError("Neutral loss count must be a nonzero integer")
        if sum(value is not None for value in (self.base_formula, self.base_mass, self.base_reference)) != 1:
            raise PaftacularError("Exactly one of formula, mass, or reference must be set")
        for value in (self.base_formula, self.base_reference):
            if value is not None and (not isinstance(value, str) or not value):
                raise PaftacularError("Formula and reference must be nonempty strings")
        if self.base_mass is not None:
            validate_number(self.base_mass)

    @property
    def reference(self) -> RefMolInfo | str | None:
        if self.base_reference is None:
            return None
        try:
            return lookup_reference(self.base_reference)
        except PafUnknownReferenceError:
            return self.base_reference

    @property
    def loss_type(self) -> Literal["mass", "formula", "reference"]:
        if self.base_mass is not None:
            return "mass"
        elif self.base_formula is not None:
            return "formula"
        elif self.base_reference is not None:
            return "reference"
        else:
            raise PaftacularError("Invalid NeutralLoss state")

    @property
    def _single_composition(self) -> Counter[ElementInfo]:
        match self.loss_type:
            case "formula":
                if self.base_formula is None:  # This shouldn't happen given __post_init__
                    raise RuntimeError("Invalid state: formula is None")
                return formula_to_composition(self.base_formula)
            case "reference":
                return lookup_reference(str(self.base_reference)).composition
            case "mass":
                raise PaftacularError(f"Cannot calculate composition for mass-based loss ({self.base_mass} Da). Use a formula or reference instead.")

    @property
    def proforma_formula(self) -> str:
        return composition_to_proforma_formula_string(self.composition)

    @property
    def _single_formula(self) -> str:
        """Get formula for a single instance of the loss (without count/sign)"""
        match self.loss_type:
            case "formula":
                if self.base_formula is None:
                    raise RuntimeError("Formula is None for formula-based loss")
                return self.base_formula
            case "reference":
                return lookup_reference(str(self.base_reference)).formula
            case "mass":
                raise PaftacularError(f"Cannot get formula for mass-based loss: {self.base_mass}")
            case _:
                raise PaftacularError(f"Invalid loss_type: {self.loss_type}")

    @property
    def formula(self) -> str:
        single_formula = self._single_formula
        return f"{self._sign_prefix}{single_formula}"

    def _mass_single(self, *, monoisotopic: bool = True) -> float:
        match self.loss_type:
            case "mass":
                if self.base_mass is None:
                    raise RuntimeError("Mass is None for mass-based loss")
                return self.base_mass
            case "formula":
                comp: Counter[ElementInfo] = self._single_composition
                if comp is None:
                    raise RuntimeError("Composition is None for formula-based loss")
                m = 0
                for elem, count in comp.items():
                    m += elem.get_mass(monoisotopic=monoisotopic) * count
                return m
            case "reference":
                return lookup_reference(str(self.base_reference)).get_mass(monoisotopic=monoisotopic)

    def get_mass(self, *, monoisotopic: bool = True) -> float:
        return self._mass_single(monoisotopic=monoisotopic) * self.count

    def serialize(self, *, loss_type: Literal["mass", "formula", "reference"] | None = None, monoisotopic: bool = True) -> str:
        if loss_type is None:
            loss_type = self.loss_type

        match loss_type:
            case "mass":
                mass = self.get_mass(monoisotopic=monoisotopic)
                return ("+" if mass >= 0 else "-") + format_number(abs(mass), minimum_places=5)
            case "formula":
                formula = self.formula
                return f"{formula}"
            case "reference":
                if self.base_reference is not None:
                    ref_name = self.base_reference
                    return f"{self._sign_prefix}[{ref_name}]"
                else:
                    raise PaftacularError("Cannot serialize reference: reference name is undefined")

        raise PaftacularError("Invalid loss_type for serialization")

    @staticmethod
    def parse(loss_str: str) -> "NeutralLoss":
        """Parse a neutral loss string into a NeutralLoss object"""
        loss_str = loss_str.strip()
        if not loss_str:
            raise PaftacularError("Empty neutral loss")
        sign = loss_str[0]
        sign_mult: int
        if sign == "+":
            sign_mult = 1
        elif sign == "-":
            sign_mult = -1
        else:
            raise PaftacularError(f"Invalid sign in neutral loss: '{loss_str}'")
        content = loss_str[1:]  # Remove sign

        # Try to parse as mass (decimal number)
        if _MASS_CONTENT.fullmatch(content):
            return NeutralLoss(sign_mult, base_mass=float(content))

        # Parse as a formula: one or more atoms, each either plain (e.g. "H2O") or an
        # isotope-labeled atom in brackets (e.g. "[18O1]"), optionally count-prefixed, and the
        # two forms may be mixed (e.g. "H2[18O1]"). Tried before the reference-group branch since
        # a reference name can never itself start with an atom/isotope-bracket token.
        match = _FORMULA_CONTENT.fullmatch(content)
        if match:
            count_str, formula = match.groups()
            count = int(count_str) if count_str else 1
            return NeutralLoss(count * sign_mult, base_formula=formula)

        # Parse as reference group [Name] or COUNT[Name]
        match = _REFERENCE_CONTENT.fullmatch(content)
        if match:
            count_str, ref_name = match.groups()
            count = int(count_str) if count_str else 1
            return NeutralLoss(count * sign_mult, base_reference=ref_name)

        raise PaftacularError(f"Could not parse neutral loss: '{loss_str}'")


@dataclass(frozen=True, slots=True)
class Adduct(Serializable, ScalableComposition, MassProvider):
    """Represents a charge-carrier adduct such as ``+Na`` or ``+2H``"""

    count: int
    base_formula: str

    def __post_init__(self):
        if type(self.count) is not int or self.count == 0:
            raise PaftacularError("Count must be a non-zero integer")
        if not isinstance(self.base_formula, str) or not self.base_formula:
            raise PaftacularError("Formula cannot be empty")

    @property
    def _single_composition(self) -> Counter[ElementInfo]:
        return formula_to_composition(self.base_formula)  # Use helper!

    @property
    def formula(self) -> str:
        return f"{self._sign_prefix}{self.base_formula}"

    @property
    def proforma_formula(self) -> str:
        return composition_to_proforma_formula_string(self._single_composition)  # Use helper!

    def serialize(self) -> str:
        return f"{self._sign_prefix}{self.base_formula}"

    @staticmethod
    def parse(s: str) -> "Adduct":
        """Parse a single adduct string like '+H', '+2Na', '-NH4'"""
        s = s.strip()
        match = re.fullmatch(ADDUCT_REGEX_PATTERN, s)
        if not match:
            raise PaftacularError(f"Invalid adduct: '{s}'")

        sign_str, count_str, formula = match.groups()
        sign = 1 if sign_str == "+" else -1
        count = (int(count_str) if count_str else 1) * sign
        return Adduct(count, formula)
