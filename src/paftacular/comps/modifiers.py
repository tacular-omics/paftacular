"""Modifier components for mzPAF annotations"""

import re
from collections import Counter
from dataclasses import dataclass
from typing import ClassVar, Literal

from tacular import ELEMENT_LOOKUP, REFMOL_LOOKUP, ElementInfo, RefMolInfo

from paftacular.constants import ADDUCT_REGEX_PATTERN, ISOTOPE_REGEX_PATTERN

from ..constants import MAX_CACHE_SIZE
from .base import CompositionProvider, MassProvider, ScalableComposition, Serializable
from .util import composition_to_proforma_formula_string, formula_to_composition


@dataclass(frozen=True, slots=True)
class MassError(Serializable):
    """Represents mass error with value and unit"""

    value: float
    unit: Literal["da", "ppm"] = "da"

    def serialize(self) -> str:
        if self.unit == "ppm":
            return f"{self.value:g}ppm"
        elif self.unit == "da":
            return f"{self.value:g}"
        else:
            raise ValueError(f"Unknown mass error unit: {self.unit}")

    @staticmethod
    def parse(s: str) -> "MassError":
        """Parse mass error string like '0.55ppm' or '0.06'"""
        s = s.strip()
        if s.endswith("ppm"):
            return MassError(value=float(s[:-3]), unit="ppm")
        else:
            return MassError(value=float(s), unit="da")


@dataclass(frozen=True, slots=True)
class IsotopeSpecification(Serializable, CompositionProvider, MassProvider):
    """Represents isotope information"""

    count: int = 0  # number of isotopes above/below monoisotope
    element: str | None = None  # e.g., "13C", "15N"
    is_average: bool = False  # True for averaged isotopomers

    _cache: ClassVar[dict[tuple, "IsotopeSpecification"]] = {}

    def __new__(cls, count: int = 0, element: str | None = None, is_average: bool = False):
        """Create or retrieve cached instance"""
        key = (count, element, is_average)
        if key not in cls._cache:
            # Evict oldest entry if cache is full
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

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
        match = re.match(ISOTOPE_REGEX_PATTERN, s)
        if not match:
            raise ValueError(f"Invalid isotope specification: '{s}'")

        sign_str, count_str, element_or_avg = match.groups()
        sign = -1 if sign_str == "-" else 1
        count = (int(count_str) if count_str else 1) * sign

        if element_or_avg == "A":
            return IsotopeSpecification(count=count, is_average=True)
        elif element_or_avg:
            return IsotopeSpecification(count=count, element=element_or_avg)
        else:
            return IsotopeSpecification(count=count)

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass contribution of isotope specification"""
        if monoisotopic is False:
            raise ValueError("Cannot calculate mass shift for average isotopomer specification")
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
            raise ValueError("Cannot calculate composition for average isotopomer specification")

        if self.element is None:
            raise ValueError("Cannot calculate composition for generic isotope specification without element")

        if self.element not in ELEMENT_LOOKUP:
            raise ValueError(f"Unknown element for isotope specification: {self.element}")

        elem_info: ElementInfo = ELEMENT_LOOKUP[self.element]
        # Get monoisotopic using the base element symbol (e.g., "C" from "13C")
        base_symbol = elem_info.symbol
        mono_info: ElementInfo = ELEMENT_LOOKUP.get_monoisotopic(base_symbol)
        comp: Counter[ElementInfo] = Counter()
        comp[elem_info] = self.count
        comp[mono_info] = -self.count
        return comp

    def as_dict(self) -> dict:
        """Convert isotope specification to dictionary representation"""
        try:
            _monoisotopic_mass = round(self.monoisotopic_mass, 5)
        except ValueError:
            _monoisotopic_mass = None

        try:
            _average_mass = round(self.average_mass, 5)
        except ValueError:
            _average_mass = None

        try:
            _dict_composition = self.dict_composition
        except ValueError:
            _dict_composition = None

        return {
            "count": self.count,
            "element": self.element,
            "is_average": self.is_average,
            "monoisotopic_mass": _monoisotopic_mass,
            "average_mass": _average_mass,
            "composition": _dict_composition,
        }


@dataclass(frozen=True, slots=True)
class NeutralLoss(
    Serializable,
    ScalableComposition,
    MassProvider,
):
    """Represents a neutral loss or gain"""

    count: int
    base_formula: str | None = None  # e.g., "H2O", "NH3"
    base_mass: float | None = None  # e.g., 17.03 for direct mass specification
    base_reference: str | None = None  # e.g., "Phospho", "iTRAQ115" (without brackets)

    _cache: ClassVar[dict[tuple, "NeutralLoss"]] = {}

    def __new__(cls, count: int, base_formula: str | None = None, base_mass: float | None = None, base_reference: str | None = None):
        """Create or retrieve cached instance"""
        key = (count, base_formula, base_mass, base_reference)
        if key not in cls._cache:
            # Evict oldest entry if cache is full
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

    def __post_init__(self):
        """Validate that exactly one of formula/mass/reference is set"""
        set_count = sum([self.base_formula is not None, self.base_mass is not None, self.base_reference is not None])
        if set_count != 1:
            raise ValueError("Exactly one of formula, mass, or reference must be set")

    @property
    def reference(self) -> RefMolInfo | str | None:
        if self.base_reference is None:
            return None
        try:
            val = REFMOL_LOOKUP[self.base_reference]
            return val
        except KeyError:
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
            raise ValueError("Invalid NeutralLoss state")

    @property
    def _single_composition(self) -> Counter[ElementInfo]:
        match self.loss_type:
            case "formula":
                if self.base_formula is None:  # This shouldn't happen given __post_init__
                    raise RuntimeError("Invalid state: formula is None")
                return formula_to_composition(self.base_formula)
            case "reference":
                refmol = self.reference
                if not isinstance(refmol, RefMolInfo):
                    raise ValueError(f"Unknown reference molecule '{self.base_reference}'. Check that it exists in REFMOL_LOOKUP.")
                return refmol.composition
            case "mass":
                raise ValueError(f"Cannot calculate composition for mass-based loss ({self.base_mass} Da). Use a formula or reference instead.")

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
                refmol: RefMolInfo | str | None = self.reference
                if isinstance(refmol, RefMolInfo):
                    return refmol.chemical_formula
                else:
                    raise ValueError(f"Cannot get formula for unknown reference molecule '{refmol}' of type: {type(refmol)}")
            case "mass":
                raise ValueError(f"Cannot get formula for mass-based loss: {self.base_mass}")
            case _:
                raise ValueError(f"Invalid loss_type: {self.loss_type}")

    @property
    def formula(self) -> str:
        single_formula = self._single_formula
        return f"{self._sign_prefix}{single_formula}"

    def _mass_single(self, monoisotopic: bool = True) -> float:
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
                    m += elem.get_mass(monoisotopic) * count
                return m
            case "reference":
                refmol: RefMolInfo | str | None = self.reference
                if isinstance(refmol, str) or refmol is None:
                    raise ValueError(f"Cannot get mass for unknown reference molecule '{refmol}'")
                return refmol.get_mass(monoisotopic)

    def mass(self, monoisotopic: bool = True) -> float:
        single_mass: float = self._mass_single(monoisotopic)
        return single_mass * self.count

    def serialize(self, loss_type: Literal["mass", "formula", "reference"] | None = None, monoisotopic: bool = True) -> str:
        if loss_type is None:
            loss_type = self.loss_type

        match loss_type:
            case "mass":
                mass = self.mass(monoisotopic=monoisotopic)
                return f"{mass:+.5f}"
            case "formula":
                formula = self.formula
                return f"{formula}"
            case "reference":
                if self.base_reference is not None:
                    ref_name = self.base_reference
                    return f"{self._sign_prefix}[{ref_name}]"
                else:
                    raise ValueError("Cannot serialize reference: reference name is undefined")

        raise ValueError("Invalid loss_type for serialization")

    def as_dict(self) -> dict:
        """Convert the neutral loss to a dictionary representation"""
        try:
            _monoisotopic_mass = round(self.monoisotopic_mass, 5)
        except ValueError:
            _monoisotopic_mass = None

        try:
            _average_mass = round(self.average_mass, 5)
        except ValueError:
            _average_mass = None

        try:
            _formula = self.formula
        except ValueError:
            _formula = None

        try:
            _dict_composition = self.dict_composition
        except ValueError:
            _dict_composition = None

        return {
            "count": self.count,
            "base_formula": self.base_formula,
            "base_mass": self.base_mass,
            "base_reference": self.base_reference,
            "monoisotopic_mass": _monoisotopic_mass,
            "average_mass": _average_mass,
            "composition": _dict_composition,
            "formula": _formula,
        }

    @staticmethod
    def parse(loss_str: str) -> "NeutralLoss":
        """Parse a neutral loss string into a NeutralLoss object"""
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
            return NeutralLoss(count=count, base_mass=float(content))

        # Parse as reference group [Name] or COUNT[Name]
        elif "[" in content:  # Changed from content.startswith('[')
            # Extract count and reference name
            match = re.match(r"^(\d*)\[([^\]]+)\]$", content)
            if match:
                count_str, ref_name = match.groups()
                count = int(count_str) if count_str else 1
                return NeutralLoss(count=count * sign_mult, base_reference=ref_name)

        # Parse as formula (with optional count prefix)
        else:
            # Extract count and formula
            match = re.match(r"^(\d*)([A-Z].*)$", content)
            if match:
                count_str, formula = match.groups()
                count = int(count_str) if count_str else 1
                return NeutralLoss(count=count * sign_mult, base_formula=formula)

        raise ValueError(f"Could not parse neutral loss: '{loss_str}'")


@dataclass(frozen=True, slots=True)
class Adduct(Serializable, ScalableComposition, MassProvider):
    count: int
    base_formula: str

    _cache: ClassVar[dict[tuple, "Adduct"]] = {}

    def __new__(cls, count: int, base_formula: str):
        """Create or retrieve cached instance"""
        key = (count, base_formula)
        if key not in cls._cache:
            # Evict oldest entry if cache is full
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

    def __post_init__(self):
        """Validate adduct"""
        if self.count == 0:
            raise ValueError(f"Count must be non-zero, got {self.count}")
        if not self.base_formula:
            raise ValueError("Formula cannot be empty")

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
        match = re.match(ADDUCT_REGEX_PATTERN, s)
        if not match:
            raise ValueError(f"Invalid adduct: '{s}'")

        sign_str, count_str, formula = match.groups()
        sign = 1 if sign_str == "+" else -1
        count = (int(count_str) if count_str else 1) * sign
        return Adduct(count=count, base_formula=formula)

    def as_dict(self) -> dict:
        """Convert the adduct to a dictionary representation"""

        try:
            _monoisotopic_mass = round(self.monoisotopic_mass, 5)
        except ValueError:
            _monoisotopic_mass = None

        try:
            _average_mass = round(self.average_mass, 5)
        except ValueError:
            _average_mass = None

        try:
            _formula = self.formula
        except ValueError:
            _formula = None

        try:
            _dict_composition = self.dict_composition
        except ValueError:
            _dict_composition = None

        return {
            "count": self.count,
            "base_formula": self.base_formula,
            "monoisotopic_mass": _monoisotopic_mass,
            "average_mass": _average_mass,
            "composition": _dict_composition,
            "formula": _formula,
        }
