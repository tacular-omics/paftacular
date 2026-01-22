"""
mzPAF Parser Module - Enhanced with backbone cleavage notation for internal fragments
Parses and serializes mzPAF (Peak Annotation Format) strings for peptide mass spectra.
Based on mzPAF specification v1.0 with extended internal fragment support
"""

from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from functools import cached_property
from typing import Literal

from tacular import ELEMENT_LOOKUP, FRAGMENT_ION_LOOKUP, REFMOL_LOOKUP, ElementInfo, RefMolInfo

from .constants import AminoAcids, IonSeries
from .util import parse_formula


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


def composition_to_formula_string(
    comp: Counter[ElementInfo],
) -> str:
    """Convert composition to standard chemical formula string"""
    # assert all are positive counts or all negative counts
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


class Serializable(ABC):
    """Base class for serializable objects"""

    @abstractmethod
    def serialize(self) -> str:
        pass

    def __str__(self) -> str:
        return self.serialize()


class MassProvider(ABC):
    """Base class for objects that can provide mass"""

    @abstractmethod
    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass"""
        pass

    @property
    def monoisotopic_mass(self) -> float:
        return self.mass(monoisotopic=True)

    @property
    def average_mass(self) -> float:
        return self.mass(monoisotopic=False)


class CompositionProvider(ABC):
    """Base class for objects that can provide composition"""

    @property
    @abstractmethod
    def composition(self) -> Counter[ElementInfo]:
        """Get elemental composition"""
        pass

    @property
    def dict_composition(self) -> dict[str, int]:
        """Get composition as a dictionary with element symbols as keys"""
        return {str(elem): count for elem, count in self.composition.items()}

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass from composition"""
        m = 0.0
        for elem, count in self.composition.items():
            m += elem.get_mass(monoisotopic) * count
        return m


@dataclass(frozen=True)
class MassError(Serializable):
    """Represents mass error with value and unit"""

    value: float
    unit: Literal["da", "ppm"] = "da"  # "da" or "ppm"

    def serialize(self) -> str:
        if self.unit == "ppm":
            return f"{self.value:g}ppm"
        elif self.unit == "da":
            return f"{self.value:g}"
        else:
            raise ValueError(f"Unknown mass error unit: {self.unit}")


@dataclass(frozen=True)
class IsotopeSpecification(Serializable, CompositionProvider, MassProvider):
    """Represents isotope information"""

    count: int = 0  # number of isotopes above/below monoisotope
    element: str | None = None  # e.g., "13C", "15N"
    is_average: bool = False  # True for averaged isotopomers

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


@dataclass(frozen=True)
class PeptideIon(Serializable, CompositionProvider, MassProvider):
    """Represents a primary peptide fragment ion"""

    series: IonSeries
    position: int
    sequence: str | None = None  # ProForma sequence

    def mass(self, monoisotopic: bool = True) -> float:
        return FRAGMENT_ION_LOOKUP[self.series].get_mass(monoisotopic)

    @property
    def formula(self) -> str:
        formula = FRAGMENT_ION_LOOKUP[self.series].formula
        if formula is None:
            raise ValueError(f"Formula not available for ion series: {self.series}")
        return formula

    @property
    def composition(self) -> Counter[ElementInfo]:
        comp: Counter[ElementInfo] = FRAGMENT_ION_LOOKUP[self.series].composition
        if comp is None:
            raise ValueError(f"Composition not available for ion series: {self.series}")
        return comp

    def serialize(self) -> str:
        result = f"{self.series}{self.position}"
        if self.sequence:
            result += f"{{{self.sequence}}}"
        return result


@dataclass(frozen=True)
class InternalFragment(Serializable, CompositionProvider, MassProvider):
    """Represents an internal fragment ion with optional backbone cleavage specification"""

    start_position: int
    end_position: int
    sequence: str | None = None

    # Optional backbone cleavage types.
    # The mzPAF documentation specifies these using neutral loss for some reason...
    nterm_ion_type: IonSeries | None = None  # e.g., IonSeries.A, IonSeries.B, IonSeries.C
    cterm_ion_type: IonSeries | None = None  # e.g., IonSeries.X, IonSeries.Y, IonSeries.Z

    def serialize(self) -> str:
        # If using default yb cleavage, just use 'm'
        result = f"m{self.start_position}:{self.end_position}"

        if self.sequence:
            result += f"{{{self.sequence}}}"
        return result

    def mass(self, monoisotopic: bool = True) -> float:
        # Start with base mass of internal fragment (no cleavage)
        return FRAGMENT_ION_LOOKUP["by"].get_mass(monoisotopic)

    @property
    def formula(self) -> str:
        formula = FRAGMENT_ION_LOOKUP["by"].formula
        if formula is None:
            raise ValueError("Formula not available for internal fragment")
        return formula

    @property
    def composition(self) -> Counter[ElementInfo]:
        comp: Counter[ElementInfo] = FRAGMENT_ION_LOOKUP["by"].composition
        if comp is None:
            raise ValueError("Composition not available for internal fragment")
        return comp


@dataclass(frozen=True)
class ImmoniumIon(Serializable, CompositionProvider, MassProvider):
    """Represents an immonium ion"""

    amino_acid: AminoAcids
    modification: str | None = None

    def serialize(self) -> str:
        result = f"I{self.amino_acid}"
        if self.modification:
            result += f"[{self.modification}]"
        return result

    def mass(self, monoisotopic: bool = True) -> float:
        if self.modification is not None:
            raise NotImplementedError("Mass calculation for modified immonium ions is not implemented")

        return FRAGMENT_ION_LOOKUP[self.amino_acid].get_mass(monoisotopic) + FRAGMENT_ION_LOOKUP["by"].get_mass(
            monoisotopic
        )

    @property
    def formula(self) -> str:
        return composition_to_proforma_formula_string(self.composition)

    @property
    def composition(self) -> Counter[ElementInfo]:
        comp: Counter[ElementInfo] = (
            FRAGMENT_ION_LOOKUP[self.amino_acid].composition + FRAGMENT_ION_LOOKUP["by"].composition
        )
        if comp is None:
            raise ValueError(f"Composition not available for immonium ion of amino acid: {self.amino_acid}")
        return comp


@dataclass(frozen=True)
class ReferenceIon(Serializable, CompositionProvider, MassProvider):
    """Represents a reference ion"""

    name: str

    @property
    def reference(self) -> RefMolInfo:
        return REFMOL_LOOKUP[self.name]

    def mass(self, monoisotopic: bool = True) -> float:
        return self.reference.get_mass(monoisotopic)

    @property
    def formula(self) -> str | None:
        return self.reference.chemical_formula

    @property
    def composition(self) -> Counter[ElementInfo]:
        return self.reference.composition

    def serialize(self) -> str:
        return f"r[{self.name}]"


@dataclass(frozen=True)
class NamedCompound(Serializable, CompositionProvider, MassProvider):
    """
    Represents a named compound.

    Example: 0@_{Urocanic Acid}
    """

    name: str

    def mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError("Mass calculation for NamedCompound is not implemented")

    @property
    def composition(self) -> Counter[ElementInfo]:
        raise NotImplementedError("Composition calculation for NamedCompound is not implemented")

    def serialize(self) -> str:
        return f"_{{{self.name}}}"


@dataclass(frozen=True)
class ChemicalFormula(Serializable, CompositionProvider, MassProvider):
    """
    Represents a chemical formula

    Example:
        f{C13H9}/-0.55ppm
        f{C12H9N}/0.06ppm
        f{C13H9N}/-2.01ppm
        f{C13H10N}/-0.11ppm
        f{C13H11N}/-0.09ppm
        f{C13H12N}/0.26ppm
        f{C14H10N}/0.19ppm
        f{C14H11N}/0.45ppm
        f{C14H10NO}/0.03ppm
    """

    formula: str

    @property
    def proforma_formula(self) -> str:
        return self.formula

    @property
    def composition(self) -> Counter[ElementInfo]:
        return formula_to_composition(self.formula)

    def serialize(self) -> str:
        return f"f{{{self.formula}}}"


@dataclass(frozen=True)
class SMILESCompound(Serializable, CompositionProvider, MassProvider):
    """
    Represents a SMILES string

    Example:
        s{CN=C=O}[M+H]/-0.55ppm
        s{COc(c1)cccc1C#N}[M+H+Na]^2/1.29ppm
    """

    smiles: str

    def serialize(self) -> str:
        return f"s{{{self.smiles}}}"

    @cached_property
    def composition(self) -> Counter[ElementInfo]:
        try:
            import pysmiles
        except ImportError as e:
            raise ImportError("pysmiles is required for SMILES parsing. Install with: pip install pysmiles") from e

        try:
            mol = pysmiles.read_smiles(self.smiles, explicit_hydrogen=True)
        except Exception as e:
            raise ValueError(f"Invalid SMILES string '{self.smiles}': {e}") from e

        elem_counts: Counter[str] = Counter()
        for node_id in mol.nodes():
            elem = mol.nodes[node_id].get("element", "*")
            if elem == "*":
                raise ValueError(
                    f"Unknown element '*' in SMILES '{self.smiles}'. Ensure all atoms are properly specified."
                )
            elem_counts[elem] += 1

        return Counter({ELEMENT_LOOKUP[elem]: count for elem, count in elem_counts.items()})

    @property
    def proforma_formula(self) -> str:
        return composition_to_proforma_formula_string(self.composition)

    @property
    def formula(self) -> str:
        return f"+{self.proforma_formula}"


@dataclass(frozen=True)
class UnknownIon(Serializable, CompositionProvider, MassProvider):
    """Represents an unknown/unannotated ion"""

    label: int | None = None

    def mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError("Mass calculation for UnknownIon is not implemented")

    @property
    def composition(self) -> Counter[ElementInfo]:
        raise NotImplementedError("Composition calculation for UnknownIon is not implemented")

    def serialize(self) -> str:
        if self.label is not None:
            return f"?{self.label}"
        return "?"


@dataclass(frozen=True)
class PrecursorIon(Serializable, CompositionProvider, MassProvider):
    """Represents a precursor ion"""

    def serialize(self) -> str:
        return "p"

    def mass(self, monoisotopic: bool = True) -> float:
        return FRAGMENT_ION_LOOKUP["p"].get_mass(monoisotopic)

    @property
    def formula(self) -> str | None:
        return FRAGMENT_ION_LOOKUP["p"].formula

    @property
    def composition(self) -> Counter[ElementInfo]:
        return FRAGMENT_ION_LOOKUP["p"].composition


class ScalableComposition(CompositionProvider):
    """Mixin for compositions that scale by count and sign"""

    sign: int
    count: int

    @property
    @abstractmethod
    def _single_composition(self) -> Counter[ElementInfo]:
        """Get composition for single instance (before scaling)"""
        pass

    @property
    def composition(self) -> Counter[ElementInfo]:
        """Get scaled composition"""
        return Counter({elem: count * self.count * self.sign for elem, count in self._single_composition.items()})

    @property
    def _sign_prefix(self) -> str:
        """Get sign prefix for serialization"""
        sign_str = "+" if self.sign > 0 else "-"
        count_str = "" if self.count == 1 else str(self.count)
        return f"{sign_str}{count_str}"


@dataclass(frozen=True, slots=True)
class NeutralLoss(
    Serializable,
    ScalableComposition,
    MassProvider,
):
    """Represents a neutral loss or gain"""

    sign: int  # 1 for gain (+), -1 for loss (-)
    count: int = 1  # Number of times this loss occurs
    _formula: str | None = None  # e.g., "H2O", "NH3"
    _mass: float | None = None  # e.g., 17.03 for direct mass specification
    _reference: str | None = None  # e.g., "Phospho", "iTRAQ115" (without brackets)

    def __post_init__(self):
        """Validate that exactly one of formula/mass/reference is set"""
        set_count = sum([self._formula is not None, self._mass is not None, self._reference is not None])
        if set_count != 1:
            raise ValueError("Exactly one of formula, mass, or reference must be set")

    @property
    def reference(self) -> RefMolInfo | str | None:
        if self._reference is None:
            return None
        try:
            val = REFMOL_LOOKUP[self._reference]
            return val
        except KeyError:
            return self._reference

    @property
    def loss_type(self) -> Literal["mass", "formula", "reference"]:
        if self._mass is not None:
            return "mass"
        elif self._formula is not None:
            return "formula"
        elif self._reference is not None:
            return "reference"
        else:
            raise ValueError("Invalid NeutralLoss state")

    @property
    def _single_composition(self) -> Counter[ElementInfo]:
        match self.loss_type:
            case "formula":
                if self._formula is None:  # This shouldn't happen given __post_init__
                    raise RuntimeError("Invalid state: formula is None")
                return formula_to_composition(self._formula)
            case "reference":
                refmol = self.reference
                if not isinstance(refmol, RefMolInfo):
                    raise ValueError(
                        f"Unknown reference molecule '{self._reference}'. Check that it exists in REFMOL_LOOKUP."
                    )
                return refmol.composition
            case "mass":
                raise ValueError(
                    f"Cannot calculate composition for mass-based loss ({self._mass} Da). "
                    f"Use a formula or reference instead."
                )

    @property
    def proforma_formula(self) -> str:
        return composition_to_proforma_formula_string(self.composition)

    @property
    def _single_formula(self) -> str:
        """Get formula for a single instance of the loss (without count/sign)"""
        match self.loss_type:
            case "formula":
                if self._formula is None:
                    raise RuntimeError("Formula is None for formula-based loss")
                return self._formula
            case "reference":
                refmol: RefMolInfo | str | None = self.reference
                if isinstance(refmol, RefMolInfo):
                    return refmol.chemical_formula
                else:
                    raise ValueError(
                        f"Cannot get formula for unknown reference molecule '{refmol}' of type: {type(refmol)}"
                    )
            case "mass":
                raise ValueError(f"Cannot get formula for mass-based loss: {self._mass}")
            case _:
                raise ValueError(f"Invalid loss_type: {self.loss_type}")

    @property
    def formula(self) -> str:
        single_formula = self._single_formula
        return f"{self._sign_prefix}{single_formula}"

    def _mass_single(self, monoisotopic: bool = True) -> float:
        match self.loss_type:
            case "mass":
                if self._mass is None:
                    raise RuntimeError("Mass is None for mass-based loss")
                return self._mass
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
        return single_mass * self.count * self.sign

    def serialize(
        self, loss_type: Literal["mass", "formula", "reference"] | None = None, monoisotopic: bool = True
    ) -> str:
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
                if self._reference is not None:
                    ref_name = self._reference
                    return f"{self._sign_prefix}[{ref_name}]"
                else:
                    raise ValueError("Cannot serialize reference: reference name is undefined")

        raise ValueError("Invalid loss_type for serialization")


@dataclass(frozen=True, slots=True)
class Adduct(Serializable, ScalableComposition, MassProvider):
    sign: int
    count: int = 1
    _formula: str = ""

    def __post_init__(self):
        """Validate adduct"""
        if self.sign not in (-1, 1):
            raise ValueError(f"Sign must be 1 or -1, got {self.sign}")
        if self.count < 1:
            raise ValueError(f"Count must be >= 1, got {self.count}")
        if not self._formula:
            raise ValueError("Formula cannot be empty")

    @property
    def _single_composition(self) -> Counter[ElementInfo]:
        return formula_to_composition(self._formula)  # Use helper!

    @property
    def formula(self) -> str:
        return f"{self._sign_prefix}{self._formula}"

    @property
    def proforma_formula(self) -> str:
        return composition_to_proforma_formula_string(self._single_composition)  # Use helper!

    def serialize(self) -> str:
        return f"{self._sign_prefix}{self._formula}"
