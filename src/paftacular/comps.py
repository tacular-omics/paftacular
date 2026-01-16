"""
mzPAF Parser Module - Enhanced with backbone cleavage notation for internal fragments
Parses and serializes mzPAF (Peak Annotation Format) strings for peptide mass spectra.
Based on mzPAF specification v1.0 with extended internal fragment support
"""

from collections import Counter
from dataclasses import dataclass
from typing import Literal

from .constants import AminoAcids, IonSeries
from .elements import ELEMENT_LOOKUP, ElementInfo
from .refmols import REFMOL_LOOKUP, RefMolInfo
from .util import parse_formula


@dataclass(frozen=True)
class MassError:
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

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True)
class IsotopeSpecification:
    """Represents isotope information"""

    count: int = 0  # number of isotopes above/below monoisotope
    element: str | None = None  # e.g., "13C", "15N"
    is_average: bool = False  # True for averaged isotopomers

    def serialize(self) -> str:
        if self.count == 0:
            return ""

        # Explicit sign: '+' for positive, '-' for negative (no sign for zero handled above)
        sign = "+" if self.count > 0 else "-"
        count_str = "" if abs(self.count) == 1 else str(abs(self.count))

        if self.is_average is True:
            return f"{sign}{count_str}iA"
        elif self.element is not None:
            return f"{sign}{count_str}i{self.element}"
        else:
            return f"{sign}{count_str}i"

    def mass_shift(self) -> float:
        """Calculate mass contribution of isotope specification"""
        if self.count == 0:
            return 0.0

        if self.is_average:
            # Average isotopomer mass shift is not well-defined; raise error
            raise ValueError("Cannot calculate mass for average isotopomer specification")

        if self.element is None:
            # Generic isotope shift not defined
            raise ValueError("Cannot calculate mass for generic isotope specification without element")

        if self.element not in ELEMENT_LOOKUP:
            raise ValueError(f"Unknown element for isotope specification: {self.element}")

        elem_info: ElementInfo = ELEMENT_LOOKUP[self.element]
        mono_info: ElementInfo = ELEMENT_LOOKUP.get_monoisotopic(self.element)
        mass_diff = elem_info.mass - mono_info.mass
        return mass_diff * self.count

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True)
class PeptideIon:
    """Represents a primary peptide fragment ion"""

    series: IonSeries
    position: int
    sequence: str | None = None  # ProForma sequence

    def serialize(self) -> str:
        result = f"{self.series}{self.position}"
        if self.sequence:
            result += f"{{{self.sequence}}}"
        return result

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class InternalFragment:
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

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class ImmoniumIon:
    """Represents an immonium ion"""

    amino_acid: AminoAcids
    modification: str | None = None

    def serialize(self) -> str:
        result = f"I{self.amino_acid}"
        if self.modification:
            result += f"[{self.modification}]"
        return result

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class ReferenceIon:
    """Represents a reference ion"""

    name: str

    @property
    def reference(self) -> RefMolInfo | None:
        try:
            val = REFMOL_LOOKUP[self.name]
            return val
        except KeyError:
            return None

    def serialize(self) -> str:
        return f"r[{self.name}]"

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class NamedCompound:
    """
    Represents a named compound.

    Example: 0@_{Urocanic Acid}
    """

    name: str

    def serialize(self) -> str:
        return f"_{{{self.name}}}"

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class ChemicalFormula:
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
    def composition(self) -> Counter[ElementInfo]:
        elem_counts: Counter[str] = parse_formula(self.formula)
        comp: Counter[ElementInfo] = Counter()
        for elem, count in elem_counts.items():
            comp[ELEMENT_LOOKUP[elem]] = count
        return comp

    def mass(self, monoisotopic: bool = True) -> float:
        m = 0.0
        for elem, count in self.composition.items():
            m += elem.get_mass(monoisotopic) * count
        return m

    @property
    def monoisotopic_mass(self) -> float:
        return self.mass(monoisotopic=True)

    @property
    def average_mass(self) -> float:
        return self.mass(monoisotopic=False)

    def serialize(self) -> str:
        return f"f{{{self.formula}}}"

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class SMILESCompound:
    """
    Represents a SMILES string

    Example:
        s{CN=C=O}[M+H]/-0.55ppm
        s{COc(c1)cccc1C#N}[M+H+Na]^2/1.29ppm
    """

    smiles: str

    def serialize(self) -> str:
        return f"s{{{self.smiles}}}"

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class UnknownIon:
    """Represents an unknown/unannotated ion"""

    label: int | None = None

    def serialize(self) -> str:
        if self.label is not None:
            return f"?{self.label}"
        return "?"

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True)
class PrecursorIon:
    """Represents a precursor ion"""

    def serialize(self) -> str:
        return "p"

    def __str__(self):
        return self.serialize()


@dataclass(frozen=True, slots=True)
class NeutralLoss:
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
    def composition(self) -> Counter[ElementInfo] | None:
        match self.loss_type:
            case "formula":
                if self._formula is None:
                    raise RuntimeError("Formula is None for formula-based loss")
                elem_counts: Counter[str] = parse_formula(self._formula)
                comp: Counter[ElementInfo] = Counter()
                for elem, count in elem_counts.items():
                    comp[ELEMENT_LOOKUP[elem]] = count * self.count * self.sign
                return comp
            case "reference":
                refmol: RefMolInfo | str | None = self.reference
                if isinstance(refmol, RefMolInfo):
                    comp = Counter()
                    for elem, count in refmol.composition.items():
                        comp[elem] = count * self.count * self.sign
                    return comp
                else:
                    raise ValueError(f"Cannot get composition for unknown reference molecule '{refmol}'")
            case "mass":
                raise ValueError("Mass-based losses do not have a defined composition")

    @property
    def formula(self) -> str | None:
        comp: Counter[ElementInfo] | None = self.composition
        if comp is None:
            return None
        s = ""
        for elem in sorted(comp.keys(), key=lambda e: e.symbol):
            count = comp[elem]
            if count == 0:
                continue
            s += f"{elem.symbol}{count if count != 1 else ''}"
        return s if s else None

    def mass(self, monoisotopic: bool = True) -> float | None:
        match self.loss_type:
            case "mass":
                if self._mass is None:
                    raise RuntimeError("Mass is None for mass-based loss")
                return self._mass * self.count * self.sign
            case "formula":
                comp: Counter[ElementInfo] | None = self.composition
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
                return refmol.get_mass(monoisotopic) * self.count * self.sign

    def __str__(self) -> str:
        """Serialize back to mzPAF format"""
        return self.serialize()

    def serialize(self, loss_type: Literal["mass", "formula", "reference"] | None = None) -> str:
        sign_str = "+" if self.sign > 0 else "-"
        count_str = "" if self.count == 1 else str(self.count)

        if loss_type is None:
            loss_type = self.loss_type

        match loss_type:
            case "mass":
                mass = self.mass()
                if mass is None:
                    raise ValueError("Cannot serialize mass: mass is undefined")
                return f"{sign_str}{count_str}{mass}"
            case "formula":
                formula = self.formula
                if formula is None:
                    raise ValueError("Cannot serialize formula: composition is undefined")
                return f"{sign_str}{count_str}{formula}"
            case "reference":
                ref = self.reference
                if isinstance(ref, RefMolInfo):
                    ref_name = ref.name
                elif isinstance(ref, str):
                    ref_name = ref
                else:
                    raise ValueError("Cannot serialize reference: reference is undefined")
                return f"{sign_str}{count_str}[{ref_name}]"

        raise ValueError("Invalid loss_type for serialization")


@dataclass(frozen=True, slots=True)
class Adduct:
    """Represents a single adduct component (e.g., H, Na, NH4)"""

    sign: int  # 1 for addition (+), -1 for removal (-)
    count: int = 1  # Number of this adduct
    formula: str = ""  # Chemical formula or element symbol (H, Na, NH4, etc.)

    def __post_init__(self):
        """Validate adduct"""
        if self.sign not in (-1, 1):
            raise ValueError(f"Sign must be 1 or -1, got {self.sign}")
        if self.count < 1:
            raise ValueError(f"Count must be >= 1, got {self.count}")
        if not self.formula:
            raise ValueError("Formula cannot be empty")

    @property
    def composition(self) -> Counter[ElementInfo]:
        """Get elemental composition of this adduct"""
        elem_counts: Counter[str] = parse_formula(self.formula)
        comp: Counter[ElementInfo] = Counter()
        for elem, count in elem_counts.items():
            comp[ELEMENT_LOOKUP[elem]] = count * self.count * self.sign
        return comp

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass contribution of this adduct"""
        comp = self.composition
        m = 0.0
        for elem, count in comp.items():
            m += elem.get_mass(monoisotopic) * count
        return m

    def __str__(self) -> str:
        """Serialize to mzPAF format component (without M prefix)"""
        sign_str = "+" if self.sign > 0 else "-"
        count_str = "" if self.count == 1 else str(self.count)
        return f"{sign_str}{count_str}{self.formula}"
