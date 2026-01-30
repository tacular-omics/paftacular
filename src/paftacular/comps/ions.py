"""Ion type definitions for mzPAF annotations"""

import re
from collections import Counter
from dataclasses import dataclass
from functools import cached_property
from typing import ClassVar

from tacular import ELEMENT_LOOKUP, FRAGMENT_ION_LOOKUP, REFMOL_LOOKUP, ElementInfo, RefMolInfo

from ..constants import MAX_CACHE_SIZE, AminoAcids, IonSeries
from .base import CompositionProvider, MassProvider, Serializable
from .util import composition_to_proforma_formula_string, formula_to_composition


@dataclass(frozen=True, slots=True)
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

    @staticmethod
    def parse(s: str) -> "PeptideIon":
        """Parse peptide ion string like 'b5', 'y10{PEPTIDE}'"""
        s = s.strip()
        match = re.match(r"((?:da|db|wa|wb)|[axbyczdwv]\.?)(\d+)(?:\{(.+)\})?", s)
        if not match:
            raise ValueError(f"Invalid peptide ion: '{s}'")

        series_str, position_str, sequence = match.groups()
        return PeptideIon(series=IonSeries(series_str), position=int(position_str), sequence=sequence)


@dataclass(frozen=True, slots=True)
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

    @staticmethod
    def parse(s: str) -> "InternalFragment":
        """Parse internal fragment string like 'm5:10', 'm5:10{PEPTIDE}'"""
        s = s.strip()
        match = re.match(r"m(\d+):(\d+)(?:\{(.+)\})?", s)
        if not match:
            raise ValueError(f"Invalid internal fragment: '{s}'")

        start_str, end_str, sequence = match.groups()
        return InternalFragment(start_position=int(start_str), end_position=int(end_str), sequence=sequence)

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


@dataclass(frozen=True, slots=True)
class ImmoniumIon(Serializable, CompositionProvider, MassProvider):
    """Represents an immonium ion"""

    amino_acid: AminoAcids
    modification: str | None = None

    _cache: ClassVar[dict[tuple, "ImmoniumIon"]] = {}

    def __new__(cls, amino_acid: AminoAcids, modification: str | None = None):
        """Create or retrieve cached instance"""
        key = (amino_acid, modification)
        if key not in cls._cache:
            # Evict oldest entry if cache is full
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

    def serialize(self) -> str:
        result = f"I{self.amino_acid}"
        if self.modification:
            result += f"[{self.modification}]"
        return result

    @staticmethod
    def parse(s: str) -> "ImmoniumIon":
        """Parse immonium ion string like 'IK', 'IM[Oxidation]'"""
        s = s.strip()
        match = re.match(r"I([A-Z])(?:\[([^\]]+)\])?", s)
        if not match:
            raise ValueError(f"Invalid immonium ion: '{s}'")

        aa_str, modification = match.groups()
        return ImmoniumIon(amino_acid=AminoAcids(aa_str), modification=modification)

    def mass(self, monoisotopic: bool = True) -> float:
        if self.modification is not None:
            raise NotImplementedError("Mass calculation for modified immonium ions is not implemented")

        return FRAGMENT_ION_LOOKUP[self.amino_acid].get_mass(monoisotopic) + FRAGMENT_ION_LOOKUP["by"].get_mass(monoisotopic)

    @property
    def formula(self) -> str:
        return composition_to_proforma_formula_string(self.composition)

    @property
    def composition(self) -> Counter[ElementInfo]:
        comp: Counter[ElementInfo] = FRAGMENT_ION_LOOKUP[self.amino_acid].composition + FRAGMENT_ION_LOOKUP["by"].composition
        if comp is None:
            raise ValueError(f"Composition not available for immonium ion of amino acid: {self.amino_acid}")
        return comp


@dataclass(frozen=True, slots=True)
class ReferenceIon(Serializable, CompositionProvider, MassProvider):
    """Represents a reference ion"""

    name: str

    _cache: ClassVar[dict[tuple, "ReferenceIon"]] = {}

    def __new__(cls, name: str):
        """Create or retrieve cached instance"""
        key = (name,)
        if key not in cls._cache:
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

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

    @staticmethod
    def parse(s: str) -> "ReferenceIon":
        """Parse reference ion string like 'r[Phospho]'"""
        s = s.strip()
        match = re.match(r"r\[([^\]]+)\]", s)
        if not match:
            raise ValueError(f"Invalid reference ion: '{s}'")
        return ReferenceIon(name=match.group(1))


@dataclass(frozen=True, slots=True)
class NamedCompound(Serializable, CompositionProvider, MassProvider):
    """
    Represents a named compound.

    Example: 0@_{Urocanic Acid}
    """

    name: str

    _cache: ClassVar[dict[tuple, "NamedCompound"]] = {}

    def __new__(cls, name: str):
        """Create or retrieve cached instance"""
        key = (name,)
        if key not in cls._cache:
            # Evict oldest entry if cache is full
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

    def mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError("Mass calculation for NamedCompound is not implemented")

    @property
    def composition(self) -> Counter[ElementInfo]:
        raise NotImplementedError("Composition calculation for NamedCompound is not implemented")

    def serialize(self) -> str:
        return f"_{{{self.name}}}"

    @staticmethod
    def parse(s: str) -> "NamedCompound":
        """Parse named compound string like '_{Urocanic Acid}'"""
        s = s.strip()
        match = re.match(r"_\{([^\}]+)\}", s)
        if not match:
            raise ValueError(f"Invalid named compound: '{s}'")
        return NamedCompound(name=match.group(1))


@dataclass(frozen=True, slots=True)
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

    @staticmethod
    def parse(s: str) -> "ChemicalFormula":
        """Parse chemical formula string like 'f{C13H9}'"""
        s = s.strip()
        match = re.match(r"f\{([^\}]+)\}", s)
        if not match:
            raise ValueError(f"Invalid chemical formula: '{s}'")
        return ChemicalFormula(formula=match.group(1))


@dataclass(frozen=True, slots=True)
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

    @staticmethod
    def parse(s: str) -> "SMILESCompound":
        """Parse SMILES compound string like 's{CN=C=O}'"""
        s = s.strip()
        match = re.match(r"s\{([^\}]+)\}", s)
        if not match:
            raise ValueError(f"Invalid SMILES compound: '{s}'")
        return SMILESCompound(smiles=match.group(1))

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
                raise ValueError(f"Unknown element '*' in SMILES '{self.smiles}'. Ensure all atoms are properly specified.")
            elem_counts[elem] += 1

        return Counter({ELEMENT_LOOKUP[elem]: count for elem, count in elem_counts.items()})

    @property
    def proforma_formula(self) -> str:
        return composition_to_proforma_formula_string(self.composition)

    @property
    def formula(self) -> str:
        return f"+{self.proforma_formula}"


@dataclass(frozen=True, slots=True)
class UnknownIon(Serializable, CompositionProvider, MassProvider):
    """Represents an unknown/unannotated ion"""

    label: int | None = None

    _cache: ClassVar[dict[tuple, "UnknownIon"]] = {}

    def __new__(cls, label: int | None = None):
        """Create or retrieve cached instance"""
        key = (label,)
        if key not in cls._cache:
            # Evict oldest entry if cache is full
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

    def mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError("Mass calculation for UnknownIon is not implemented")

    @property
    def composition(self) -> Counter[ElementInfo]:
        raise NotImplementedError("Composition calculation for UnknownIon is not implemented")

    def serialize(self) -> str:
        if self.label is not None:
            return f"?{self.label}"
        return "?"

    @staticmethod
    def parse(s: str) -> "UnknownIon":
        """Parse unknown ion string like '?' or '?5'"""
        s = s.strip()
        if s == "?":
            return UnknownIon(label=None)
        match = re.match(r"\?(\d+)", s)
        if not match:
            raise ValueError(f"Invalid unknown ion: '{s}'")
        return UnknownIon(label=int(match.group(1)))


@dataclass(frozen=True, slots=True)
class PrecursorIon(Serializable, CompositionProvider, MassProvider):
    """Represents a precursor ion"""

    _cache: ClassVar[dict[tuple, "PrecursorIon"]] = {}

    def __new__(cls):
        """Create or retrieve cached instance - singleton pattern"""
        key = ()
        if key not in cls._cache:
            # Evict oldest entry if cache is full (won't happen for singleton but keeping pattern consistent)
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
            instance = object.__new__(cls)
            cls._cache[key] = instance
        return cls._cache[key]

    def serialize(self) -> str:
        return "p"

    @staticmethod
    def parse(s: str) -> "PrecursorIon":
        """Parse precursor ion string 'p'"""
        s = s.strip()
        if s != "p":
            raise ValueError(f"Invalid precursor ion: '{s}'")
        return PrecursorIon()

    def mass(self, monoisotopic: bool = True) -> float:
        return FRAGMENT_ION_LOOKUP["p"].get_mass(monoisotopic)

    @property
    def formula(self) -> str | None:
        return FRAGMENT_ION_LOOKUP["p"].formula

    @property
    def composition(self) -> Counter[ElementInfo]:
        return FRAGMENT_ION_LOOKUP["p"].composition


# Type aliases for cleaner code
IonType = PeptideIon | InternalFragment | ImmoniumIon | ReferenceIon | NamedCompound | ChemicalFormula | SMILESCompound | UnknownIon | PrecursorIon
