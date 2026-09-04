"""Ion type definitions for mzPAF annotations"""

import re
from collections import Counter
from dataclasses import dataclass
from functools import cached_property
from typing import ClassVar

try:
    import peptacular as pt
except ImportError:
    pt = None  # type: ignore[assignment]
from tacular import AA_LOOKUP, ELEMENT_LOOKUP, FRAGMENT_ION_LOOKUP, REFMOL_LOOKUP, ElementInfo, RefMolInfo

from ..constants import MAX_CACHE_SIZE, AminoAcids, IonSeries
from ..util import validate_integer
from .base import CompositionProvider, MassProvider, Serializable
from .util import composition_to_formula_string, composition_to_proforma_formula_string, formula_to_composition


def _require_peptacular() -> None:
    if pt is None:
        raise ImportError("peptacular is required for this feature. Install it with: pip install paftacular[peptacular]")


@dataclass(frozen=True, slots=True)
class PeptideIon(Serializable, CompositionProvider, MassProvider):
    """Represents a primary peptide fragment ion"""

    series: IonSeries
    position: int
    sequence: str | None = None  # ProForma sequence

    def __post_init__(self):
        validate_integer(self.position, "Position", 1)
        IonSeries(self.series)

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

    def serialize(self, include_sequence: bool = True) -> str:
        result = f"{self.series}{self.position}"
        if include_sequence and self.sequence:
            result += f"{{{self.sequence}}}"
        return result

    @staticmethod
    def parse(s: str) -> "PeptideIon":
        """Parse peptide ion string like 'b5', 'y10{PEPTIDE}'"""
        from ..annotation import PafAnnotation
        from ..parser import parse_single

        try:
            annotation = parse_single(s)
        except ValueError as error:
            raise ValueError(f"Invalid peptide ion: {s!r}") from error
        if not isinstance(annotation.ion_type, PeptideIon) or annotation != PafAnnotation(annotation.ion_type):
            raise ValueError(f"Invalid peptide ion component: {s!r}")
        return annotation.ion_type


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

    def __post_init__(self):
        """Validate that backbone cleavage types are set together or not at all"""
        if (self.nterm_ion_type is None) != (self.cterm_ion_type is None):
            raise ValueError(
                "nterm_ion_type and cterm_ion_type must both be set or both be None, "
                f"got nterm_ion_type={self.nterm_ion_type!r}, cterm_ion_type={self.cterm_ion_type!r}"
            )
        validate_integer(self.start_position, "Start position", 1)
        validate_integer(self.end_position, "End position", self.start_position)
        if self.nterm_ion_type is not None and (
            self.nterm_ion_type not in (IonSeries.A, IonSeries.B, IonSeries.C) or self.cterm_ion_type not in (IonSeries.X, IonSeries.Y, IonSeries.Z)
        ):
            raise ValueError("Internal cleavage types must be a/b/c and x/y/z")

    @property
    def _fragment_ion_key(self) -> str:
        """tacular FRAGMENT_ION_LOOKUP key for this fragment's backbone cleavage type, e.g. 'by', 'ax'"""
        nterm = self.nterm_ion_type if self.nterm_ion_type is not None else IonSeries.B
        cterm = self.cterm_ion_type if self.cterm_ion_type is not None else IonSeries.Y
        return f"{nterm}{cterm}"

    def serialize(self, include_sequence: bool = True) -> str:
        # If using default yb cleavage, just use 'm'
        result = f"m{self.start_position}:{self.end_position}"

        if include_sequence and self.sequence:
            result += f"{{{self.sequence}}}"
        return result + self.cleavage_correction

    @property
    def cleavage_correction(self) -> str:
        """Encode the actual backbone composition as mzPAF gains and losses."""
        positive = Counter({element: count for element, count in self.composition.items() if count > 0})
        negative = Counter({element: -count for element, count in self.composition.items() if count < 0})
        gain = "+" + composition_to_formula_string(positive) if positive else ""
        loss = "-" + composition_to_formula_string(negative) if negative else ""
        return gain + loss

    @staticmethod
    def parse(s: str) -> "InternalFragment":
        """Parse internal fragment string like 'm5:10', 'm5:10{PEPTIDE}'"""
        from ..annotation import PafAnnotation
        from ..constants import InternalSeries
        from ..parser import parse_single

        try:
            annotation = parse_single(s)
        except ValueError as error:
            raise ValueError(f"Invalid internal fragment: {s!r}") from error
        ion = annotation.ion_type
        if not isinstance(ion, InternalFragment) or annotation != PafAnnotation(ion, neutral_losses=annotation.neutral_losses):
            raise ValueError(f"Invalid internal fragment component: {s!r}")
        if not annotation.neutral_losses:
            return ion
        correction: Counter[ElementInfo] = Counter()
        for loss in annotation.neutral_losses:
            correction.update(loss.composition)
        for series in InternalSeries:
            candidate = InternalFragment(ion.start_position, ion.end_position, ion.sequence, IonSeries(series[0]), IonSeries(series[1]))
            if candidate.composition == correction:
                return candidate
        raise ValueError("Neutral correction does not describe a supported internal cleavage")

    def mass(self, monoisotopic: bool = True) -> float:
        return FRAGMENT_ION_LOOKUP[self._fragment_ion_key].get_mass(monoisotopic)

    @property
    def formula(self) -> str:
        formula = FRAGMENT_ION_LOOKUP[self._fragment_ion_key].formula
        if formula is None:
            raise ValueError("Formula not available for internal fragment")
        return formula

    @property
    def composition(self) -> Counter[ElementInfo]:
        comp: Counter[ElementInfo] = FRAGMENT_ION_LOOKUP[self._fragment_ion_key].composition
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
        AminoAcids(amino_acid)
        if modification is not None and (not isinstance(modification, str) or not modification):
            raise ValueError("Modification must be a nonempty string")
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
        match = re.fullmatch(r"I([A-Z])(?:\[([^\]]+)\])?", s)
        if not match:
            raise ValueError(f"Invalid immonium ion: '{s}'")

        aa_str, modification = match.groups()
        return ImmoniumIon(amino_acid=AminoAcids(aa_str), modification=modification)

    def mass(self, monoisotopic: bool = True) -> float:
        m = 0.0
        if self.modification is not None:
            _require_peptacular()
            mod_tag: pt.ModificationTags = pt.ModificationTags.from_string(self.modification)
            m += mod_tag.get_mass(monoisotopic)

        aa_mass = AA_LOOKUP[self.amino_acid].get_mass(monoisotopic)
        if aa_mass is None:
            raise ValueError(f"Mass not available for amino acid: {self.amino_acid}")
        else:
            m += aa_mass

        m += FRAGMENT_ION_LOOKUP["i"].get_mass(monoisotopic)

        return m

    @property
    def formula(self) -> str:
        return composition_to_proforma_formula_string(self.composition)

    @property
    def composition(self) -> Counter[ElementInfo]:
        # Counter's `+`/`+=` (and unary `+`) drop any element whose total is <= 0, which would
        # make this composition silently disagree with mass() whenever a modification removes more
        # of an element than the residue+immonium supply (net-negative) or exactly cancels it
        # (net-zero). Accumulate with `.update()` (which never filters), then strip only the
        # exact-zero entries at the end -- negatives are kept so comp() stays consistent with mass().
        c: Counter[ElementInfo] = Counter()
        if self.modification is not None:
            _require_peptacular()
            mod_tag: pt.ModificationTags = pt.ModificationTags.from_string(self.modification)
            mod_comp = mod_tag.get_composition()
            if mod_comp is None:
                raise ValueError(f"Composition not available for modification: {self.modification}")
            c.update(mod_comp)

        aa_comp = AA_LOOKUP[self.amino_acid].composition
        if aa_comp is None:
            raise ValueError(f"Composition not available for amino acid: {self.amino_acid}")
        c.update(aa_comp)
        c.update(FRAGMENT_ION_LOOKUP["i"].composition)
        return Counter({el: n for el, n in c.items() if n != 0})


@dataclass(frozen=True, slots=True)
class ReferenceIon(Serializable, CompositionProvider, MassProvider):
    """Represents a reference ion"""

    name: str

    _cache: ClassVar[dict[tuple, "ReferenceIon"]] = {}

    def __new__(cls, name: str):
        """Create or retrieve cached instance"""
        if not isinstance(name, str) or not name:
            raise ValueError("Reference name must be a nonempty string")
        key = (name,)
        if key not in cls._cache:
            if len(cls._cache) >= MAX_CACHE_SIZE:
                cls._cache.pop(next(iter(cls._cache)))
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
        match = re.fullmatch(r"r\[([^\]]+)\]", s)
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
        if not isinstance(name, str) or not name:
            raise ValueError("Compound name must be a nonempty string")
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
        match = re.fullmatch(r"_\{([^\}]+)\}", s)
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
        match = re.fullmatch(r"f\{([^\}]+)\}", s)
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
        match = re.fullmatch(r"s\{([^\}]+)\}", s)
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
        if sum(mol.nodes[node_id].get("charge", 0) for node_id in mol.nodes()) != 0:
            raise ValueError("mzPAF SMILES must describe a neutral molecule. Specify charge carriers using adducts")
        for node_id in mol.nodes():
            elem = mol.nodes[node_id].get("element", "*")
            if elem == "*":
                raise ValueError(f"Unknown element '*' in SMILES '{self.smiles}'. Ensure all atoms are properly specified.")
            isotope = mol.nodes[node_id].get("isotope")
            if isotope is not None:
                elem = f"{isotope}{elem}"
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
        if label is not None:
            validate_integer(label, "Unknown ion label")
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
        match = re.fullmatch(r"\?(\d+)", s)
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
