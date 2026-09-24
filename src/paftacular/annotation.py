from __future__ import annotations

import warnings
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, TypedDict, Unpack

if TYPE_CHECKING:
    import peptacular as pt
else:
    try:
        import peptacular as pt
    except ImportError:
        pt = None
from tacular import AA_LOOKUP, ELEMENT_LOOKUP, ElementInfo

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
    composition_to_formula_string,
    composition_to_proforma_formula_string,
    formula_to_composition,
)
from .comps.ions import SIDE_CHAIN_SERIES
from .constants import INTERNAL_SERIES_TO_DIFF, AminoAcids, InternalSeries, IonSeries
from .util import format_number, validate_integer, validate_number

# mzPAF 1.0.1 section 4.4.3 side-chain ions keep the other n-1 residues, the backbone
# part of residue n and the beta-carbon substituent that is not lost. d keeps C2H3N,
# w keeps C3H3O2 and v keeps C2H3NO2 (the whole side chain is lost).
_SIDE_CHAIN_BACKBONE = {"d": "C2H3N", "w": "C3H3O2"}
_BETA_SUBSTITUENT: dict[tuple[str, str], str] = {
    **{(series, residue): "H" for series in ("d", "w") for residue in "CDEFHKLMNOQRSUWY"},
    ("d", "V"): "CH3",
    ("w", "V"): "CH3",
    **{(f"{series}a", "T"): "OH" for series in ("d", "w")},
    **{(f"{series}b", "T"): "CH3" for series in ("d", "w")},
    **{(f"{series}a", "I"): "C2H5" for series in ("d", "w")},
    **{(f"{series}b", "I"): "CH3" for series in ("d", "w")},
}
_V_ION_RESIDUES = "ACDEFGHIKLMNOPQRSTUVWY"


def _require_peptacular() -> None:
    if pt is None:
        raise ImportError("peptacular is required for this feature. Install it with: pip install paftacular[peptacular]")


class CommonAnnotationParams(TypedDict, total=False):
    """Common parameters shared across all annotation factory methods"""

    analyte_reference: int | None
    is_auxiliary: bool
    neutral_losses: list[NeutralLoss | str] | None
    isotopes: list[IsotopeSpecification | str | int] | None
    adducts: list[Adduct | str] | None
    charge: int
    mass_error: float | None
    mass_error_unit: Literal["da", "ppm"]
    confidence: float | None


@dataclass(frozen=True, slots=True)
class PafAnnotation:
    """Fragment ion annotation following mzPAF specification"""

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
    resolved_sequence: str | None = None

    def __post_init__(self):
        """Validate annotation constraints"""
        validate_integer(self.charge, "Charge", 1)
        if self.analyte_reference is not None:
            validate_integer(self.analyte_reference, "Analyte reference")
        if self.confidence is not None:
            validate_number(self.confidence)
        if self.resolved_sequence is not None:
            if not isinstance(self.ion_type, PeptideIon | InternalFragment | PrecursorIon):
                raise ValueError("Resolved sequence requires a peptide, internal, or precursor ion")
            if not isinstance(self.resolved_sequence, str) or not self.resolved_sequence:
                raise ValueError("Resolved sequence must be a nonempty string")
            if isinstance(self.ion_type, PeptideIon | InternalFragment):
                embedded = self.ion_type.sequence
                if embedded is not None and embedded != self.resolved_sequence:
                    raise ValueError("Embedded and resolved sequences must agree")
        if self.confidence is not None and not (0.0 <= self.confidence <= 1.0):
            raise ValueError(f"Confidence must be between 0.0 and 1.0, got {self.confidence}")

    @property
    def peptacular_ion_type(self) -> pt.IonType | None:
        """Map to peptacular IonType if applicable, else None"""
        _require_peptacular()
        if isinstance(self.ion_type, PeptideIon):
            series_map = {
                IonSeries.A: pt.IonType.A,
                IonSeries.B: pt.IonType.B,
                IonSeries.C: pt.IonType.C,
                IonSeries.X: pt.IonType.X,
                IonSeries.Y: pt.IonType.Y,
                IonSeries.Z: pt.IonType.Z_RADICAL,  # mzPAF z is the z-dot radical
            }
            return series_map.get(self.ion_type.series, None)
        elif isinstance(self.ion_type, PrecursorIon):
            return pt.IonType.PRECURSOR
        elif isinstance(self.ion_type, ImmoniumIon):
            return pt.IonType.IMMONIUM
        elif isinstance(self.ion_type, InternalFragment):
            return pt.IonType(self.ion_type._fragment_ion_key)
        return None

    @staticmethod
    def _create_annotation(ion_type: IonType, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Internal factory method that handles common parameter conversion"""

        # Convert neutral losses from strings if necessary
        nl_objs: list[NeutralLoss] = []
        if neutral_losses := kwargs.get("neutral_losses"):
            for nl in neutral_losses:
                nl_objs.append(nl if isinstance(nl, NeutralLoss) else NeutralLoss.parse(nl))

        # Convert isotopes from strings if necessary
        iso_objs: list[IsotopeSpecification] = []
        if isotopes := kwargs.get("isotopes"):
            for iso in isotopes:
                match iso:
                    case int():
                        iso_objs.append(IsotopeSpecification(count=iso))
                    case str():
                        iso_objs.append(IsotopeSpecification.parse(iso))
                    case IsotopeSpecification():
                        iso_objs.append(iso)
                    case _:
                        raise ValueError(f"Invalid isotope specification: {iso}")

        # Convert adducts from strings if necessary
        adduct_objs: list[Adduct] = []
        if adducts := kwargs.get("adducts"):
            for ad in adducts:
                adduct_objs.append(ad if isinstance(ad, Adduct) else Adduct.parse(ad))

        # Create MassError object if necessary
        mass_error_obj: MassError | None = None
        if (mass_error_val := kwargs.get("mass_error")) is not None:
            mass_error_obj = MassError(value=mass_error_val, unit=kwargs.get("mass_error_unit", "da"))

        return PafAnnotation(
            ion_type=ion_type,
            analyte_reference=kwargs.get("analyte_reference"),
            is_auxiliary=kwargs.get("is_auxiliary", False),
            neutral_losses=tuple(nl_objs),
            isotopes=tuple(iso_objs),
            adducts=tuple(adduct_objs),
            charge=kwargs.get("charge", 1),
            mass_error=mass_error_obj,
            confidence=kwargs.get("confidence"),
        )

    @staticmethod
    def make_precursor(**kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for the precursor ion"""
        return PafAnnotation._create_annotation(PrecursorIon(), **kwargs)

    @staticmethod
    def make_peptide(ion_type: str | IonSeries, position: int, sequence: str | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a peptide fragment ion"""
        ion_series: IonSeries = ion_type if isinstance(ion_type, IonSeries) else IonSeries(ion_type)
        return PafAnnotation._create_annotation(PeptideIon(series=ion_series, position=position, sequence=sequence), **kwargs)

    @staticmethod
    def make_internal(
        start_position: int, end_position: int, ion_type: str | InternalSeries = "by", sequence: str | None = None, **kwargs: Unpack[CommonAnnotationParams]
    ) -> PafAnnotation:
        """Create a PafAnnotation for an internal fragment"""
        internal_ion = InternalFragment(
            start_position=start_position,
            end_position=end_position,
            sequence=sequence,
        )

        ion_type_enum = InternalSeries(ion_type) if isinstance(ion_type, str) else ion_type

        # Add series-specific neutral loss if applicable
        if series_loss := INTERNAL_SERIES_TO_DIFF[ion_type_enum]:
            neutral_losses = list(kwargs.get("neutral_losses") or [])
            neutral_losses.append(NeutralLoss.parse(series_loss))
            kwargs["neutral_losses"] = neutral_losses

        return PafAnnotation._create_annotation(internal_ion, **kwargs)

    @staticmethod
    def make_immonium(amino_acid: str | AminoAcids, modification: str | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for an immonium ion"""
        aa_enum: AminoAcids = amino_acid if isinstance(amino_acid, AminoAcids) else AminoAcids(amino_acid)
        return PafAnnotation._create_annotation(ImmoniumIon(amino_acid=aa_enum, modification=modification), **kwargs)

    @staticmethod
    def make_reference(name: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a reference ion"""
        return PafAnnotation._create_annotation(ReferenceIon(name=name), **kwargs)

    @staticmethod
    def make_named_compound(name: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a named compound"""
        return PafAnnotation._create_annotation(NamedCompound(name=name), **kwargs)

    @staticmethod
    def make_formula(formula: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a chemical formula"""
        return PafAnnotation._create_annotation(ChemicalFormula(formula=formula), **kwargs)

    @staticmethod
    def make_smiles(smiles: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a SMILES compound"""
        return PafAnnotation._create_annotation(SMILESCompound(smiles=smiles), **kwargs)

    @staticmethod
    def make_unknown(label: int | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for an unknown/unannotated ion"""
        return PafAnnotation._create_annotation(UnknownIon(label=label), **kwargs)

    def _parse_sequence(self, sequence: str) -> pt.ProFormaAnnotation:
        """Parse the fragment sequence and warn when its length disagrees with the ion position."""
        _require_peptacular()
        annot = pt.parse(sequence)
        if annot.has_charge:
            raise ValueError("Sequence in annotation should not have charge for mass calculation")
        ion = self.ion_type
        if isinstance(ion, PeptideIon):
            expected = ion.position
        elif isinstance(ion, InternalFragment):
            expected = ion.end_position - ion.start_position + 1
        else:
            return annot
        residues = len(annot)
        if residues != expected:
            # mzPAF 1.0.1 section 4.4.3: the sequence length MUST NOT be less than, and SHOULD NOT be greater than, the ordinal.
            warnings.warn(
                f"The embedded sequence {sequence!r} has {residues} residues but {ion.serialize(include_sequence=False)} spans {expected}. "
                f"The calculation uses all {residues} residues.",
                UserWarning,
                stacklevel=4,
            )
        return annot

    def _side_chain_sequence(self, annot: pt.ProFormaAnnotation) -> tuple[pt.ProFormaAnnotation, Counter[ElementInfo]]:
        """Split a side-chain ion into the residues it keeps whole and the composition of the rest.

        Returns the sequence to sum (residue n included, and for v its modifications removed)
        and the composition to add to it, which replaces residue n by the part the ion keeps.
        """
        ion = self.ion_type
        assert isinstance(ion, PeptideIon)
        series = str(ion.series)
        index = len(annot) - 1 if series.startswith("d") else 0
        residue = annot.stripped_sequence[index]
        label = f"{series}{ion.position}"
        if series == "v":
            if residue not in _V_ION_RESIDUES:
                raise ValueError(f"{label} is not defined for residue {residue}")
            # The side chain leaves whole, so a modification on it leaves too.
            annot = annot.copy()
            annot.clear_internal_mod_at_index(index)
            kept = formula_to_composition("C2H3NO2")
        else:
            substituent = _BETA_SUBSTITUENT.get((series, residue))
            if substituent is None:
                hint = f" Use {series}a or {series}b." if residue in "TI" and len(series) == 1 else ""
                raise ValueError(f"{label} is not defined for residue {residue}.{hint}")
            if annot.has_internal_mods_at_index(index):
                raise ValueError(f"{label} is not defined when residue {residue} carries a modification")
            kept = formula_to_composition(_SIDE_CHAIN_BACKBONE[series[0]])
            kept.update(formula_to_composition(substituent))
        kept.subtract(AA_LOOKUP[residue].composition)
        return annot, kept

    def _ion_parts(self, calculate_sequence: bool) -> tuple[pt.ProFormaAnnotation | None, Counter[ElementInfo] | None]:
        """Return the sequence to add, and the ion composition when a side-chain ion replaces the ion offset."""
        if calculate_sequence is not True or self.sequence is None:
            return None, None
        annot = self._parse_sequence(self.sequence)
        if isinstance(self.ion_type, PeptideIon) and self.ion_type.series in SIDE_CHAIN_SERIES:
            return self._side_chain_sequence(annot)
        return annot, None

    def mass(self, monoisotopic: bool = True, calculate_sequence: bool = True) -> float:
        """Calculate the mass of the annotated ion including modifications"""
        annot, ion_comp = self._ion_parts(calculate_sequence)
        if ion_comp is not None:
            base_mass = sum(element.get_mass(monoisotopic) * count for element, count in ion_comp.items())
        else:
            base_mass = self.ion_type.mass(monoisotopic=monoisotopic)

        # Apply neutral losses/gains
        for loss in self.neutral_losses:
            base_mass += loss.mass(monoisotopic=monoisotopic)

        if isinstance(self.ion_type, ChemicalFormula):
            # A ChemicalFormula's atom count already represents the fully charged species
            # (mzPAF section 4.4.9): any adduct only labels which atoms carry the charge and
            # MUST NOT add mass, and the theoretical m/z needs only an electron-mass correction
            # per charge, not a full proton per charge like the other (neutral-basis) ion types.
            base_mass -= self.charge * 0.000548579909
        else:
            # Apply adducts
            for adduct in self.adducts:
                base_mass += adduct.mass(monoisotopic=monoisotopic)
            if self.adducts:
                base_mass -= self.charge * 0.000548579909

            # Adjust for charge state (if no adducts specified) default protonation/deprotonation
            if self.charge != 0 and len(self.adducts) == 0:
                base_mass += self.charge * 1.007276466812

        # Apply isotopes
        for isotope in self.isotopes:
            base_mass += isotope.mass(monoisotopic=monoisotopic)

        if annot is not None:
            base_mass += annot.mass(monoisotopic=monoisotopic, ion_type="n")

        return base_mass

    def mz(self, monoisotopic: bool = True, calculate_sequence: bool = True) -> float:
        """Calculate the m/z of the annotated ion"""
        total_mass = self.mass(monoisotopic=monoisotopic, calculate_sequence=calculate_sequence)
        return total_mass / self.charge

    def comp(self, calculate_sequence: bool = True) -> Counter[ElementInfo]:
        """Calculate the elemental composition of the annotated ion including modifications"""
        comp: Counter[ElementInfo] = Counter()
        annot, ion_comp = self._ion_parts(calculate_sequence)

        # Base ion composition
        comp.update(self.ion_type.composition if ion_comp is None else ion_comp)

        # Apply neutral losses/gains
        for loss in self.neutral_losses:
            comp.update(loss.composition)

        if isinstance(self.ion_type, ChemicalFormula):
            # A ChemicalFormula's atom count already represents the fully charged species
            # (mzPAF section 4.4.9): any adduct only labels which atoms carry the charge and
            # MUST NOT add atoms on top, and there's no separate proton to add for the charge
            # itself (electrons aren't atoms, so there's nothing to add to the composition).
            pass
        else:
            # Apply adducts
            for adduct in self.adducts:
                comp.update(adduct.composition)

            # Adjust for charge state (if no adducts specified) default protonation/deprotonation
            if self.charge != 0 and len(self.adducts) == 0:
                proton = ELEMENT_LOOKUP["H"]
                comp[proton] += self.charge

        # Apply isotopes
        for isotope in self.isotopes:
            comp.update(isotope.composition)

        if annot is not None:
            comp.update(annot.comp(ion_type="n"))

        # Consume ordinary atoms when an isotope delta removes the monoisotope.
        # Keep genuine deficits when the annotation provides insufficient context.
        for element, count in tuple(comp.items()):
            if count < 0 and element == ELEMENT_LOOKUP.get_monoisotopic(element.symbol):
                ordinary = ELEMENT_LOOKUP[element.symbol]
                if ordinary != element:
                    available = min(-count, max(0, comp[ordinary]))
                    comp[element] += available
                    comp[ordinary] -= available

        return comp

    def dict_composition(self, calculate_sequence: bool = True) -> dict[str, int]:
        """Get the elemental composition as a dictionary of element symbols to counts"""
        comp_counter = self.comp(calculate_sequence=calculate_sequence)
        return {str(elem): count for elem, count in comp_counter.items()}

    @property
    def sequence(self) -> str | None:
        """Get the peptide sequence if applicable, else None"""
        if self.resolved_sequence is not None:
            return self.resolved_sequence
        if isinstance(self.ion_type, PeptideIon):
            return self.ion_type.sequence
        elif isinstance(self.ion_type, InternalFragment):
            return self.ion_type.sequence
        return None

    def resolve(self, analytes: str | Mapping[int, str]) -> PafAnnotation:
        """Resolve this peptide, internal, or precursor ion against analyte context."""
        from .resolution import resolve

        return resolve(self, analytes)

    def to_dict(self) -> dict:
        """Export the versioned, reversible structured representation."""
        from .serialization import to_dict

        return to_dict(self)

    @staticmethod
    def from_dict(data: Mapping[str, object]) -> PafAnnotation:
        """Validate and reconstruct a versioned structured representation."""
        from .serialization import from_dict

        return from_dict(data)

    def formula(self, calculate_sequence: bool = True) -> str:
        """Get the chemical formula string of the annotated ion"""
        return composition_to_formula_string(self.comp(calculate_sequence=calculate_sequence))

    def proforma_formula(self, calculate_sequence: bool = True) -> str:
        """Get the ProForma-style chemical formula string of the annotated ion"""
        return composition_to_proforma_formula_string(self.comp(calculate_sequence=calculate_sequence))

    def serialize(self, include_sequence: bool = True) -> str:
        """Serialize the annotation back to mzPAF string format"""
        parts: list[str] = []

        # Auxiliary marker
        if self.is_auxiliary:
            parts.append("&")

        # Analyte reference
        if self.analyte_reference is not None:
            parts.append(f"{self.analyte_reference}@")

        # Ion type
        # check if ion type has a sequence attribute
        match self.ion_type:
            case PeptideIon() | InternalFragment():
                parts.append(self.ion_type.serialize(include_sequence=include_sequence))
            case _:
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
            parts.append(f"*{format_number(self.confidence)}")

        return "".join(parts)

    @staticmethod
    def parse(annotation_str: str) -> PafAnnotation:
        """Parse a single mzPAF annotation string into a FragmentAnnotation object"""
        from .parser import MZ_PAF_PARSER

        return MZ_PAF_PARSER.parse(annotation_str)

    def as_dict(self) -> dict:
        """Convert the annotation to a dictionary representation"""
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
