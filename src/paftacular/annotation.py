from collections import Counter
from dataclasses import asdict, dataclass
from typing import Literal, TypedDict, Unpack

from tacular import ELEMENT_LOOKUP, ElementInfo

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
)
from .constants import INTERNAL_SERIES_TO_DIFF, AminoAcids, InternalSeries, IonSeries


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

    def __post_init__(self):
        """Validate annotation constraints"""
        if self.charge < 1:
            raise ValueError(f"Charge must be >= 1, got {self.charge}")
        if self.confidence is not None and not (0.0 <= self.confidence <= 1.0):
            raise ValueError(f"Confidence must be between 0.0 and 1.0, got {self.confidence}")

    @staticmethod
    def _create_annotation(ion_type: IonType, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
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
    def make_precursor(**kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for the precursor ion"""
        return PafAnnotation._create_annotation(PrecursorIon(), **kwargs)

    @staticmethod
    def make_peptide(ion_type: str | IonSeries, position: int, sequence: str | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for a peptide fragment ion"""
        ion_series: IonSeries = ion_type if isinstance(ion_type, IonSeries) else IonSeries(ion_type)
        return PafAnnotation._create_annotation(PeptideIon(series=ion_series, position=position, sequence=sequence), **kwargs)

    @staticmethod
    def make_internal(
        start_position: int, end_position: int, ion_type: str | InternalSeries = "by", sequence: str | None = None, **kwargs: Unpack[CommonAnnotationParams]
    ) -> "PafAnnotation":
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
    def make_immonium(amino_acid: str | AminoAcids, modification: str | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for an immonium ion"""
        aa_enum: AminoAcids = amino_acid if isinstance(amino_acid, AminoAcids) else AminoAcids(amino_acid)
        return PafAnnotation._create_annotation(ImmoniumIon(amino_acid=aa_enum, modification=modification), **kwargs)

    @staticmethod
    def make_reference(name: str, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for a reference ion"""
        return PafAnnotation._create_annotation(ReferenceIon(name=name), **kwargs)

    @staticmethod
    def make_named_compound(name: str, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for a named compound"""
        return PafAnnotation._create_annotation(NamedCompound(name=name), **kwargs)

    @staticmethod
    def make_formula(formula: str, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for a chemical formula"""
        return PafAnnotation._create_annotation(ChemicalFormula(formula=formula), **kwargs)

    @staticmethod
    def make_smiles(smiles: str, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for a SMILES compound"""
        return PafAnnotation._create_annotation(SMILESCompound(smiles=smiles), **kwargs)

    @staticmethod
    def make_unknown(label: int | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> "PafAnnotation":
        """Create a PafAnnotation for an unknown/unannotated ion"""
        return PafAnnotation._create_annotation(UnknownIon(label=label), **kwargs)

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
            import peptacular as pt

            annot = pt.parse(self.sequence)

            if annot.has_charge:  # ty: ignore
                raise ValueError("Sequence in annotation should not have charge for mass calculation")

            sequence_mass = annot.mass(monoisotopic=monoisotopic, ion_type="n")  # type: ignore
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
            import peptacular as pt

            annot = pt.parse(self.sequence)

            if annot.has_charge:  # ty: ignore
                raise ValueError("Sequence in annotation should not have charge for mass calculation")

            seq_comp = annot.comp()  # type: ignore
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
        from .parser import MZ_PAF_PARSER

        return MZ_PAF_PARSER.parse(annotation_str)

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
