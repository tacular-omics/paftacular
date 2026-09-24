from __future__ import annotations

import warnings
from collections import Counter
from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass
from functools import lru_cache
from typing import TYPE_CHECKING, Literal, TypedDict, Unpack

if TYPE_CHECKING:
    import peptacular as pt
else:
    try:
        import peptacular as pt
    except ImportError:
        pt = None
from tacular import AA_LOOKUP, ELEMENT_LOOKUP, ElementInfo
from tacular.constants import ELECTRON_MASS, PROTON_MASS

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
from .constants import _INTERNAL_SERIES_TO_DIFF, AminoAcids, InternalSeries, IonSeries
from .errors import PaftacularError, reraise_as_paftacular
from .util import format_number, to_enum, validate_integer, validate_number

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


_HYDROGEN = ELEMENT_LOOKUP["H"]
# The average charge carrier is a natural-abundance H atom less an electron.
_AVERAGE_PROTON_MASS = _HYDROGEN.get_mass(monoisotopic=False) - ELECTRON_MASS


def _require_peptacular() -> None:
    if pt is None:
        raise ImportError("peptacular is required for this feature. Install it with: pip install paftacular[peptacular]")


@lru_cache(maxsize=1024)
def _parse_proforma_cached(sequence: str) -> pt.ProFormaAnnotation:
    try:
        return pt.parse(sequence)
    except ValueError as error:
        raise PaftacularError(f"Invalid ProForma sequence {sequence!r}: {error}") from error


def _parse_proforma(sequence: str) -> pt.ProFormaAnnotation:
    """Parse a ProForma sequence with peptacular.

    The result is shared through a cache. Copy it before editing it in place.
    """
    _require_peptacular()
    return _parse_proforma_cached(sequence)


def _drop_static_mods_at(annot: pt.ProFormaAnnotation, static_mods: Mapping[int, list], index: int) -> None:
    """Write the global fixed modifications of ``annot`` in place at every site except ``index``."""
    annot.set_static_mods(None, validate=False)
    for site, mods in static_mods.items():
        if site == index:
            continue
        for mod in mods:
            for _ in range(mod.count):
                if site == -1:
                    annot.append_nterm_mod(mod.value)
                elif site == -2:
                    annot.append_cterm_mod(mod.value)
                else:
                    annot.append_internal_mod_at_index(site, mod.value)


def _isotope_label_map(annot: pt.ProFormaAnnotation | None) -> dict[ElementInfo, ElementInfo]:
    """The global isotope labels (<13C>) of a sequence as {ordinary element: isotope}."""
    if annot is None or not annot.has_isotope_mods:
        return {}
    return annot._map_isotopes()


def _sequence_mass(annot: pt.ProFormaAnnotation, monoisotopic: bool) -> float:
    """The neutral mass of a sequence's residues and modifications, at full precision.

    peptacular's fast mass path adds a Unimod modification by its tabulated mass, rounded to
    6 decimals (Oxidation 15.994915, not 15.99491462). Summing the composition instead gives
    the exact mass. A mass-only modification ([+42.010565]) has no composition, so it keeps the
    fast path.
    """
    if annot.has_mods():
        try:
            return annot.mass(monoisotopic=monoisotopic, ion_type="n", calculate_with_composition=True)
        except pt.CompositionError:
            pass
    return annot.mass(monoisotopic=monoisotopic, ion_type="n")


@lru_cache(maxsize=4096)
def _cached_sequence_mass(sequence: str, monoisotopic: bool) -> float:
    return _sequence_mass(_parse_proforma(sequence), monoisotopic)


def _removed_carrier_atoms(charge: int, adducts: tuple[Adduct, ...]) -> Counter[ElementInfo]:
    """The atoms charge carriers remove from a neutral ion, as negative counts.

    That is one H per charge for a default negative charge, and the atoms of each negative
    adduct (``[M-H]``). A global isotope label covers these atoms, so ``<2H>`` at ``^-1``
    loses a deuteron.
    """
    removed: Counter[ElementInfo] = Counter()
    if not adducts:
        if charge < 0:
            removed[_HYDROGEN] += charge
        return removed
    for adduct in adducts:
        if adduct.count < 0:
            removed.update(adduct.composition)
    return removed


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
    """Fragment ion annotation following mzPAF specification.

    ``charge`` is any nonzero integer. A negative charge means a deprotonated (negative mode)
    ion: without adducts it removes ``abs(charge)`` protons, and ``mz()`` divides by ``abs(charge)``.
    """

    # Core ion description
    ion_type: IonType
    _: KW_ONLY

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
        if type(self.charge) is not int or self.charge == 0:
            raise PaftacularError(f"Charge must be a nonzero integer, got {self.charge!r}")
        if self.analyte_reference is not None:
            validate_integer(self.analyte_reference, "Analyte reference")
        if self.confidence is not None:
            validate_number(self.confidence)
        if self.resolved_sequence is not None:
            if not isinstance(self.ion_type, PeptideIon | InternalFragment | PrecursorIon):
                raise PaftacularError("Resolved sequence requires a peptide, internal, or precursor ion")
            if not isinstance(self.resolved_sequence, str) or not self.resolved_sequence:
                raise PaftacularError("Resolved sequence must be a nonempty string")
            if isinstance(self.ion_type, PeptideIon | InternalFragment):
                embedded = self.ion_type.sequence
                if embedded is not None and embedded != self.resolved_sequence:
                    raise PaftacularError("Embedded and resolved sequences must agree")
        if self.confidence is not None and not (0.0 <= self.confidence <= 1.0):
            raise PaftacularError(f"Confidence must be between 0.0 and 1.0, got {self.confidence}")

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
                        iso_objs.append(IsotopeSpecification(iso))
                    case str():
                        iso_objs.append(IsotopeSpecification.parse(iso))
                    case IsotopeSpecification():
                        iso_objs.append(iso)
                    case _:
                        raise PaftacularError(f"Invalid isotope specification: {iso}")

        # Convert adducts from strings if necessary
        adduct_objs: list[Adduct] = []
        if adducts := kwargs.get("adducts"):
            for ad in adducts:
                adduct_objs.append(ad if isinstance(ad, Adduct) else Adduct.parse(ad))

        # Create MassError object if necessary
        mass_error_obj: MassError | None = None
        if (mass_error_val := kwargs.get("mass_error")) is not None:
            mass_error_obj = MassError(mass_error_val, unit=kwargs.get("mass_error_unit", "da"))

        return PafAnnotation(
            ion_type,
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
    def make_peptide(ion_type: str | IonSeries, position: int, *, sequence: str | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a peptide fragment ion"""
        return PafAnnotation._create_annotation(PeptideIon(to_enum(IonSeries, ion_type, "ion series"), position, sequence=sequence), **kwargs)

    @staticmethod
    def make_internal(
        start_position: int,
        end_position: int,
        *,
        ion_type: str | InternalSeries = "by",
        sequence: str | None = None,
        **kwargs: Unpack[CommonAnnotationParams],
    ) -> PafAnnotation:
        """Create a PafAnnotation for an internal fragment"""
        internal_ion = InternalFragment(start_position, end_position, sequence=sequence)
        ion_type_enum = to_enum(InternalSeries, ion_type, "internal series")

        # Add series-specific neutral loss if applicable
        if series_loss := _INTERNAL_SERIES_TO_DIFF[ion_type_enum]:
            neutral_losses = list(kwargs.get("neutral_losses") or [])
            neutral_losses.append(NeutralLoss.parse(series_loss))
            kwargs["neutral_losses"] = neutral_losses

        return PafAnnotation._create_annotation(internal_ion, **kwargs)

    @staticmethod
    def make_immonium(amino_acid: str | AminoAcids, *, modification: str | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for an immonium ion"""
        return PafAnnotation._create_annotation(ImmoniumIon(to_enum(AminoAcids, amino_acid, "immonium amino acid"), modification=modification), **kwargs)

    @staticmethod
    def make_reference(name: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a reference ion"""
        return PafAnnotation._create_annotation(ReferenceIon(name), **kwargs)

    @staticmethod
    def make_named_compound(name: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a named compound"""
        return PafAnnotation._create_annotation(NamedCompound(name), **kwargs)

    @staticmethod
    def make_formula(formula: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a chemical formula"""
        return PafAnnotation._create_annotation(ChemicalFormula(formula), **kwargs)

    @staticmethod
    def make_smiles(smiles: str, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for a SMILES compound"""
        return PafAnnotation._create_annotation(SMILESCompound(smiles), **kwargs)

    @staticmethod
    def make_unknown(*, label: int | None = None, **kwargs: Unpack[CommonAnnotationParams]) -> PafAnnotation:
        """Create a PafAnnotation for an unknown/unannotated ion"""
        return PafAnnotation._create_annotation(UnknownIon(label=label), **kwargs)

    def _parse_sequence(self, sequence: str) -> pt.ProFormaAnnotation:
        """Parse the fragment sequence and warn when its length disagrees with the ion position."""
        annot = _parse_proforma(sequence)
        if annot.has_charge:
            raise PaftacularError("Sequence in annotation should not have charge for mass calculation")
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
        # A global fixed modification (<[Carbamidomethyl]@C>) sits on the side chain just like
        # an explicit one, so both follow the same rule.
        static_mods = annot.map_static_mods_to_indexes() if annot.has_static_mods else {}
        if series == "v":
            if residue not in _V_ION_RESIDUES:
                raise PaftacularError(f"{label} is not defined for residue {residue}")
            # The side chain leaves whole, so a modification on it leaves too.
            annot = annot.copy()
            annot.clear_internal_mod_at_index(index)
            if index in static_mods:
                _drop_static_mods_at(annot, static_mods, index)
            kept = formula_to_composition("C2H3NO2")
        else:
            substituent = _BETA_SUBSTITUENT.get((series, residue))
            if substituent is None:
                hint = f" Use {series}a or {series}b." if residue in "TI" and len(series) == 1 else ""
                raise PaftacularError(f"{label} is not defined for residue {residue}.{hint}")
            if annot.has_internal_mods_at_index(index) or index in static_mods:
                raise PaftacularError(f"{label} is not defined when residue {residue} carries a modification")
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

    def _offset_and_delta_comp(self, ion_comp: Counter[ElementInfo] | None, *, strict: bool = False) -> Counter[ElementInfo]:
        """The composition of the ion offset (or side-chain remnant) plus formula and reference deltas.

        Mass-only deltas have no atoms. They are skipped, or raise with ``strict``.
        """
        comp: Counter[ElementInfo] = Counter()
        comp.update(self.ion_type.composition if ion_comp is None else ion_comp)
        for loss in self.neutral_losses:
            if strict or loss.loss_type != "mass":
                comp.update(loss.composition)
        return comp

    def _removed_carrier_comp(self) -> Counter[ElementInfo]:
        """The atoms the charge carriers remove from the neutral ion. A formula ion has none."""
        if isinstance(self.ion_type, ChemicalFormula):
            return Counter()
        return _removed_carrier_atoms(self.charge, self.adducts)

    @reraise_as_paftacular
    def get_mass(self, *, monoisotopic: bool = True, calculate_sequence: bool = True) -> float:
        """Calculate the charged-species mass of the annotated ion including modifications.

        Without an embedded or resolved sequence, peptide, internal and precursor ions give only
        the ion offset plus modifiers. Pass ``calculate_sequence=False`` to get that offset even
        when a sequence is present.
        """
        annot, ion_comp = self._ion_parts(calculate_sequence)
        if ion_comp is not None:
            base_mass = sum(element.get_mass(monoisotopic=monoisotopic) * count for element, count in ion_comp.items())
        else:
            base_mass = self.ion_type.get_mass(monoisotopic=monoisotopic)

        # Apply neutral losses/gains
        for loss in self.neutral_losses:
            base_mass += loss.get_mass(monoisotopic=monoisotopic)

        # A global isotope label (<13C>) replaces its element everywhere in the neutral ion,
        # the ion offset and formula deltas included, like peptacular. Atoms a charge carrier
        # removes come out of the labelled ion, so <2H> at ^-1 loses a deuteron.
        if labels := _isotope_label_map(annot):
            labelled = self._offset_and_delta_comp(ion_comp)
            labelled.update(self._removed_carrier_comp())
            for element, count in labelled.items():
                if (isotope := labels.get(element)) is not None:
                    base_mass += (isotope.get_mass(monoisotopic=monoisotopic) - element.get_mass(monoisotopic=monoisotopic)) * count

        proton = PROTON_MASS if monoisotopic else _AVERAGE_PROTON_MASS
        if isinstance(self.ion_type, ChemicalFormula):
            # A ChemicalFormula's atom count already represents the fully charged species
            # (mzPAF section 4.4.9): any adduct only labels which atoms carry the charge and
            # MUST NOT add mass, and the theoretical m/z needs only an electron-mass correction
            # per charge, not a full proton per charge like the other (neutral-basis) ion types.
            base_mass -= self.charge * ELECTRON_MASS
        elif self.adducts:
            # An H carrier with the sign of the charge is a proton (or a removed proton), so
            # [M+H] and [M-H]^-1 give exactly the mass of the default charge. [M+H]^-1 adds a
            # hydrogen atom and an electron.
            protons = 0
            for adduct in self.adducts:
                if adduct.base_formula == "H" and (adduct.count > 0) == (self.charge > 0):
                    protons += adduct.count
                else:
                    base_mass += adduct.get_mass(monoisotopic=monoisotopic)
            base_mass += protons * proton - (self.charge - protons) * ELECTRON_MASS
        else:
            # Default protonation (positive charge) or deprotonation (negative charge)
            base_mass += self.charge * proton

        # Apply isotopes
        for isotope in self.isotopes:
            base_mass += isotope.get_mass(monoisotopic=monoisotopic)

        if annot is not None:
            if ion_comp is None and self.sequence is not None:
                base_mass += _cached_sequence_mass(self.sequence, monoisotopic)
            else:
                base_mass += _sequence_mass(annot, monoisotopic)

        return base_mass

    def mz(self, *, monoisotopic: bool = True) -> float:
        """Calculate the m/z of the annotated ion.

        Peptide, internal and precursor ions need a sequence (embedded or from ``resolve()``).
        Without one only the ion offset is known, so this raises :class:`PaftacularError`.
        """
        if self.sequence is None and isinstance(self.ion_type, PeptideIon | InternalFragment | PrecursorIon):
            raise PaftacularError(
                f"m/z of {self.serialize()!r} needs a sequence. Embed one ({{PEPTIDE}}) or call resolve() with the analyte. "
                "get_mass() still returns the ion offset without a sequence."
            )
        return self.get_mass(monoisotopic=monoisotopic) / abs(self.charge)

    @reraise_as_paftacular
    def comp(self, *, calculate_sequence: bool = True) -> Counter[ElementInfo]:
        """Calculate the elemental composition of the annotated ion including modifications"""
        comp: Counter[ElementInfo] = Counter()
        annot, ion_comp = self._ion_parts(calculate_sequence)

        # Base ion composition and neutral losses/gains. A global isotope label (<13C>)
        # replaces its element in both, like peptacular.
        neutral = self._offset_and_delta_comp(ion_comp, strict=True)
        if labels := _isotope_label_map(annot):
            for element, count in neutral.items():
                comp[labels.get(element, element)] += count
        else:
            comp.update(neutral)

        if isinstance(self.ion_type, ChemicalFormula):
            # A ChemicalFormula's atom count already represents the fully charged species
            # (mzPAF section 4.4.9): any adduct only labels which atoms carry the charge and
            # MUST NOT add atoms on top, and there's no separate proton to add for the charge
            # itself (electrons aren't atoms, so there's nothing to add to the composition).
            pass
        else:
            # Apply adducts
            for adduct in self.adducts:
                if labels and adduct.count < 0:
                    # A removed carrier atom comes out of the labelled ion.
                    for element, count in adduct.composition.items():
                        comp[labels.get(element, element)] += count
                else:
                    comp.update(adduct.composition)

            # Adjust for charge state (if no adducts specified) default protonation/deprotonation
            if not self.adducts:
                comp[labels.get(_HYDROGEN, _HYDROGEN) if self.charge < 0 else _HYDROGEN] += self.charge

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

    def formula(self, *, calculate_sequence: bool = True) -> str:
        """Get the chemical formula string of the annotated ion"""
        return composition_to_formula_string(self.comp(calculate_sequence=calculate_sequence))

    def proforma_formula(self, *, calculate_sequence: bool = True) -> str:
        """Get the ProForma-style chemical formula string of the annotated ion"""
        return composition_to_proforma_formula_string(self.comp(calculate_sequence=calculate_sequence))

    def serialize(self, *, include_sequence: bool = True, signed_charge: bool = True) -> str:
        """Serialize the annotation back to mzPAF string format.

        A negative charge is written as ``^-n`` by default so the text round trips. mzPAF 1.0.1
        section 4.8 says the charge MUST NOT include the minus sign (negative mode is a property
        of the spectrum), so pass ``signed_charge=False`` to write only the magnitude.
        """
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

        # Charge state (charge 1 is implied)
        charge = self.charge if signed_charge else abs(self.charge)
        if charge != 1:
            parts.append(f"^{charge}")

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
        from .parser import parse

        return parse(annotation_str)

    def __str__(self) -> str:
        return self.serialize()
