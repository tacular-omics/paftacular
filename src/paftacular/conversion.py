from collections.abc import Mapping
from typing import Literal

import peptacular as pt

from paftacular import IonSeries

from .annotation import PafAnnotation
from .comps import (
    Adduct,
    ImmoniumIon,
    InternalFragment,
    IsotopeSpecification,
    MassError,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    UnknownIon,
)
from .constants import INTERNAL_MASS_DIFFS, AminoAcids


def to_mzpaf(
    frag: pt.Fragment,
    confidence: float | None = None,
    mass_error: float | None = None,
    mass_error_type: Literal["ppm", "da"] = "ppm",
    include_annotation: bool = True,
) -> PafAnnotation:
    """Convert fragment to mzPAF format string."""

    ion = None
    internal_loss = None
    if frag.ion_type is None:
        ion = UnknownIon()
    else:
        ion_info: pt.FragmentIonInfo = pt.FRAGMENT_ION_LOOKUP[frag.ion_type]

        match ion_info.properties:
            case pt.IonTypeProperty.FORWARD | pt.IonTypeProperty.BACKWARD:
                if not isinstance(frag.position, int):
                    position = -1
                else:
                    position = frag.position

                # some peptacular ion_types mapt to different mzPAF ion types
                match ion_info.ion_type:
                    case pt.IonType.W_VALINE:
                        ion_type = pt.IonType.W
                    case pt.IonType.D_VALINE:
                        ion_type = pt.IonType.D
                    case pt.IonType.WB_ISOLEUCINE | pt.IonType.WB_THREONINE:
                        ion_type = pt.IonType.WB
                    case pt.IonType.WA_ISOLEUCINE | pt.IonType.WA_THREONINE:
                        ion_type = pt.IonType.WA
                    case pt.IonType.DB_ISOLEUCINE | pt.IonType.DB_THREONINE:
                        ion_type = pt.IonType.DB
                    case pt.IonType.DA_ISOLEUCINE | pt.IonType.DA_THREONINE:
                        ion_type = pt.IonType.DA
                    case _:
                        ion_type = ion_info.ion_type

                sequence = None
                if include_annotation:
                    sequence: str | None = pt.parse(frag.sequence).serialize(exclude_charge=True) if frag.sequence is not None else None

                ion = PeptideIon(
                    series=IonSeries(ion_type),
                    position=position,
                    sequence=sequence,
                )
            case pt.IonTypeProperty.INTERNAL:
                if ion_info.id == pt.IonType.IMMONIUM:
                    annot = pt.parse(frag.sequence) if frag.sequence is not None else None
                    if annot is None:
                        raise ValueError("Immonium ion must have a sequence annotation.")

                    if len(annot.sequence) != 1:
                        raise ValueError(f"Immonium ion sequence must be a single amino acid, got {annot.sequence}")

                    # get internal mods at position 0 if they exist
                    mods_str = None
                    if annot.has_internal_mods_at_index(0):
                        internal_mods: pt.Mods[pt.ModificationTags] = annot.get_internal_mods_at_index(0)

                        if len(internal_mods) > 1:
                            raise ValueError(f"Multiple internal mods on immonium ion not supported in mzPAF, got {internal_mods}")
                        if len(internal_mods) == 1 and internal_mods.mods[0].count > 1:
                            raise ValueError(f"Multiple occurrences of internal mod on immonium ion not supported in mzPAF, got {internal_mods}")
                        mods_str = internal_mods.serialize()[1:-1]  # remove surrounding brackets

                    if mods_str == "":
                        raise ValueError(f"Empty modification string for immonium ion is not valid in mzPAF. Internal mods: {internal_mods}")

                    ion = ImmoniumIon(amino_acid=AminoAcids(annot.sequence), modification=mods_str)

                else:
                    if not isinstance(frag.position, tuple) or len(frag.position) != 2:
                        start = -1
                        end = -1
                    else:
                        # Assuming frag.position is Tuple[int, int]
                        start, end = frag.position

                    internal_ion_key = tuple(list(ion_info.ion_type.value))
                    if internal_ion_key not in INTERNAL_MASS_DIFFS:
                        raise ValueError(f"Internal ion type {ion_info.ion_type} not supported in mzPAF.")
                    internal_loss = INTERNAL_MASS_DIFFS[internal_ion_key]

                    sequence = None
                    if include_annotation:
                        sequence: str | None = pt.parse(frag.sequence).serialize(exclude_charge=True) if frag.sequence is not None else None

                    ion = InternalFragment(
                        start_position=start,
                        end_position=end,
                        sequence=sequence,
                    )
            case pt.IonTypeProperty.INTACT:
                if ion_info.ion_type == pt.IonType.PRECURSOR:
                    ion = PrecursorIon()
                else:
                    raise ValueError(f"Cannot convert intact ion type {ion_info.id} to mzPAF.")
            case _:
                pass

    if ion is None:
        raise ValueError(f"Cannot convert fragment with ion type {frag.ion_type} to mzPAF.")

    # handle isotope
    isotopes: list[IsotopeSpecification] = []
    match frag.isotopes:
        case dict():
            _isotopes: Mapping[pt.ElementInfo, int] = frag.isotopes
            if _isotopes is None:
                raise ValueError("Isotopes property is None despite dict type.")
            for element, count in _isotopes.items():
                isotopes.append(IsotopeSpecification(element=str(element), count=count))
        case int() as iso_count:
            isotopes.append(IsotopeSpecification(count=iso_count))
        case None:
            pass
        case _:
            raise TypeError(f"Invalid isotopes type: {type(frag.isotopes)}")

    # handle losses
    losses: list[NeutralLoss] = []
    match frag_losses := frag.losses:
        case dict():
            for loss, count in frag_losses.items():
                if isinstance(loss, float):
                    nloss = NeutralLoss(count=count, base_mass=loss)
                    losses.append(nloss)
                elif isinstance(loss, pt.ChargedFormula):
                    paf_formula = loss.to_mz_paf()
                    sign = paf_formula[0]
                    if sign not in ("+", "-"):
                        raise ValueError(f"Invalid formula sign in loss: {paf_formula}")
                    mult = 1 if sign == "+" else -1
                    nloss = NeutralLoss(count=count * mult, base_formula=paf_formula[1:])  # skip +/- sign
                    losses.append(nloss)
                else:
                    raise TypeError(f"Invalid loss type in tuple: {type(loss)}")
        case None:
            pass
        case _:
            raise TypeError(f"Invalid losses type: {type(frag.losses)}")

    # Add internal loss if applicable
    if internal_loss is not None:
        if internal_loss.startswith("-"):
            cnt = -1
        elif internal_loss.startswith("+"):
            cnt = 1
        else:
            raise ValueError(f"Invalid internal loss format: {internal_loss}")
        losses.append(NeutralLoss(count=cnt, base_formula=internal_loss[1:]))

    pzpaf_adducts: list[Adduct] = []
    if frag._charge_adducts is not None:
        _adducts: tuple[pt.GlobalChargeCarrier, ...] = tuple(mod.value for mod in frag.charge_adducts.mods)
        pzpaf_adducts = [Adduct(count=a.occurance, base_formula=a.to_mz_paf()[2:]) for a in _adducts]

    mz_paf_mass_error = None
    if mass_error is not None:
        mz_paf_mass_error = MassError(value=mass_error, unit=mass_error_type)

    return PafAnnotation(
        ion_type=ion,
        neutral_losses=tuple(losses),
        isotopes=tuple(isotopes),
        adducts=tuple(pzpaf_adducts),
        charge=frag.charge_state,
        mass_error=mz_paf_mass_error,
        confidence=confidence,
    )
