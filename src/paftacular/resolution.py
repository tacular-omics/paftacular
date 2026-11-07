from __future__ import annotations

import warnings
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path

from .annotation import PafAnnotation, _require_peptacular, pt
from .comps import InternalFragment, PeptideIon, PrecursorIon
from .constants import IonSeries
from .errors import PaftacularError, reraise_as_paftacular

_PACKAGE_DIR = Path(__file__).parent


def _drop_unassignable_mods(fragment: pt.ProFormaAnnotation, *, warn: bool = True) -> pt.ProFormaAnnotation:
    """Drop the analyte modifications a fragment slice cannot claim, in place.

    peptacular's slice carries an analyte's unlocalized and labile modifications into every
    slice. That is right for a sequence slice and wrong for a fragment: an unlocalized
    modification is assigned to no residue, so no slice can be said to carry it, and a
    fragment loses a labile one. Comparing an embedded sequence against a slice drops them on
    both sides, without the warning, because an embedded sequence is used as written.
    """
    if fragment.has_unknown_mods:
        if warn:
            warnings.warn(
                f"The analyte's unlocalized modification {fragment.unknown_mods_str} is assigned to no residue, so the fragment and its mass exclude it.",
                UserWarning,
                skip_file_prefixes=(str(_PACKAGE_DIR),),
            )
        fragment.clear_unknown_mods()
    if fragment.has_labile_mods:
        fragment.clear_labile_mods()
    return fragment


@reraise_as_paftacular
def resolve(annotation: PafAnnotation, analytes: str | Mapping[int, str]) -> PafAnnotation:
    """Attach the fragment sequence selected from a full ProForma analyte.

    Mapping keys are mzPAF analyte references. An omitted reference selects 1.
    The returned annotation retains context in to_dict(), but not in mzPAF text.
    """
    _require_peptacular()
    ion = annotation.ion_type
    if not isinstance(ion, PeptideIon | InternalFragment | PrecursorIon):
        raise PaftacularError("Analyte resolution supports peptide, internal, and precursor ions")
    reference = annotation.analyte_reference if annotation.analyte_reference is not None else 1
    if isinstance(analytes, str):
        source = analytes
    else:
        try:
            source = analytes[reference]
        except KeyError as error:
            raise PaftacularError(f"Missing analyte reference {reference}") from error
    if not isinstance(source, str) or not source:
        raise PaftacularError("Analyte must be a nonempty ProForma string")
    try:
        analyte = pt.parse(source)
    except ValueError as error:
        raise PaftacularError(f"Invalid analyte {source!r}: {error}") from error
    length = len(analyte.sequence)
    # an embedded sequence is used as written, so there the slice only checks agreement
    warn_dropped = annotation.sequence is None
    if isinstance(ion, PeptideIon):
        if ion.position > length:
            raise PaftacularError(f"Fragment position {ion.position} exceeds analyte length {length}")
        if ion.series in (IonSeries.A, IonSeries.B, IonSeries.C, IonSeries.D, IonSeries.DA, IonSeries.DB):
            fragment = _drop_unassignable_mods(analyte.slice(0, ion.position), warn=warn_dropped)
        else:
            fragment = _drop_unassignable_mods(analyte.slice(length - ion.position, length), warn=warn_dropped)
    elif isinstance(ion, InternalFragment):
        if ion.start_position <= 1 or ion.end_position >= length:
            raise PaftacularError("An internal fragment must exclude both analyte termini")
        fragment = _drop_unassignable_mods(analyte.slice(ion.start_position - 1, ion.end_position), warn=warn_dropped)
    else:
        fragment = analyte
    sequence = fragment.serialize(exclude_charge=True)
    if annotation.sequence is not None:
        try:
            embedded = pt.parse(annotation.sequence)
        except ValueError as error:
            raise PaftacularError(f"Invalid embedded sequence {annotation.sequence!r}: {error}") from error
        if embedded.has_charge:
            raise PaftacularError("Embedded sequence must not specify charge")
        if _drop_unassignable_mods(embedded, warn=False).serialize() != sequence:
            raise PaftacularError("Embedded or resolved sequence disagrees with the selected analyte fragment")
    return replace(annotation, resolved_sequence=annotation.sequence or sequence)
