from collections.abc import Mapping
from dataclasses import replace

from .annotation import PafAnnotation, _require_peptacular, pt
from .comps import InternalFragment, PeptideIon, PrecursorIon
from .constants import IonSeries


def resolve(annotation: PafAnnotation, analytes: str | Mapping[int, str]) -> PafAnnotation:
    """Attach the fragment sequence selected from a full ProForma analyte.

    Mapping keys are mzPAF analyte references. An omitted reference selects 1.
    The returned annotation retains context in to_dict(), but not in mzPAF text.
    """
    _require_peptacular()
    ion = annotation.ion_type
    if not isinstance(ion, PeptideIon | InternalFragment | PrecursorIon):
        raise ValueError("Analyte resolution supports peptide, internal, and precursor ions")
    reference = annotation.analyte_reference if annotation.analyte_reference is not None else 1
    if isinstance(analytes, str):
        source = analytes
    else:
        try:
            source = analytes[reference]
        except KeyError as error:
            raise ValueError(f"Missing analyte reference {reference}") from error
    if not isinstance(source, str) or not source:
        raise ValueError("Analyte must be a nonempty ProForma string")
    analyte = pt.parse(source)
    length = len(analyte.sequence)
    if isinstance(ion, PeptideIon):
        if ion.position > length:
            raise ValueError(f"Fragment position {ion.position} exceeds analyte length {length}")
        if ion.series in (IonSeries.A, IonSeries.B, IonSeries.C):
            fragment = analyte.slice(0, ion.position)
        elif ion.series in (IonSeries.X, IonSeries.Y, IonSeries.Z):
            fragment = analyte.slice(length - ion.position, length)
        else:
            raise ValueError(f"Analyte resolution is not supported for ion series {ion.series}")
    elif isinstance(ion, InternalFragment):
        if ion.start_position <= 1 or ion.end_position >= length:
            raise ValueError("An internal fragment must exclude both analyte termini")
        fragment = analyte.slice(ion.start_position - 1, ion.end_position)
    else:
        fragment = analyte
    sequence = fragment.serialize(exclude_charge=True)
    if annotation.sequence is not None:
        embedded = pt.parse(annotation.sequence)
        if embedded.has_charge:
            raise ValueError("Embedded sequence must not specify charge")
        if embedded.serialize() != sequence:
            raise ValueError("Embedded or resolved sequence disagrees with the selected analyte fragment")
    return replace(annotation, resolved_sequence=annotation.sequence or sequence)
