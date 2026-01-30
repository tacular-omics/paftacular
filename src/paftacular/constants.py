# Table from the specification showing differences from yb
import re
from enum import StrEnum


class InternalSeries(StrEnum):
    """Enumeration of internal ion series types"""

    AX = "ax"
    BX = "bx"
    CX = "cx"
    AY = "ay"
    BY = "by"
    CY = "cy"
    AZ = "az"
    BZ = "bz"
    CZ = "cz"


INTERNAL_SERIES_TO_DIFF: dict[InternalSeries, str | None] = {
    InternalSeries.AX: None,
    InternalSeries.BX: "+CO",
    InternalSeries.CX: "+CHNO",
    InternalSeries.AY: "-CO",
    InternalSeries.BY: None,
    InternalSeries.CY: "+NH",
    InternalSeries.AZ: "-CHNO",
    InternalSeries.BZ: "-NH",
    InternalSeries.CZ: None,
}

INTERNAL_MASS_DIFFS: dict[tuple[str, str], None | str] = {
    ("a", "x"): None,  #  Default, no difference
    ("b", "x"): "+CO",
    ("c", "x"): "+CHNO",
    ("a", "y"): "-CO",
    ("b", "y"): None,  # Default, no difference
    ("c", "y"): "+NH",
    ("a", "z"): "-CHNO",
    ("b", "z"): "-NH",
    ("c", "z"): None,  # No difference
}


class IonSeries(StrEnum):
    """Enumeration of ion series types"""

    A = "a"
    B = "b"
    C = "c"
    D = "d"
    V = "v"
    W = "w"
    X = "x"
    Y = "y"
    Z = "z"
    DA = "da"
    DB = "db"
    WA = "wa"
    WB = "wb"


class BackboneCleavageType(StrEnum):
    """Types of backbone cleavages for internal fragments"""

    A = "a"  # C-CO bond cleavage
    B = "b"  # CO-NH bond cleavage
    C = "c"  # NH-CH bond cleavage
    X = "x"  # CH-CO bond cleavage
    Y = "y"  # CO-NH bond cleavage
    Z = "z"  # NH-CH bond cleavage


class AnnotationName(StrEnum):
    PRECURSOR = "precursor"
    IMMONIUM = "immonium"
    REFERENCE = "reference"
    NAMED_COMPOUND = "named_compound"
    FORMULA = "formula"
    SMILES = "smiles"
    UNANNOTATED = "unannotated"
    SERIES = "series"
    INTERNAL = "internal"


class SeriesName(StrEnum):
    A = "a"
    B = "b"
    C = "c"
    X = "x"
    Y = "y"
    Z = "z"
    D = "d"
    W = "w"
    V = "v"
    DA = "da"
    DB = "db"
    WA = "wa"
    WB = "wb"


class AminoAcids(StrEnum):
    """Standard amino acids"""

    A = "A"
    C = "C"
    D = "D"
    E = "E"
    F = "F"
    G = "G"
    H = "H"
    I = "I"
    K = "K"
    L = "L"
    M = "M"
    N = "N"
    P = "P"
    Q = "Q"
    R = "R"
    S = "S"
    T = "T"
    V = "V"
    W = "W"
    Y = "Y"


ISOTOPE_REGEX_PATTERN = r"([+-]?)(\d*)i((?:\d+)?(?:[A-Z][a-z]*)?|A)?"
NEUTRAL_LOSS_REGEX_PATTERN = (
    r"[+-](?:\d+(?:\.\d+)?(?!\[)|\d*(?:(?:\[[0-9]+[A-Z][A-Za-z0-9]*\])|(?:[A-Z][A-Za-z0-9]*))+|\d*\[(?:[A-Za-z0-9:\.]+)(?:\[[A-Za-z0-9\.:\-]+\])?\])"
)
ADDUCT_REGEX_PATTERN = r"([+-])(\d*)([A-Z][A-Za-z0-9]*)"


MAX_CACHE_SIZE = 10_000


# Regex components for better readability
_AUXILIARY = r"(?P<is_auxiliary>&)?"
_ANALYTE_REF = r"(?:(?P<analyte_reference>\d+)@)?"

# Ion type patterns
_PEPTIDE_SERIES = r"(?:(?P<series>(?:da|db|wa|wb)|[axbyczdwv]\.?)(?P<ordinal>\d+)(?:\{(?P<sequence_ordinal>.+)\})?)"
_INTERNAL = r"(?P<series_internal>m(?P<internal_start>\d+):(?P<internal_end>\d+)(?:\{(?P<sequence_internal>.+)\})?)"
_PRECURSOR = r"(?P<precursor>p)"
_IMMONIUM = r"(?:I(?P<immonium>[A-Z])(?:\[(?P<immonium_modification>(?:[^\]]+))\])?)"
_REFERENCE = r"(?P<reference>r(?:(?:\[(?P<reference_label>[^\]]+)\])))"
_FORMULA = r"(?:f\{(?P<formula>[A-Za-z0-9\[\]]+)\})"
_NAMED = r"(?:_\{(?P<named_compound>[^\{\}\s,/]+)\})"
_SMILES = r"(?:s\{(?P<smiles>[^\}]+)\})"
_UNKNOWN = r"(?:(?P<unannotated>\?)(?P<unannotated_label>\d+)?)"

# Combine all ion types
_ION_TYPES = f"(?:{_PEPTIDE_SERIES}|{_INTERNAL}|{_PRECURSOR}|{_IMMONIUM}|{_REFERENCE}|{_FORMULA}|{_NAMED}|{_SMILES}|{_UNKNOWN})"

# Modifiers
_NEUTRAL_LOSSES = r"(?P<neutral_losses>(?:[+-](?:\d+(?:\.\d+)?|\d*(?:(?:(?:\[[0-9]+[A-Z][A-Za-z0-9]*\])\
    |(?:[A-Z][A-Za-z0-9]*))+)|(?:\d*\[(?:(?:[A-Za-z0-9:\.]+)(?:\[(?:[A-Za-z0-9\.:\-]+)\])?)\])))+)?"
_ISOTOPE = r"(?P<isotope>(?:(?:[+-]\d*)i(?:(?:\d+)?(?:[A-Z][a-z]*)?|A)?)+)?"
_ADDUCTS = r"(?:\[(?P<adducts>M(?:[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?"
_CHARGE = r"(?:\^(?P<charge>[+-]?\d+))?"
_MASS_ERROR = r"(?:/(?P<mass_error>-?\d+(?:\.\d+)?)(?P<mass_error_unit>ppm)?)?"
_CONFIDENCE = r"(?:\*(?P<confidence>\d*(?:\.\d+)?))?"

# Full pattern
FULL_PAF_PATTERN = re.compile(f"^{_AUXILIARY}{_ANALYTE_REF}{_ION_TYPES}{_NEUTRAL_LOSSES}{_ISOTOPE}{_ADDUCTS}{_CHARGE}{_MASS_ERROR}{_CONFIDENCE}$")
