# Table from the specification showing differences from yb
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
