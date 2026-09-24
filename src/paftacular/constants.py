"""Enums, grammar regexes and the internal-fragment correction table."""

import re
from enum import StrEnum
from types import MappingProxyType

from tacular import AminoAcid


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


# Table from the specification (section 4.4.4) showing differences from by.
_INTERNAL_SERIES_TO_DIFF: MappingProxyType[InternalSeries, str | None] = MappingProxyType(
    {
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
)

_INTERNAL_MASS_DIFFS: MappingProxyType[tuple[str, str], None | str] = MappingProxyType(
    {
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
)


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


# mzPAF immonium ions use the 20 standard amino acids. tacular.AminoAcid also lists
# ambiguous (B, J, X, Z) and rare (O, U) codes, which paftacular does not accept here.
IMMONIUM_AMINO_ACIDS: frozenset[AminoAcid] = frozenset(AminoAcid(code) for code in "ACDEFGHIKLMNPQRSTVWY")


# A single chemical-formula "atom" token: either a plain element+count (e.g. "H2") or an
# isotope-labeled atom in brackets (e.g. "[18O1]"). Shared by neutral losses and adducts so a
# formula segment can freely mix both forms (e.g. "H2[18O1]", per mzPAF's formula notation).
# The trailing element/count runs are POSSESSIVE (`*+`, py3.11+): [A-Z] overlaps [A-Za-z0-9], so a
# plain `[A-Za-z0-9]*` under the surrounding `_ATOM_TOKEN+` is a classic `(a+)+`-style catastrophic-
# backtracking (ReDoS) shape -- an anchored non-match on a long single-letter run (e.g. "y1+HHHH...!"
# via parse) would hang. Nothing that legitimately follows an atom run starts with an
# alnum char, so refusing to give characters back never rejects a valid annotation.
_ATOM_TOKEN = r"(?:\[[0-9]+[A-Z][A-Za-z0-9]*+\]|[A-Z][A-Za-z0-9]*+)"

# Isotope-nucleon-count is mandatory once an element is specified (mzPAF: "+iN" with no count is
# invalid); at most one lowercase letter follows the element symbol (real element symbols are 1-2
# letters). A bare "i" with nothing after it (generic isotope, no element) remains valid.
_ISOTOPE_ELEMENT = r"(?:(?:[0-9]+[A-Z][a-z]?)|A)"
ISOTOPE_REGEX_PATTERN = rf"([+-]?)([0-9]*)i({_ISOTOPE_ELEMENT})?"

# A single signed neutral-loss/gain token. Order matters: try "count? + formula" (which may embed
# isotope-bracket atoms) and "count? + [reference name]" before the bare-mass fallback, so a
# count-prefixed formula like "-2H2O" isn't misread as a bare mass of "-2" with "H2O" dropped. The
# bare-mass alternative also excludes being followed by "i" so e.g. "+2i13C" is left whole for the
# isotope component instead of being split into a bare-mass loss of "+2" plus a dangling "i13C".
# Bracketed names also allow "_" and "-": section 4.5 permits any reference molecule name there
# (Appendix B has TMTpro_zero, sidechain_A, TMT126-ETD) and the section 6.2 grammar allows both,
# although the section 6.1 regex omits them. They also allow balanced, unnested parentheses
# for Unimod names such as HexNAc(2) (section 4.4.7). The two alternatives of _REFERENCE_NAME
# start with different characters, so the repetition cannot backtrack catastrophically.
_REFERENCE_NAME_CHAR = r"[A-Za-z0-9:\._\-]"
_REFERENCE_NAME = rf"(?:{_REFERENCE_NAME_CHAR}|\({_REFERENCE_NAME_CHAR}*\))+"
NEUTRAL_LOSS_REGEX_PATTERN = rf"[+-](?:[0-9]*{_ATOM_TOKEN}+|[0-9]*\[{_REFERENCE_NAME}(?:\[[A-Za-z0-9\.:\-]+\])?\]|[0-9]+(?:\.[0-9]+)?(?!i))"
ADDUCT_REGEX_PATTERN = rf"([+-])([0-9]*)({_ATOM_TOKEN}+)"


# Bound for the parser's component caches (keyed by annotation substring).
MAX_CACHE_SIZE = 10_000


# Regex components for better readability
_AUXILIARY = r"(?P<is_auxiliary>&)?"
_ANALYTE_REF = r"(?:(?P<analyte_reference>[0-9]+)@)?"

# Ion type patterns
_PEPTIDE_SERIES = r"(?:(?P<series>(?:da|db|wa|wb)|[axbyczdwv]\.?)(?P<ordinal>[0-9]+)(?:\{(?P<sequence_ordinal>.+)\})?)"
_INTERNAL = r"(?P<series_internal>m(?P<internal_start>[0-9]+):(?P<internal_end>[0-9]+)(?:\{(?P<sequence_internal>.+)\})?)"
_PRECURSOR = r"(?P<precursor>p)"
# Adduct text inside the brackets, ``M+H+Na``. An immonium modification never matches it, so
# ``IK[M+K]`` is the K immonium ion with a K+ adduct.
_ADDUCT_BODY = rf"M(?:[+-][0-9]*{_ATOM_TOKEN}+)+"
_IMMONIUM = rf"(?:I(?P<immonium>[A-Z])(?:\[(?!{_ADDUCT_BODY}\])(?P<immonium_modification>(?:[^\]]+))\])?)"
_REFERENCE = r"(?P<reference>r(?:(?:\[(?P<reference_label>[^\]]+)\])))"
_FORMULA = r"(?:f\{(?P<formula>[A-Za-z0-9\[\]]+)\})"
_NAMED = r"(?:_\{(?P<named_compound>[^\{\}/]+)\})"
_SMILES = r"(?:s\{(?P<smiles>[^\}]+)\})"
_UNKNOWN = r"(?:(?P<unannotated>\?)(?P<unannotated_label>[0-9]+)?)"

# Combine all ion types
_ION_TYPES = f"(?P<ion>{_PEPTIDE_SERIES}|{_INTERNAL}|{_PRECURSOR}|{_IMMONIUM}|{_REFERENCE}|{_FORMULA}|{_NAMED}|{_SMILES}|{_UNKNOWN})"

# Modifiers
_NEUTRAL_LOSSES = rf"(?P<neutral_losses>(?:{NEUTRAL_LOSS_REGEX_PATTERN})+)?"
_ISOTOPE = rf"(?P<isotope>(?:(?:[+-][0-9]*)i(?:{_ISOTOPE_ELEMENT})?)+)?"
_ADDUCTS = rf"(?:\[(?P<adducts>{_ADDUCT_BODY})\])?"
_CHARGE = r"(?:\^(?P<charge>[+-]?[0-9]+))?"
_MASS_ERROR = r"(?:/(?P<mass_error>[+-]?[0-9]+(?:\.[0-9]+)?)(?P<mass_error_unit>ppm)?)?"
_CONFIDENCE = r"(?:\*(?P<confidence>[0-9]*(?:\.[0-9]+)?))?"

# Full pattern
_ANNOTATION_PATTERN_BODY = f"{_AUXILIARY}{_ANALYTE_REF}{_ION_TYPES}{_NEUTRAL_LOSSES}{_ISOTOPE}{_ADDUCTS}{_CHARGE}{_MASS_ERROR}{_CONFIDENCE}"
FULL_PAF_PATTERN = re.compile(f"^{_ANNOTATION_PATTERN_BODY}$")
# Same pattern with no end anchor, for matching one annotation out of a comma-separated string
# (mzPAF spec Appendix A: greedily match one annotation, then require a comma or end-of-string).
PARTIAL_PAF_PATTERN = re.compile(_ANNOTATION_PATTERN_BODY)
