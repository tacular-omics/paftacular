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


# mzPAF immonium ions name one amino acid by its one-letter code (section 4.4.5). paftacular
# accepts the 20 standard codes, selenocysteine (U), pyrrolysine (O) and J (I or L, which share
# one mass). B, X and Z have no single mass, so they are rejected.
IMMONIUM_AMINO_ACIDS: frozenset[AminoAcid] = frozenset(AminoAcid(code) for code in "ACDEFGHIJKLMNOPQRSTUVWY")


# A single chemical-formula "atom" token: either a plain element+count (e.g. "H2") or an
# isotope-labeled atom in brackets (e.g. "[18O1]"). Shared by neutral losses and adducts so a
# formula segment can freely mix both forms (e.g. "H2[18O1]", per mzPAF's formula notation).
# The trailing element/count runs are POSSESSIVE (`*+`, py3.11+): [A-Z] overlaps [A-Za-z0-9], so a
# plain `[A-Za-z0-9]*` under the surrounding `_ATOM_TOKEN+` is a classic `(a+)+`-style catastrophic-
# backtracking (ReDoS) shape -- an anchored non-match on a long single-letter run (e.g. "y1+HHHH...!"
# via parse) would hang. Nothing that legitimately follows an atom run starts with an
# alnum char, so refusing to give characters back never rejects a valid annotation.
# The bracketed form is the section 6.2 grammar's isotope atom (ATOM_COUNT): a nucleon count, one
# element symbol and an optional count. A looser form would read a Unimod name that starts with
# a digit (``-[2HPG]``) as a malformed isotope instead of a reference name.
_ATOM_TOKEN = r"(?:\[[0-9]+[A-Z][a-z]?[0-9]*\]|[A-Z][A-Za-z0-9]*+)"

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
# A bracketed name is any reference molecule or Unimod entry name (sections 4.4.5, 4.4.7 and
# 4.5). Unimod names use characters the section 6.1 regex omits (Met->Hse, Hex(1)HexNAc(1),
# Myristoyl+Delta:H(-4), names with spaces or "/"), and some end in a bracketed group
# (Xlink:DSS[156], Cation:Fe[III]). The section 6.2 grammar (BRACE_ENCLOSED_CONTENT) allows one
# nested bracket group, so a name is any run of characters other than brackets and line breaks,
# with balanced one-level bracket groups. A neutral-loss name also keeps its parentheses
# balanced and unnested, as every Unimod name does. The alternatives of each name start with
# different characters, so the repetition cannot backtrack catastrophically.
BRACKETED_NAME = r"(?:[^\[\]\r\n]|\[[^\[\]\r\n]+\])+"
_LOSS_NAME_CHAR = r"[^\[\]()\r\n]"
_LOSS_NAME = rf"(?:{_LOSS_NAME_CHAR}|\({_LOSS_NAME_CHAR}*\)|\[{_LOSS_NAME_CHAR}+\])+"
NEUTRAL_LOSS_REGEX_PATTERN = rf"[+-](?:[0-9]*{_ATOM_TOKEN}+|[0-9]*\[{_LOSS_NAME}\]|[0-9]+(?:\.[0-9]+)?(?!i))"
# An adduct carrier is a formula or an electron, "e" (section 4.4.10: [M-e], [M+2e]).
ELECTRON_CARRIER = "e"
ADDUCT_REGEX_PATTERN = rf"([+-])([0-9]*)({_ATOM_TOKEN}+|{ELECTRON_CARRIER})"


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
_ADDUCT_BODY = rf"M(?:[+-][0-9]*(?:{_ATOM_TOKEN}+|{ELECTRON_CARRIER}))+"
_IMMONIUM = rf"(?:I(?P<immonium>[A-Z])(?:\[(?!{_ADDUCT_BODY}\])(?P<immonium_modification>{BRACKETED_NAME})\])?)"
_REFERENCE = rf"(?P<reference>r(?:(?:\[(?P<reference_label>{BRACKETED_NAME})\])))"
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
