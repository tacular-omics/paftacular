"""Invariant tests: generated round trips and one check for every entry of each lookup table.

The Hypothesis strategies here extend tests/test_properties.py with negative charges, signed and
electron adducts, every Unimod entry name, explicit internal cleavages and non-canonical
spellings. The table loops make sure each ion series, immonium residue, side-chain substituent,
internal cleavage, reference molecule and isotope that paftacular accepts is exercised at least once.
"""

import math
import warnings
from collections import Counter

import pytest

pytest.importorskip("hypothesis")
pytest.importorskip("peptacular")

from hypothesis import given
from hypothesis import strategies as st
from tacular import AA_LOOKUP, ELEMENT_LOOKUP, REFMOL_LOOKUP, UNIMOD_LOOKUP
from tacular.constants import ELECTRON_MASS, PROTON_MASS

import paftacular as pft
from paftacular import (
    Adduct,
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    IonSeries,
    IsotopeSpecification,
    MassError,
    NamedCompound,
    NeutralLoss,
    PafAnnotation,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)
from paftacular.annotation import _BETA_SUBSTITUENT, _V_ION_RESIDUES
from paftacular.comps.util import formula_to_composition
from paftacular.constants import _INTERNAL_MASS_DIFFS, IMMONIUM_AMINO_ACIDS, InternalSeries
from paftacular.util import format_number

RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
RESIDUE_MODS = [None, None, None, "Oxidation", "Phospho", "Carbamidomethyl", "+15.995", "UNIMOD:35"]
UNIMOD_NAMES = sorted({entry.name for entry in UNIMOD_LOOKUP.values()})
REFERENCE_NAMES = sorted(REFMOL_LOOKUP.keys())
LOSS_FORMULAS = ["H2O", "NH3", "CO", "CO2", "H3PO4", "HPO3", "CH4OS", "C2H3NO", "H", "H2[18O1]", "[2H1]"]
ADDUCT_FORMULAS = ["H", "Na", "K", "NH4", "Li", "[2H2]", "[15N1]H4"]
ISOTOPE_ELEMENTS = [None, "13C", "15N", "18O", "2H", "34S"]
SMILES = ["CN=C=O", "CCO", "c1ccccc1", "OCC(O)CO", "COc(c1)cccc1C#N"]


def _formula_mass(formula: str) -> float:
    """Monoisotopic mass of a plain formula such as ``C2H4N``, from tacular element masses only."""
    return sum(element.get_mass(monoisotopic=True) * count for element, count in formula_to_composition(formula).items())


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, **kwargs)


def _mass_or_error(annotation: PafAnnotation) -> float | str:
    """The mass, or the error type when it cannot be calculated (``+iA``, named compounds)."""
    try:
        return _quiet(annotation.get_mass)
    except ValueError as error:
        return type(error).__name__


# Strategies


def decimals(low: int, high: int, places: int) -> st.SearchStrategy[float]:
    """Floats with a short decimal form, like the numbers written in annotations."""
    return st.integers(low * 10**places, high * 10**places).map(lambda value: value / 10**places)


def residue_strings(length: int) -> st.SearchStrategy[str]:
    residue = st.tuples(st.sampled_from(RESIDUES), st.sampled_from(RESIDUE_MODS))
    body = st.lists(residue, min_size=length, max_size=length).map(lambda items: "".join(aa if mod is None else f"{aa}[{mod}]" for aa, mod in items))
    nterm = st.sampled_from(["", "", "[Acetyl]-"])
    cterm = st.sampled_from(["", "", "-[Amidated]"])
    return st.builds(lambda n, b, c: n + b + c, nterm, body, cterm)


@st.composite
def peptide_ions(draw):
    series = draw(st.sampled_from(list(IonSeries)))
    position = draw(st.integers(1, 30))
    sequence = draw(st.none() if position > 6 else st.none() | residue_strings(position))
    return PeptideIon(series=series, position=position, sequence=sequence)


@st.composite
def internal_fragments(draw, explicit_cleavage: bool = False):
    start = draw(st.integers(1, 20))
    end = draw(st.integers(start, start + 6))
    sequence = draw(st.none() | residue_strings(end - start + 1))
    if not explicit_cleavage:
        return InternalFragment(start, end, sequence=sequence)
    series = draw(st.sampled_from(list(InternalSeries)))
    return InternalFragment(start, end, sequence=sequence, nterm_ion_type=IonSeries(series[0]), cterm_ion_type=IonSeries(series[1]))


immonium_ions = st.builds(
    ImmoniumIon,
    amino_acid=st.sampled_from(sorted(IMMONIUM_AMINO_ACIDS)),
    modification=st.none() | st.sampled_from(UNIMOD_NAMES) | st.sampled_from(["+58.005", "-17.0265", "UNIMOD:4", "MOD:00046"]),
)
reference_ions = st.builds(ReferenceIon, name=st.sampled_from(REFERENCE_NAMES + UNIMOD_NAMES))
formula_ions = st.builds(
    lambda c, h, n, o, label: ChemicalFormula("".join(f"{e}{k}" for e, k in (("C", c), ("H", h), ("N", n), ("O", o)) if k) + label),
    st.integers(1, 20),
    st.integers(1, 30),
    st.integers(0, 4),
    st.integers(0, 6),
    st.sampled_from(["", "", "[13C1]", "[15N2]"]),
)
named_compounds = st.builds(NamedCompound, name=st.from_regex(r"[A-Za-z][A-Za-z0-9 ]{0,10}[A-Za-z0-9]", fullmatch=True))
smiles_ions = st.builds(SMILESCompound, smiles=st.sampled_from(SMILES))
unknown_ions = st.builds(UnknownIon, label=st.none() | st.integers(1, 99))
precursor_ions = st.just(PrecursorIon())

any_ion = st.one_of(
    peptide_ions(),
    internal_fragments(),
    immonium_ions,
    reference_ions,
    formula_ions,
    named_compounds,
    smiles_ions,
    unknown_ions,
    precursor_ions,
)

signed_counts = st.integers(1, 3).flatmap(lambda n: st.sampled_from([n, -n]))
neutral_losses = st.one_of(
    st.builds(lambda count, formula: NeutralLoss(count=count, base_formula=formula), signed_counts, st.sampled_from(LOSS_FORMULAS)),
    st.builds(lambda count, name: NeutralLoss(count=count, base_reference=name), signed_counts, st.sampled_from(REFERENCE_NAMES + UNIMOD_NAMES)),
    st.builds(lambda sign, mass: NeutralLoss(count=sign, base_mass=mass), st.sampled_from([1, -1]), decimals(1, 500, 5).filter(lambda m: m > 0)),
)
isotopes = st.one_of(
    st.builds(lambda count, element: IsotopeSpecification(count=count, element=element), signed_counts, st.sampled_from(ISOTOPE_ELEMENTS)),
    st.builds(lambda count: IsotopeSpecification(count=count, is_average=True), signed_counts),
)
formula_adducts = st.builds(lambda count, formula: Adduct(count=count, base_formula=formula), signed_counts, st.sampled_from(ADDUCT_FORMULAS))
mass_errors = st.builds(MassError, value=decimals(-50, 50, 4), unit=st.sampled_from(["da", "ppm"]))
charges = st.integers(1, 4).flatmap(lambda n: st.sampled_from([n, -n]))
# Carriers outside paftacular's list of known charges, so any charge is allowed with them.
unknown_adducts = st.builds(lambda count, formula: Adduct(count=count, base_formula=formula), signed_counts, st.sampled_from(["Fe", "Zn", "Cu"]))


def _carrier_charge(adducts) -> int:
    return sum(adduct.charge for adduct in adducts)


@st.composite
def adducts_and_charge(draw, ion) -> tuple[tuple[Adduct, ...], int]:
    """A charge with adducts that fit it: none, formula carriers, electrons only, both, or a
    carrier of unknown charge (any charge).

    Sections 4.7 and 4.8: known carriers fix the charge magnitude ([M+2Na]^2), and the written
    charge has no minus sign, so a negative net carrier charge may be written unsigned ([M-2H]^2)
    and is stored negative. The paftacular ^-n form must carry the carriers' own sign.
    """
    if isinstance(ion, ImmoniumIon) and ion.modification is None:
        # IA[M+H] reads as an immonium modification named "M+H" (section 6.1), see test_properties.
        return (), draw(charges)
    kind = draw(st.sampled_from(["none", "formula", "electron", "mixed", "unknown"]))
    if kind == "none":
        return (), draw(charges)
    if kind == "unknown":
        return (draw(unknown_adducts), *draw(st.lists(formula_adducts, max_size=1))), draw(charges)
    if kind == "formula":
        adducts = tuple(draw(st.lists(formula_adducts, min_size=1, max_size=3).filter(lambda a: _carrier_charge(a) != 0)))
    elif kind == "electron":
        # Section 4.4.10: [M-e] is the 1+ ion and [M+2e] the 2- ion, so electrons gained = -charge.
        electrons = -draw(charges)
        split = draw(st.integers(-2, 2).filter(lambda k: k not in (0, electrons)))
        parts = [electrons] if draw(st.booleans()) else [split, electrons - split]
        adducts = tuple(Adduct(count=count, base_formula="e") for count in parts)
    else:
        adducts = draw(
            st.tuples(formula_adducts, st.builds(lambda n: Adduct(count=n, base_formula="e"), signed_counts)).filter(lambda a: _carrier_charge(a) != 0)
        )
    net = _carrier_charge(adducts)
    return adducts, draw(st.sampled_from([net, abs(net)]))


@st.composite
def annotations(draw, ions=any_ion):
    ion = draw(ions)
    adducts, charge = draw(adducts_and_charge(ion))
    return PafAnnotation(
        ion_type=ion,
        analyte_reference=draw(st.none() | st.integers(0, 12)),
        is_auxiliary=draw(st.booleans()),
        neutral_losses=tuple(draw(st.lists(neutral_losses, max_size=3))),
        isotopes=tuple(draw(st.lists(isotopes, max_size=2))),
        adducts=adducts,
        charge=charge,
        mass_error=draw(st.none() | mass_errors),
        confidence=draw(st.none() | decimals(0, 1, 3)),
    )


# Generated round trips


@given(annotations())
def test_round_trip_is_identity(annotation):
    canonical = annotation.serialize()
    parsed = pft.parse(canonical)
    assert parsed == annotation
    assert parsed.serialize() == canonical
    assert pft.parse(parsed.serialize()) == parsed


@given(annotations(ions=internal_fragments(explicit_cleavage=True)))
def test_explicit_internal_cleavage_round_trip(annotation):
    # The cleavage is written as signed gains and losses, so parsing gives a default m{..} ion
    # with those losses. Text and mass survive, and the second parse is identical to the first.
    canonical = annotation.serialize()
    parsed = pft.parse(canonical)
    assert parsed.serialize() == canonical
    assert pft.parse(parsed.serialize()) == parsed
    assert _mass_or_error(parsed) == pytest.approx(_mass_or_error(annotation), rel=0, abs=1e-9)


@given(st.lists(annotations(), min_size=1, max_size=5))
def test_multi_annotation_round_trip(items):
    text = ",".join(item.serialize() for item in items)
    parsed = pft.parse_multi(text)
    assert parsed == items
    assert ",".join(item.serialize() for item in parsed) == text


def _variant(annotation: PafAnnotation, draw) -> str:
    """Write ``annotation`` with optional spellings the grammar allows: explicit counts of one,
    explicit positive signs, a written charge of one and trailing zeros."""

    def count(n: int) -> str:
        sign = "+" if n > 0 else "-"
        magnitude = abs(n)
        return sign + (str(magnitude) if magnitude != 1 or draw(st.booleans()) else "")

    def padded(value: float) -> str:
        text = format_number(value)
        return text + ("0" if "." in text else ".0") if draw(st.booleans()) else text

    parts = ["&" if annotation.is_auxiliary else ""]
    if annotation.analyte_reference is not None:
        parts.append(f"{annotation.analyte_reference}@")
    parts.append(annotation.ion_type.serialize())
    for loss in annotation.neutral_losses:
        if loss.base_mass is not None:
            parts.append(("+" if loss.count > 0 else "-") + f"{loss.base_mass:.{draw(st.integers(5, 8))}f}")
        elif loss.base_reference is not None:
            parts.append(count(loss.count) + f"[{loss.base_reference}]")
        else:
            parts.append(count(loss.count) + str(loss.base_formula))
    for isotope in annotation.isotopes:
        element = "A" if isotope.is_average else (isotope.element or "")
        parts.append(count(isotope.count) + "i" + element)
    if annotation.adducts:
        parts.append("[M" + "".join(count(adduct.count) + adduct.base_formula for adduct in annotation.adducts) + "]")
    charge = annotation.charge
    if charge != 1 or draw(st.booleans()):
        parts.append("^" + ("+" if charge > 0 and draw(st.booleans()) else "") + str(charge))
    if (error := annotation.mass_error) is not None:
        sign = "+" if error.value >= 0 and draw(st.booleans()) else ""
        parts.append("/" + sign + padded(error.value) + ("ppm" if error.unit == "ppm" else ""))
    if annotation.confidence is not None:
        parts.append("*" + padded(annotation.confidence))
    return "".join(parts)


@given(annotations(), st.data())
def test_equivalent_spellings_parse_to_the_canonical_annotation(annotation, data):
    variant = _variant(annotation, data.draw)
    parsed = pft.parse(variant)
    assert parsed == annotation, variant
    assert parsed.serialize() == annotation.serialize()


@given(annotations())
def test_signed_charge_false_drops_only_the_sign(annotation):
    unsigned = annotation.serialize(signed_charge=False)
    magnitude = abs(annotation.charge)
    assert unsigned == annotation.serialize().replace(f"^{annotation.charge}", f"^{magnitude}" if magnitude != 1 else "")


# Every ion series at every charge

SERIES_SEQUENCE = {
    # d keeps residue n, the last one, and w and v keep the first one. Each residue chosen has
    # a section 4.4.3 substituent for the series.
    "d": "PEK",
    "da": "PET",
    "db": "PEI",
    "w": "KPE",
    "wa": "TPE",
    "wb": "IPE",
    "v": "KPE",
}
SIGNED_CHARGES = [1, 2, 3, -1, -2, -3]


def _charged(text: str, charge: int) -> PafAnnotation:
    return pft.parse(f"{text}^{charge}")


def _check_charge_states(text: str, analyte: str | None = None, *, step: float = PROTON_MASS) -> None:
    masses = {}
    for charge in SIGNED_CHARGES:
        annotation = _charged(text, charge)
        if analyte is not None:
            annotation = annotation.resolve(analyte)
        assert annotation.charge == charge
        mz = annotation.mz()
        assert math.isfinite(mz) and mz > 0, (text, charge, mz)
        masses[charge] = annotation.get_mass()
        assert pft.parse(annotation.serialize()) == (annotation if analyte is None else pft.parse(f"{text}^{charge}"))
    for n in (1, 2, 3):
        # A positive charge adds n carriers and a negative charge removes them.
        assert masses[n] - masses[-n] == pytest.approx(2 * n * step, rel=0, abs=1e-9), (text, n)
        assert masses[n] - masses[1] == pytest.approx((n - 1) * step, rel=0, abs=1e-9), (text, n)


def test_charge_state_absolute_anchor():
    # Section 4.4.3: y = sum(AA) + H2O + (H+)z. Hand-typed monoisotopic masses: P 97.05276384,
    # E 129.04259309, K 128.09496302, H2O 18.01056468, proton 1.00727646688.
    neutral = 97.05276384 + 129.04259309 + 128.09496302 + 18.01056468
    for charge in SIGNED_CHARGES:
        expected = (neutral + charge * 1.00727646688) / abs(charge)
        assert pft.parse(f"y3{{PEK}}^{charge}").mz() == pytest.approx(expected, rel=0, abs=1e-6)


@pytest.mark.parametrize("series", list(IonSeries))
def test_every_ion_series_every_charge(series):
    sequence = SERIES_SEQUENCE.get(str(series), "PEK")
    _check_charge_states(f"{series}3{{{sequence}}}")


@pytest.mark.parametrize("series", list(InternalSeries))
def test_every_internal_cleavage_every_charge(series):
    fragment = InternalFragment(2, 4, sequence="EPT", nterm_ion_type=IonSeries(series[0]), cterm_ion_type=IonSeries(series[1]))
    _check_charge_states(PafAnnotation(fragment).serialize())


@pytest.mark.parametrize(
    ("text", "analyte", "step"),
    [
        ("p", "PEPTIDE", PROTON_MASS),
        ("IK", None, PROTON_MASS),
        ("r[TMT126]", None, PROTON_MASS),
        ("s{CCO}", None, PROTON_MASS),
        ("f{C13H9N}", None, -ELECTRON_MASS),
    ],
)
def test_other_ion_types_every_charge(text, analyte, step):
    pytest.importorskip("pysmiles")
    _check_charge_states(text, analyte, step=step)


@pytest.mark.parametrize("charge", SIGNED_CHARGES)
def test_negative_charge_composition_removes_hydrogens(charge):
    positive = pft.parse(f"y3{{PEK}}^{abs(charge)}").comp()
    signed = pft.parse(f"y3{{PEK}}^{charge}").comp()
    hydrogen = ELEMENT_LOOKUP["H"]
    expected = Counter(positive)
    if charge < 0:
        expected[hydrogen] -= 2 * abs(charge)
    assert signed == expected


# Immonium residues

_CO = _formula_mass("CO")


@pytest.mark.parametrize("amino_acid", sorted(IMMONIUM_AMINO_ACIDS))
def test_every_immonium_residue(amino_acid):
    # Section 4.4.5: the immonium ion is the residue less CO, plus the charge carrier.
    annotation = pft.parse(f"I{amino_acid}")
    assert annotation.serialize() == f"I{amino_acid}"
    expected = AA_LOOKUP[amino_acid].monoisotopic_mass - _CO + PROTON_MASS
    assert annotation.mz() == pytest.approx(expected, rel=0, abs=1e-6)
    modified = pft.parse(f"I{amino_acid}[Oxidation]")
    assert modified.mz() - annotation.mz() == pytest.approx(UNIMOD_LOOKUP["Oxidation"].monoisotopic_mass, rel=0, abs=1e-6)
    assert pft.parse(modified.serialize()) == modified


@pytest.mark.parametrize("code", ["B", "X", "Z"])
def test_ambiguous_immonium_residues_rejected(code):
    with pytest.raises(ValueError):
        pft.parse(f"I{code}")


# Side-chain ions (section 4.4.3)

_RESIDUE_MASS = {code: AA_LOOKUP[code].monoisotopic_mass for code in "ACDEFGHIKLMNOPQRSTUVWY"}


# The section 4.4.3 table, typed from the specification: d is sum(n-1 AA) + C2H4N, v is
# sum(c-1 AA) + C2H3NO2 and w is sum(c-1 AA) + C3H4O2, each plus (H+)z.
SPEC_SIDE_CHAIN_OFFSET = {"d": "C2H4N", "v": "C2H3NO2", "w": "C3H4O2"}

# Remnant of residue n for each residue the section 4.4.3 remarks define. The table offsets
# carry one H on the beta carbon. Valine replaces it by CH3, Thr da by OH and db by CH3, Ile
# da by C2H5 and db by CH3. Glycine, alanine and proline have no d or w ion.
_D_GENERIC = "CDEFHKLMNOQRSUWY"
SPEC_SIDE_CHAIN_REMNANT = {
    **{("d", residue): "C2H4N" for residue in _D_GENERIC},
    **{("w", residue): "C3H4O2" for residue in _D_GENERIC},
    ("d", "V"): "C3H6N",
    ("w", "V"): "C4H6O2",
    ("da", "T"): "C2H4NO",
    ("db", "T"): "C3H6N",
    ("da", "I"): "C4H8N",
    ("db", "I"): "C3H6N",
    ("wa", "T"): "C3H4O3",
    ("wb", "T"): "C4H6O2",
    ("wa", "I"): "C5H8O2",
    ("wb", "I"): "C4H6O2",
}


@pytest.mark.parametrize(("series", "formula"), sorted(SPEC_SIDE_CHAIN_OFFSET.items()))
def test_side_chain_offset_table(series, formula):
    # Without a sequence, d, v and w give their section 4.4.3 offset plus a proton.
    assert pft.parse(f"{series}4").get_mass() == pytest.approx(_formula_mass(formula) + PROTON_MASS, rel=0, abs=1e-9)


def test_side_chain_remnants_cover_the_implementation():
    assert set(SPEC_SIDE_CHAIN_REMNANT) == set(_BETA_SUBSTITUENT)


@pytest.mark.parametrize(("series", "residue"), sorted(SPEC_SIDE_CHAIN_REMNANT))
def test_every_side_chain_substituent(series, residue):
    # d_n keeps residues 1..n-1 whole and the remnant of residue n. w_n does the same from the C
    # terminus. Expected values come from the spec remnants above with tacular residue and element
    # masses. The table offsets already hold the terminal groups.
    sequence = f"GA{residue}" if series.startswith("d") else f"{residue}GA"
    kept = sum(_RESIDUE_MASS[aa] for aa in "GA")
    remnant = _formula_mass(SPEC_SIDE_CHAIN_REMNANT[(series, residue)])
    annotation = pft.parse(f"{series}3{{{sequence}}}")
    assert annotation.mz() == pytest.approx(kept + remnant + PROTON_MASS, rel=0, abs=1e-6)


@pytest.mark.parametrize("residue", sorted(_V_ION_RESIDUES))
def test_every_v_ion_residue(residue):
    # v_n keeps residues 2..n whole and the C2H3NO2 table offset in place of residue 1.
    expected = sum(_RESIDUE_MASS[aa] for aa in "GA") + _formula_mass("C2H3NO2") + PROTON_MASS
    assert pft.parse(f"v3{{{residue}GA}}").mz() == pytest.approx(expected, rel=0, abs=1e-6)


@pytest.mark.parametrize("series", ["d", "w", "da", "db", "wa", "wb"])
def test_side_chain_undefined_residues_raise(series):
    defined = {residue for (name, residue) in SPEC_SIDE_CHAIN_REMNANT if name == series}
    for residue in sorted(set(_RESIDUE_MASS) - defined):
        sequence = f"GA{residue}" if series.startswith("d") else f"{residue}GA"
        with pytest.raises(ValueError):
            pft.parse(f"{series}3{{{sequence}}}").get_mass()


# Internal fragments (section 4.4.4)


# The section 4.4.4 table, typed from the specification: rows a/b/c, columns x/y/z.
SPEC_INTERNAL_TABLE = {
    "ax": None,
    "ay": "-CO",
    "az": "-CHNO",
    "bx": "+CO",
    "by": None,
    "bz": "-NH",
    "cx": "+CHNO",
    "cy": "+NH",
    "cz": None,
}


def test_internal_table_covers_the_implementation():
    assert {"".join(pair) for pair in _INTERNAL_MASS_DIFFS} == set(SPEC_INTERNAL_TABLE)


@pytest.mark.parametrize(("pair", "correction"), sorted(SPEC_INTERNAL_TABLE.items()))
def test_internal_specification_table(pair, correction):
    # make_internal(ion_type=) follows the section 4.4.4 table: the by mass plus the listed change.
    base = PafAnnotation.make_internal(2, 4, sequence="EPT")
    annotation = PafAnnotation.make_internal(2, 4, ion_type=pair, sequence="EPT")
    expected = 0.0 if correction is None else (1 if correction[0] == "+" else -1) * _formula_mass(correction[1:])
    assert annotation.get_mass() - base.get_mass() == pytest.approx(expected, rel=0, abs=1e-9)
    assert pft.parse(annotation.serialize()) == annotation


@pytest.mark.parametrize("series", list(InternalSeries))
def test_every_physical_internal_cleavage(series):
    fragment = InternalFragment(2, 4, sequence="EPT", nterm_ion_type=IonSeries(series[0]), cterm_ion_type=IonSeries(series[1]))
    annotation = PafAnnotation(fragment)
    text = annotation.serialize()
    parsed = pft.parse(text)
    assert parsed.serialize() == text
    assert parsed.get_mass() == pytest.approx(annotation.get_mass(), rel=0, abs=1e-9)
    assert parsed.comp() == annotation.comp()
    # by and cz share a composition, so the component parser recovers the cleavage up to that.
    assert InternalFragment.parse(text).composition == fragment.composition


# Reference molecules and isotopes


@pytest.mark.parametrize("name", REFERENCE_NAMES)
def test_every_reference_molecule(name):
    # Section 4.4.7: r[name] is the listed molecule plus the charge carrier.
    annotation = pft.parse(f"r[{name}]")
    assert annotation.serialize() == f"r[{name}]"
    assert annotation.mz() == pytest.approx(REFMOL_LOOKUP[name].monoisotopic_mass + PROTON_MASS, rel=0, abs=1e-6)
    loss = pft.parse(f"y1{{K}}-[{name}]")
    assert pft.parse("y1{K}").get_mass() - loss.get_mass() == pytest.approx(REFMOL_LOOKUP[name].monoisotopic_mass, rel=0, abs=1e-6)


def _isotopes() -> list[tuple[str, float]]:
    cases = []
    for symbol in sorted({str(element) for element, _ in ELEMENT_LOOKUP.keys()}):
        isotopes = ELEMENT_LOOKUP.get_all_isotopes(symbol)
        mono = next((iso for iso in isotopes if iso.is_monoisotopic), None)
        if mono is None:
            continue
        cases.extend((f"{iso.mass_number}{symbol}", iso.mass - mono.mass) for iso in isotopes if not iso.is_monoisotopic)
    return cases


ISOTOPE_CASES = _isotopes()


def test_every_isotope_shift():
    # Section 4.6: +iNX replaces one monoisotopic atom of X by isotope NX.
    base = pft.parse("y1{K}").get_mass()
    failed = []
    for label, shift in ISOTOPE_CASES:
        annotation = pft.parse(f"y1{{K}}+i{label}")
        if annotation.serialize() != f"y1{{K}}+i{label}" or not math.isclose(annotation.get_mass() - base, shift, rel_tol=0, abs_tol=1e-9):
            failed.append(label)
    assert len(ISOTOPE_CASES) > 200
    assert failed == []


@pytest.mark.parametrize("text", ["y5+0i", "y5-0i13C", "y5+00iA", "y5+i+0i"])
def test_zero_isotope_count_rejected(text):
    # Section 4.6: the monoisotopic ion MUST NOT carry an isotope component. +0i used to parse
    # and vanish on serialization.
    with pytest.raises(pft.PafParseError):
        pft.parse(text)


def test_zero_isotope_component_rejected():
    with pytest.raises(ValueError, match="nonzero"):
        IsotopeSpecification(0)
    with pytest.raises(ValueError):
        PafAnnotation.make_peptide("y", 5, isotopes=[0])


# Adduct charge must match the charge state (sections 4.7 and 4.8)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        # One Na+ each: two of them make a 2+ ion, which MUST be written ^2.
        ("y3{PEK}[M+2Na]", r"\[M\+2Na\] carry charge \+2, which does not match charge 1"),
        # One proton cannot make a 2+ ion.
        ("y3{PEK}[M+H]^2", r"\[M\+H\] carry charge \+1, which does not match charge 2"),
        ("y3{PEK}[M+H+Na]", r"carry charge \+2, which does not match charge 1"),
        # An explicit ^-n must agree with the carriers' sign.
        ("y3{PEK}[M+2H]^-2", r"carry charge \+2, which does not match charge -2"),
    ],
)
def test_adduct_charge_mismatch_raises(text, message):
    with pytest.raises(pft.PafParseError, match=message):
        pft.parse(text)


@pytest.mark.parametrize(
    ("text", "charge", "mz"),
    [
        ("y3{PEK}[M+2Na]^2", 2, None),
        ("y3{PEK}[M+H]", 1, "y3{PEK}"),
        ("y3{PEK}[M+H+Na]^2", 2, None),
        ("y3{PEK}[M+2H]^2", 2, "y3{PEK}^2"),
        # Section 4.8: negative-mode charges carry no minus sign, so [M-2H]^2 is the 2- ion.
        ("y3{PEK}[M-2H]^2", -2, "y3{PEK}^-2"),
        ("y3{PEK}[M-H]", -1, "y3{PEK}^-1"),
        ("y3{PEK}[M-2H]^-2", -2, "y3{PEK}^-2"),
        ("y3{PEK}[M+HCOO]", -1, None),
        # A carrier with no known charge is not checked.
        ("y3{PEK}[M+Fe]^3", 3, None),
    ],
)
def test_adduct_charge_matching_charge_parses(text, charge, mz):
    annotation = pft.parse(text)
    assert annotation.charge == charge
    assert pft.parse(annotation.serialize()) == annotation
    assert pft.parse(annotation.serialize(signed_charge=False)) == annotation
    if mz is not None:
        assert annotation.mz() == pytest.approx(pft.parse(mz).mz(), rel=0, abs=1e-9)


def test_adduct_charge_check_applies_to_construction():
    ion = pft.PeptideIon(series=pft.IonSeries.Y, position=3, sequence="PEK")
    with pytest.raises(pft.PaftacularError, match="does not match charge 1"):
        PafAnnotation(ion, adducts=(Adduct(2, "Na"),))
    assert PafAnnotation(ion, adducts=(Adduct(-2, "H"),), charge=2).charge == -2


def test_protonated_two_plus_mz():
    # The real 2+ y3{PEK} ion, which [M+H]^2 used to misreport as 186.60.
    assert pft.parse("y3{PEK}[M+2H]^2").mz() == pytest.approx(187.108, abs=1e-3)
