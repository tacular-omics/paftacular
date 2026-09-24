"""Property tests for mzPAF annotations: round trip, charge consistency and parse errors."""

import dataclasses
import warnings

import pytest

pytest.importorskip("hypothesis")

from hypothesis import assume, example, given
from hypothesis import strategies as st
from tacular import REFMOL_LOOKUP

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
    PafParseError,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)

pytest.importorskip("peptacular")

PROTON = 1.007276466812
ELECTRON = 0.000548579909
C13_SHIFT = 1.00335483507

RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
RESIDUE_MODS = [None, None, None, "Oxidation", "Phospho", "Carbamidomethyl", "+15.995"]
LOSS_FORMULAS = ["H2O", "NH3", "CO", "CO2", "H3PO4", "HPO3", "CH4OS", "C2H3NO", "H"]
REFERENCE_NAMES = sorted(REFMOL_LOOKUP.keys()) if hasattr(REFMOL_LOOKUP, "keys") else ["TMT126", "iTRAQ114"]
UNIMOD_NAMES = ["Hex", "HexNAc", "Phospho", "Carbamidomethyl"]
ADDUCT_FORMULAS = ["H", "Na", "K", "NH4", "Li"]
ISOTOPE_ELEMENTS = [None, "13C", "15N", "18O", "2H", "34S"]
SMILES = ["CN=C=O", "CCO", "c1ccccc1", "OCC(O)CO"]


def decimals(low: int, high: int, places: int) -> st.SearchStrategy[float]:
    """Floats with a short decimal form, like the numbers written in annotations."""
    return st.integers(low * 10**places, high * 10**places).map(lambda value: value / 10**places)


def residue_strings(length: int) -> st.SearchStrategy[str]:
    residue = st.tuples(st.sampled_from(RESIDUES), st.sampled_from(RESIDUE_MODS))
    return st.lists(residue, min_size=length, max_size=length).map(lambda items: "".join(aa if mod is None else f"{aa}[{mod}]" for aa, mod in items))


@st.composite
def peptide_ions(draw, with_sequence: bool | None = None):
    series = draw(st.sampled_from([s for s in IonSeries]))
    position = draw(st.integers(1, 6))
    has_sequence = draw(st.booleans()) if with_sequence is None else with_sequence
    sequence = draw(residue_strings(position)) if has_sequence else None
    return PeptideIon(series=series, position=position, sequence=sequence)


@st.composite
def internal_fragments(draw, with_sequence: bool | None = None):
    start = draw(st.integers(2, 12))
    end = draw(st.integers(start, start + 5))
    has_sequence = draw(st.booleans()) if with_sequence is None else with_sequence
    sequence = draw(residue_strings(end - start + 1)) if has_sequence else None
    return InternalFragment(start_position=start, end_position=end, sequence=sequence)


immonium_ions = st.builds(
    ImmoniumIon,
    amino_acid=st.sampled_from(RESIDUES),
    modification=st.sampled_from([None, "Oxidation", "Phospho", "Carbamidomethyl", "+58.005", "-17.0265"]),
)
reference_ions = st.builds(ReferenceIon, name=st.sampled_from(REFERENCE_NAMES + UNIMOD_NAMES))
formula_ions = st.builds(
    lambda c, h, n, o: ChemicalFormula("".join(f"{e}{k}" for e, k in (("C", c), ("H", h), ("N", n), ("O", o)) if k)),
    st.integers(1, 20),
    st.integers(1, 30),
    st.integers(0, 4),
    st.integers(0, 6),
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
    st.builds(lambda count, name: NeutralLoss(count=count, base_reference=name), signed_counts, st.sampled_from(REFERENCE_NAMES[:5] + UNIMOD_NAMES)),
    st.builds(lambda sign, mass: NeutralLoss(count=sign, base_mass=mass), st.sampled_from([1, -1]), decimals(1, 500, 5).filter(lambda m: m > 0)),
)
isotopes = st.one_of(
    st.builds(lambda count, element: IsotopeSpecification(count=count, element=element), signed_counts, st.sampled_from(ISOTOPE_ELEMENTS)),
    st.builds(lambda count: IsotopeSpecification(count=count, is_average=True), signed_counts),
)
adducts = st.builds(lambda count, formula: Adduct(count=count, base_formula=formula), st.integers(1, 3), st.sampled_from(ADDUCT_FORMULAS))
mass_errors = st.builds(MassError, value=decimals(-50, 50, 4), unit=st.sampled_from(["da", "ppm"]))


def _adducts_allowed(ion) -> bool:
    # IA[M+H] reads as an immonium modification named "M+H" under both the section 6.1
    # regex and the lark grammar, so an unmodified immonium ion cannot carry adducts.
    # See test_unmodified_immonium_with_adduct_round_trip.
    return not (isinstance(ion, ImmoniumIon) and ion.modification is None)


@st.composite
def annotations(draw, ions=any_ion, with_adducts: bool = True):
    ion = draw(ions)
    return PafAnnotation(
        ion_type=ion,
        analyte_reference=draw(st.none() | st.integers(0, 9)),
        is_auxiliary=draw(st.booleans()),
        neutral_losses=tuple(draw(st.lists(neutral_losses, max_size=3))),
        isotopes=tuple(draw(st.lists(isotopes, max_size=2))),
        adducts=tuple(draw(st.lists(adducts, max_size=2))) if with_adducts and _adducts_allowed(ion) else (),
        charge=draw(st.integers(1, 4)),
        mass_error=draw(st.none() | mass_errors),
        confidence=draw(st.none() | decimals(0, 1, 3)),
    )


def _quiet(function, *args):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args)


# Round trip


@given(annotations())
def test_parse_serialize_round_trip(annotation):
    text = annotation.serialize()
    parsed = pft.parse_single(text)
    assert parsed == annotation
    assert parsed.serialize() == text


@given(st.lists(annotations(), min_size=1, max_size=4))
def test_multi_annotation_round_trip(items):
    text = ",".join(item.serialize() for item in items)
    assert pft.parse_multi(text) == items


@pytest.mark.xfail(strict=True, reason="spec ambiguity: [M+H] after an unmodified immonium ion parses as its modification")
def test_unmodified_immonium_with_adduct_round_trip():
    annotation = PafAnnotation(ion_type=ImmoniumIon("A"), adducts=(Adduct(count=1, base_formula="H"),))
    assert pft.parse_single(annotation.serialize()) == annotation


@given(annotations())
def test_dict_round_trip(annotation):
    assert PafAnnotation.from_dict(annotation.to_dict()) == annotation


# Mass and m/z consistency

computable_ions = st.one_of(
    peptide_ions(with_sequence=True),
    internal_fragments(with_sequence=True),
    immonium_ions,
    reference_ions,
    formula_ions,
)
formula_losses = st.builds(lambda count, formula: NeutralLoss(count=count, base_formula=formula), signed_counts, st.sampled_from(LOSS_FORMULAS))


def _mass_or_skip(annotation):
    try:
        return _quiet(annotation.mass)
    except ValueError:
        # Side-chain ions on residues that have none, for example d on G. Covered by explicit tests.
        assume(False)


@given(computable_ions, st.lists(formula_losses, max_size=2), st.integers(1, 5))
def test_charge_adds_one_proton_per_charge(ion, losses, charge):
    single = PafAnnotation(ion_type=ion, neutral_losses=tuple(losses))
    charged = dataclasses.replace(single, charge=charge)
    base = _mass_or_skip(single)
    mass = _quiet(charged.mass)
    step = -ELECTRON if isinstance(ion, ChemicalFormula) else PROTON
    assert mass - base == pytest.approx((charge - 1) * step, rel=0, abs=1e-9)
    assert _quiet(charged.mz) * charge == pytest.approx(mass, rel=1e-12, abs=1e-9)


@given(computable_ions, st.lists(adducts, min_size=1, max_size=2), st.integers(1, 5))
def test_adduct_mass_independent_of_charge_except_electrons(ion, ion_adducts, charge):
    assume(not isinstance(ion, ChemicalFormula))
    single = PafAnnotation(ion_type=ion, adducts=tuple(ion_adducts))
    charged = dataclasses.replace(single, charge=charge)
    base = _mass_or_skip(single)
    assert _quiet(charged.mass) - base == pytest.approx(-(charge - 1) * ELECTRON, rel=0, abs=1e-9)


@given(computable_ions, st.integers(-3, 3).filter(bool), st.integers(1, 4))
def test_generic_isotope_shifts_mz(ion, count, charge):
    plain = PafAnnotation(ion_type=ion, charge=charge)
    shifted = dataclasses.replace(plain, isotopes=(IsotopeSpecification(count=count),))
    base = _mass_or_skip(plain)
    assert _quiet(shifted.mz) - base / charge == pytest.approx(count * C13_SHIFT / charge, rel=0, abs=1e-9)


@given(st.one_of(peptide_ions(with_sequence=True), internal_fragments(with_sequence=True), reference_ions, formula_ions), st.integers(1, 4))
def test_mass_matches_composition(ion, charge):
    sequence = getattr(ion, "sequence", None) or ""
    assume("+" not in sequence)  # a mass-only modification has no composition
    annotation = PafAnnotation(ion_type=ion, charge=charge)
    mass = _mass_or_skip(annotation)
    composition = _quiet(annotation.comp)
    composition_mass = sum(element.get_mass(True) * count for element, count in composition.items())
    # Reference masses and Unimod masses are stored to 6 or more decimals. Each named modification adds up to 5e-7.
    tolerance = 1e-6 + 5e-7 * sequence.count("[")
    assert mass == pytest.approx(composition_mass - charge * ELECTRON, rel=0, abs=tolerance)


# Invalid input raises clean errors

MZPAF_ALPHABET = "abcdmpxyzvwIrf_s?&@{}[]()^/*+-:,.0123456789HCNOPSKMeiAppm "


@example("")
@example("y")
@example("b0")
@example("y1^0")
@example("m3:2")
@example("r[]")
@example("y1[M+]")
@example("y1*2")
@example("y1/abc")
@example("f{}")
@example("s{}")
@example("_{}")
@example("y1-")
@example("0@")
@given(st.text(alphabet=MZPAF_ALPHABET, max_size=25))
def test_arbitrary_text_parses_or_raises_parse_error(text):
    try:
        result = pft.parse_multi(text)
    except (PafParseError, ValueError):
        return
    assert result or not text.strip()  # blank input is an empty annotation list by design


@given(annotations(), st.data())
def test_mutated_annotation_parses_or_raises_parse_error(annotation, data):
    text = annotation.serialize()
    index = data.draw(st.integers(0, len(text)))
    action = data.draw(st.sampled_from(["insert", "delete", "replace"]))
    char = data.draw(st.sampled_from(MZPAF_ALPHABET))
    if action == "insert":
        text = text[:index] + char + text[index:]
    elif action == "delete":
        text = text[:index] + text[index + 1 :]
    else:
        text = text[:index] + char + text[index + 1 :]
    try:
        pft.parse_multi(text)
    except (PafParseError, ValueError):
        pass


@given(annotations())
def test_mass_raises_only_documented_errors(annotation):
    try:
        _quiet(annotation.mz)
    except (ValueError, NotImplementedError, ImportError):
        pass
