"""The 2.0 API contract: error hierarchy, parse arity, get_mass/mz, negative charge and pickling."""

import copy
import pickle

import pytest

import paftacular as pft
from paftacular import PafAnnotation, PafParseError, PaftacularError, PafUnknownReferenceError

PROTON = 1.007276466812


# Error hierarchy


@pytest.mark.parametrize(
    "call",
    [
        lambda: pft.parse("y5^0"),
        lambda: pft.parse("not an annotation"),
        lambda: pft.parse("b2,y3"),
        lambda: pft.parse(""),
        lambda: pft.parse("y5+iA").get_mass(),
        lambda: pft.parse("y5-17.03").formula(),
        lambda: pft.parse("r[NotAMolecule]").get_mass(),
        lambda: pft.parse("y5").mz(),
        lambda: pft.parse("y2").resolve("P"),
        lambda: pft.parse("y2").resolve({2: "PEPTIDE"}),
        lambda: pft.PeptideIon("q", 2),
        lambda: pft.PeptideIon("y", 0),
        lambda: pft.InternalFragment(4, 2),
        lambda: pft.ImmoniumIon("Z"),
        lambda: pft.NeutralLoss(0, base_formula="H2O"),
        lambda: pft.NeutralLoss(-1),
        lambda: pft.IsotopeSpecification(1, element="C"),
        lambda: pft.Adduct(0, "H"),
        lambda: pft.MassError(1.0, unit="mDa"),  # type: ignore[arg-type]
        lambda: PafAnnotation(pft.PrecursorIon(), charge=0),
        lambda: PafAnnotation(pft.PrecursorIon(), confidence=2.0),
        lambda: PafAnnotation.make_peptide("q", 2),
        lambda: PafAnnotation.make_internal(2, 4, ion_type="qq"),
        lambda: PafAnnotation.from_dict({"schema_version": 99}),
        lambda: pft.ChemicalFormula(""),
    ],
)
def test_user_input_errors_are_paftacular_errors(call):
    with pytest.raises(PaftacularError):
        call()


def test_error_classes_subclass_the_base():
    assert issubclass(PaftacularError, ValueError)
    assert issubclass(PafParseError, PaftacularError)
    assert issubclass(PafUnknownReferenceError, PaftacularError)
    assert issubclass(PafUnknownReferenceError, KeyError)


def test_non_string_input_is_a_type_error():
    with pytest.raises(TypeError):
        pft.parse(5)  # type: ignore[arg-type]


# parse() arity


def test_parse_returns_exactly_one_annotation():
    assert isinstance(pft.parse("y5"), PafAnnotation)
    with pytest.raises(PafParseError, match="parse_multi"):
        pft.parse("y5,b3")
    with pytest.raises(PafParseError, match="got 0"):
        pft.parse("  ")


def test_parse_multi_returns_a_list():
    assert pft.parse_multi("") == []
    assert [str(a) for a in pft.parse_multi("y5, b3")] == ["y5", "b3"]
    assert [str(a) for a in pft.parse_multi("y5")] == ["y5"]


def test_removed_names_are_gone():
    for name in ("parse_single", "parse_batch", "mzPAFParser"):
        assert not hasattr(pft, name)
    for name in ("mass", "dict_composition", "as_dict"):
        assert not hasattr(PafAnnotation, name)


# Peptacular 5 neutral-delta labels


@pytest.mark.parametrize(
    ("text", "delta"),
    [
        ("b3{PEP}-NH3", -17.026549),
        ("b3{PEP}+2HCOOH", 2 * 46.005479),
        ("b3{PEP}-HCONH2", -45.021464),
        ("y2{DE}-H2O-NH3", -18.010565 - 17.026549),
    ],
)
def test_peptacular5_labels_parse(text, delta):
    annotation = pft.parse(text)
    base = pft.parse(text.split("}")[0] + "}")
    assert annotation.get_mass() - base.get_mass() == pytest.approx(delta, rel=0, abs=1e-5)
    assert pft.parse(annotation.serialize()) == annotation


# get_mass / mz


@pytest.mark.parametrize("text", ["y5", "b3^2", "m2:4", "p^2", "y5-H2O"])
def test_mz_needs_a_sequence(text):
    annotation = pft.parse(text)
    annotation.get_mass()  # offset only, still allowed
    with pytest.raises(PaftacularError, match="needs a sequence"):
        annotation.mz()


@pytest.mark.parametrize("text", ["r[TMT126]", "f{C6H12O6}", "IK", "y2{DE}^2"])
def test_mz_without_peptide_context(text):
    annotation = pft.parse(text)
    assert annotation.mz() == pytest.approx(annotation.get_mass() / annotation.charge)


def test_get_mass_is_keyword_only():
    annotation = pft.parse("y2{DE}")
    with pytest.raises(TypeError):
        annotation.get_mass(False)  # type: ignore[misc]
    assert annotation.get_mass(monoisotopic=False) != annotation.get_mass()


# Negative charge


def test_negative_charge_parses_and_round_trips():
    annotation = pft.parse("y2{DE}^-2")
    assert annotation.charge == -2
    assert annotation.serialize() == "y2{DE}^-2"
    assert annotation.serialize(signed_charge=False) == "y2{DE}^2"
    assert pft.parse(annotation.serialize()) == annotation
    assert PafAnnotation.from_dict(annotation.to_dict()) == annotation


def test_negative_charge_mass_removes_protons():
    positive = pft.parse("y2{DE}^2")
    negative = pft.parse("y2{DE}^-2")
    assert positive.get_mass() - negative.get_mass() == pytest.approx(4 * PROTON, rel=0, abs=1e-9)
    assert negative.mz() == pytest.approx(negative.get_mass() / 2)


def test_negative_charge_with_adducts():
    annotation = pft.parse("f{C6H12O6}[M-H]^-1")
    assert annotation.charge == -1
    assert pft.parse(annotation.serialize()) == annotation


# Pickle and copy with keyword-only fields


@pytest.mark.parametrize(
    "text",
    [
        "&2@y5{PEPTI}-H2O+i13C[M+H+Na]^2/-0.55ppm*0.85",
        "m2:4{EPT}",
        "IM[Oxidation]",
        "?42",
        "r[TMT126]-[Adenine]",
        "p^-2",
        "f{C6H12O6}+2iA",
    ],
)
def test_pickle_and_copy_round_trip(text):
    annotation = pft.parse(text)
    assert pickle.loads(pickle.dumps(annotation)) == annotation
    assert copy.deepcopy(annotation) == annotation
    for component in (annotation.ion_type, *annotation.neutral_losses, *annotation.isotopes, *annotation.adducts):
        assert pickle.loads(pickle.dumps(component)) == component


def test_explicit_internal_fields_pickle():
    ion = pft.InternalFragment(2, 4, sequence="EPT", nterm_ion_type=pft.IonSeries.A, cterm_ion_type=pft.IonSeries.X)
    assert pickle.loads(pickle.dumps(ion)) == ion


def test_internal_mass_diffs_is_private():
    assert not hasattr(pft, "INTERNAL_MASS_DIFFS")
    assert pft.PafAnnotation.make_internal(2, 4, ion_type="bx").serialize() == "m2:4+CO"


# Error contract: user-input failures are PaftacularError, never a bare ValueError


@pytest.mark.parametrize(
    "call",
    [
        lambda: pft.parse("y2{PE[Foo]}").get_mass(),
        lambda: pft.parse("IK[+42.010565]").comp(),
        lambda: pft.parse("m1:3{D[Formula:Zz2]IR}").comp(),
        lambda: pft.parse("x6{NVZNAI}").get_mass(),
        lambda: pft.parse("IT[B:Foo]").get_mass(),
        lambda: pft.parse("b2{P[+79.966331]E}").comp(),
        lambda: pft.ChemicalFormula("H²O").get_mass(),
    ],
    ids=["unknown-mod", "immonium-mass-mod", "bad-formula-mod", "bad-residue", "bad-immonium-mod", "mass-mod-comp", "superscript-digit"],
)
def test_calculation_errors_are_paftacular_errors(call):
    with pytest.raises(PaftacularError):
        call()


@pytest.mark.parametrize("text", ["_{foo}", "?", "?42"])
def test_unsupported_calculations_raise_paftacular_error(text):
    annotation = pft.parse(text)
    with pytest.raises(pft.PafUnsupportedCalculationError):
        annotation.get_mass()
    with pytest.raises(pft.PafUnsupportedCalculationError):
        annotation.comp()
    assert issubclass(pft.PafUnsupportedCalculationError, PaftacularError)


@pytest.mark.parametrize("text", ["IK[M+K]", "IK[M+H+Na]^2", "IK[Acetyl][M+K]", "IM[Oxidation]", "IK"])
def test_immonium_adducts_round_trip(text):
    annotation = pft.parse(text)
    assert annotation.serialize() == text
    assert pft.parse(annotation.serialize()) == annotation


def test_immonium_adduct_is_an_adduct_not_a_modification():
    annotation = pft.parse("IK[M+K]")
    assert annotation.ion_type.modification is None
    assert len(annotation.adducts) == 1
    with pytest.raises(PaftacularError):
        pft.ImmoniumIon("K", modification="M+K")


# Global isotope labels, fixed modifications on side-chain ions, charge carriers


def test_global_label_covers_the_ion_offset():
    # a2 of RY: 6 + 9 residue C, less the CO the a ion loses, is 14 labelled C.
    annotation = pft.parse("a2{<13C>RY}")
    assert annotation.formula() == "[13C14]H22N5O2"
    assert annotation.get_mass() == pytest.approx(pft.parse("a2{RY}").get_mass() + 14 * 1.00335483507, rel=0, abs=1e-9)


def test_global_label_covers_formula_deltas():
    labelled = pft.parse("b2{<15N>KE}-NH3").get_mass() - pft.parse("b2{<15N>KE}").get_mass()
    assert labelled == pytest.approx(-(15.00010889888 + 3 * 1.00782503223), rel=0, abs=1e-9)
    # A mass-only delta has no atoms to label.
    assert pft.parse("b2{<15N>KE}-17.026549").get_mass() == pytest.approx(pft.parse("b2{<15N>KE}").get_mass() - 17.026549, rel=0, abs=1e-9)


def test_v_ion_drops_a_global_fixed_modification():
    assert pft.parse("v3{<[Carbamidomethyl]@C>CFQ}").get_mass() == pytest.approx(pft.parse("v3{CFQ}").get_mass(), rel=0, abs=1e-9)
    assert pft.parse("v3{<[Carbamidomethyl]@C>CFC}").get_mass() == pytest.approx(pft.parse("v3{CFC[Carbamidomethyl]}").get_mass(), rel=0, abs=1e-9)
    assert pft.parse("v3{<[Carbamidomethyl]@C>CFQ}").formula() == pft.parse("v3{CFQ}").formula()


@pytest.mark.parametrize("text", ["w3{<[Oxidation]@M>MFQ}", "d3{<[Oxidation]@M>FQM}"])
def test_w_and_d_ions_refuse_a_global_fixed_modification(text):
    with pytest.raises(PaftacularError, match="carries a modification"):
        pft.parse(text).get_mass()


@pytest.mark.parametrize(("adduct", "default"), [("y2{DE}[M+H]", "y2{DE}"), ("y2{DE}[M-H]^-1", "y2{DE}^-1"), ("y2{DE}[M+2H]^2", "y2{DE}^2")])
def test_h_carrier_equals_the_default_charge(adduct, default):
    for monoisotopic in (True, False):
        assert pft.parse(adduct).get_mass(monoisotopic=monoisotopic) == pft.parse(default).get_mass(monoisotopic=monoisotopic)


def test_average_charge_is_natural_hydrogen_less_an_electron():
    from tacular import ELEMENT_LOOKUP
    from tacular.constants import ELECTRON_MASS

    carrier = pft.parse("y2{DE}^2").get_mass(monoisotopic=False) - pft.parse("y2{DE}").get_mass(monoisotopic=False)
    assert carrier == pytest.approx(ELEMENT_LOOKUP["H"].get_mass(monoisotopic=False) - ELECTRON_MASS, rel=0, abs=1e-12)


def test_modification_masses_are_full_precision():
    # Unimod tabulates Oxidation as 15.994915. The exact mass of O is 15.99491461957.
    delta = pft.parse("y2{M[Oxidation]K}").get_mass() - pft.parse("y2{MK}").get_mass()
    assert delta == pytest.approx(15.99491461957, rel=0, abs=1e-10)


@pytest.mark.parametrize("text", ["IK[M+Methyl]", "y2{DE}[M+Methyl]", "y2{DE}-Methyl", "y2{DE}[M+Xx]"])
def test_non_element_formulas_fail_at_parse_time(text):
    with pytest.raises(PafParseError, match="formula"):
        pft.parse(text)


@pytest.mark.parametrize("text", ["y\u0662{DE}", "y2{DE}^\u0662", "y2{DE}-\u0661\u0667.0", "m\u0661:2{DE}"])
def test_non_ascii_digits_are_rejected(text):
    with pytest.raises(PafParseError):
        pft.parse(text)


def test_removed_deuteron_is_labelled():
    # A negative charge removes a hydrogen atom. Under <2H> that atom is a deuteron.
    removed = pft.parse("b2{<2H>PE}^-1")
    assert removed.formula() == "C10[2H13]N2O4"
    assert pft.parse("b2{<2H>PE}[M-H]^-1").get_mass() == pytest.approx(removed.get_mass(), rel=0, abs=1e-12)
    assert pft.parse("b2{<2H>PE}[M-H]^-1").formula() == removed.formula()


def test_h_adduct_of_the_opposite_sign_is_a_hydride():
    from tacular import ELEMENT_LOOKUP
    from tacular.constants import ELECTRON_MASS, PROTON_MASS

    # [M+H]^-1 adds an H atom and an electron. The default +1 charge adds a proton.
    shift = pft.parse("y2{DE}[M+H]^-1").get_mass() - pft.parse("y2{DE}").get_mass()
    assert shift == pytest.approx(ELEMENT_LOOKUP["H"].get_mass() + ELECTRON_MASS - PROTON_MASS, rel=0, abs=1e-12)
