"""mzPAF 1.0.1 worked examples and independent reference masses.

The fixture tests/reference/mzpaf_reference.json is produced by
tests/reference/generate_reference.py from the specification formulas and pyteomics
atomic data. It does not use paftacular, tacular or peptacular.
"""

import json
from pathlib import Path

import pytest

import paftacular as pft
from paftacular import PafAnnotation, PafParseError

pytest.importorskip("peptacular")

REFERENCE = json.loads((Path(__file__).parent / "reference" / "mzpaf_reference.json").read_text())

# Every annotation string printed in the mzPAF 1.0.1 specification that the grammar accepts.
SPEC_EXAMPLES = [
    "b2-H2O/3.2ppm,b4-H2O^2/3.2ppm",
    "b2-H2O/3.2ppm*0.75,b4-H2O^2/3.2ppm*0.25",
    "1@y12/0.13,2@b9-NH3/0.23",
    "0@y1{K}",
    "0@y1{K}-NH3",
    "y1/-1.4ppm",
    "y1/-0.0002",
    "y4-H2O+2i[M+H+Na]^2",
    "?",
    "?^3",
    "?+2i^4",
    "?17",
    "?17+i/1.45ppm",
    "?17-H2O/-0.87ppm",
    "0@b2{LL}",
    "0@b2{LC[Carbamidomethyl]}",
    "0@b1{[Acetyl]-M}",
    "0@y4{M[Oxidation]ACK}-CH4OS[M+H+Na]^2",
    "m3:6",
    "b3-C2H3NO",
    "b3-[Carbamidomethyl]",
    "m3:6-CO",
    "m3:6-CO-H2O^2",
    "m3:4/1.1ppm,m4:5/1.1ppm",
    "m3:5",
    "IY",
    "IH",
    "IL-CH2",
    "IC[Carbamidomethyl]",
    "IY[Phospho]",
    "IC[+58.005]",
    "p^2",
    "p-H3PO4^2",
    "p^4",
    "p+H^3",
    "p^3",
    "p+2H^2",
    "p+H^2",
    "p+3H",
    "p+2H",
    "p+H",
    "p",
    "r[TMT127N]",
    "r[iTRAQ114]",
    "r[TMT6plex]",
    "r[Hex]",
    "r[Adenine]",
    "r[HexNAc(2)]",
    "0@_{Urocanic Acid}",
    "0@_{Urocanic Acid}+HPO3^2",
    "f{C13H9}/-0.55ppm",
    "f{C12H9N}/0.06ppm",
    "f{C13H9N}/-2.01ppm",
    "f{C13H10N}/-0.11ppm",
    "f{C13H11N}/-0.09ppm",
    "f{C13H12N}/0.26ppm",
    "f{C14H10N}/0.19ppm",
    "f{C14H11N}/0.45ppm",
    "f{C14H10NO}/0.03ppm",
    "f{C16H22O}+i^3",
    "f{C15[13C1]H22O}^3",
    "s{CN=C=O}[M+H]/-0.55ppm",
    "s{COc(c1)cccc1C#N}[M+H+Na]^2/1.29ppm",
    "s{OCCCC=OOH}-H2O[M+H]",
    "s{OCCCC=OOH}[M-H2O+H]",
    "p-[Hex]",
    "y2+CO-H2O",
    "y2-H2O-NH3",
    "p-[TMT6plex]-2H2O-HPO3",
    "p-2[iTRAQ115]",
    "p-[iTRAQ116]-CO-H2O-HPO3",
    "y2-[2H1]-NH3",
    "y5-H2[18O1][M+Na]",
    "y1+i",
    "y1+2i",
    "y1+3i",
    "y1-i",
    "y1-2i",
    "y1+i13C",
    "y1+3i15N",
    "y1+6i13C+2i15N",
    "y1+i15N",
    "y1+2i13C+i15N",
    "y1+iA",
    "y1+2iA",
    "y4[M+Na]",
    "y4[M+NH4]",
    "y4[M+2Na]^2",
    "y4[M+2H+Na]^3",
    "y5-H2O[M+H+Na]^2",
    "y6[M+[2H2]]^2",
    "y5[M+[15N1]H4]",
    "&1@y7/-0.002",
    "&y7/-0.001",
    "&y7/0.001",
    "b6-H2O/-0.005,&y7/0.003",
    "y12/3.4ppm*0.85,b9-NH3/5.2ppm*0.05",
    "1@y7-H2O+i[M+NH4]^2/-0.2ppm*0.5",
    "m5:8-H2O/14.4ppm",
    "p/-1.7ppm",
]

# Spec examples whose numbers are written in a non-canonical way. Serialisation keeps the value.
CANONICAL_FORM = {
    "y7/0.000*0.95": "y7/0*0.95",
    "y12-H2O^2/7.4ppm*0.70": "y12-H2O^2/7.4ppm*0.7",
}


def _serialize(text: str) -> str:
    return ",".join(annotation.serialize() for annotation in pft.parse_multi(text))


@pytest.mark.parametrize("text", SPEC_EXAMPLES)
def test_spec_example_round_trips(text):
    assert _serialize(text) == text


@pytest.mark.parametrize(("text", "canonical"), CANONICAL_FORM.items())
def test_spec_example_canonical_form(text, canonical):
    assert _serialize(text) == canonical
    assert _serialize(canonical) == canonical


@pytest.mark.parametrize("case", REFERENCE, ids=lambda case: f"{case['annotation']}|{case['analyte']}")
def test_reference_mz(case):
    annotation = PafAnnotation.parse(case["annotation"])
    if case["analyte"]:
        annotation = annotation.resolve(case["analyte"])
    assert annotation.mz() == pytest.approx(case["mz"], rel=0, abs=case["tolerance"])
    assert annotation.mass() == pytest.approx(case["mz"] * case["charge"], rel=0, abs=case["tolerance"] * case["charge"])


def test_reference_fixture_covers_every_family():
    notes = {case["note"].split()[0] for case in REFERENCE}
    assert {"4.4.3", "4.4.4", "4.4.5", "4.4.6", "4.4.7", "4.4.9", "4.4.10", "4.5", "4.6", "4.7", "Appendix"} <= notes


def test_spec_stated_hexnac2_reference_mz():
    # Section 4.4.7 states r[HexNAc(2)] is 407.1660 (Unimod HexNAc(2) plus a proton).
    assert PafAnnotation.parse("r[HexNAc(2)]").mz() == pytest.approx(407.1660, abs=5e-5)


def test_spec_stated_generic_isotope_shift():
    # Section 4.6: +i is the 13C minus 12C mass difference, 1.003355 Da.
    base = PafAnnotation.parse("y1{K}").mz()
    assert PafAnnotation.parse("y1{K}+i").mz() - base == pytest.approx(1.003355, abs=1e-6)


@pytest.mark.parametrize(
    ("formula", "printed"),
    [
        ("H", 1.007825),
        ("NH3", 17.026549),
        ("H2O", 18.010565),
        ("CO", 27.994915),
        ("CO2", 43.989829),
        ("HCONH2", 45.021464),
        ("HCOOH", 46.005479),
        ("HPO3", 79.966331),
        ("H3PO4", 97.976896),
    ],
)
def test_spec_neutral_loss_table(formula, printed):
    # The sulfur rows of the section 4.5 table are printed about 1e-5 Da off their formulas
    # and are covered by the reference fixture instead.
    base = PafAnnotation.parse("y1{K}").mass()
    assert base - PafAnnotation.parse(f"y1{{K}}-{formula}").mass() == pytest.approx(printed, abs=1e-6)


def test_formula_ion_rule_c13h9():
    # Section 4.4.9 prints 165.069988 for f{C13H9}, but its own rule (all atoms minus
    # one electron) gives 165.069877. paftacular follows the rule.
    assert PafAnnotation.parse("f{C13H9}").mz() == pytest.approx(165.069877, abs=1e-6)


@pytest.mark.xfail(strict=True, raises=PafParseError, reason="spec 5.2 example puts ^charge before [adduct], which the 6.1 grammar forbids")
def test_spec_example_charge_before_adduct():
    pft.parse_multi("1@y7-H2O+i^2[M+NH4]/-0.2ppm*0.5")


@pytest.mark.xfail(strict=True, raises=PafParseError, reason="spec 4.4.10 electron adducts are outside the 6.1 regex; deferred")
@pytest.mark.parametrize("text", ["y1[M-e]", "s{CN=C=O}[M-e]", "s{CN=C=O}[M+2e]^2"])
def test_spec_electron_adducts(text):
    pft.parse_multi(text)


@pytest.mark.parametrize(
    "text",
    ["r[NotAMolecule]", "p-[NotAMolecule]", "y1{K}-[NotAMolecule]"],
)
def test_unknown_reference_raises_value_error(text):
    with pytest.raises(ValueError, match="NotAMolecule"):
        PafAnnotation.parse(text).mass()


@pytest.mark.parametrize("text", ["d1{G}", "d2{PA}", "w1{P}", "d3{PET}", "w3{IEK}", "da3{PEL}", "wb3{LEK}", "d3{PEK[Acetyl]}"])
def test_side_chain_ion_undefined_residue_raises(text):
    with pytest.raises(ValueError):
        PafAnnotation.parse(text).mass()


def test_side_chain_v_ion_drops_side_chain_modification():
    plain = PafAnnotation.parse("v3{SEK}").mass()
    assert PafAnnotation.parse("v3{S[Phospho]EK}").mass() == pytest.approx(plain, abs=1e-9)


@pytest.mark.parametrize("text", ["da3", "db3", "wa3", "wb3"])
def test_residue_specific_side_chain_ion_needs_sequence(text):
    with pytest.raises(ValueError, match="needs a sequence"):
        PafAnnotation.parse(text).mass()


@pytest.mark.parametrize(("text", "formula"), [("d3", "C2H4N"), ("v3", "C2H3NO2"), ("w3", "C3H4O2"), ("z3", "N-1O")])
def test_offset_only_ion_formula(text, formula):
    assert PafAnnotation.parse(text).ion_type.formula == formula


@pytest.mark.parametrize("text", ["r[Deamidated]", "p-[Deamidated]"])
def test_unimod_composition_change_is_not_a_reference(text):
    with pytest.raises(ValueError, match="not a molecule"):
        PafAnnotation.parse(text).mass()


def test_unimod_reference_formula_and_loss():
    assert PafAnnotation.parse("r[Hex]").formula() == "C6H11O5"
    loss = PafAnnotation.parse("p-[Hex]").neutral_losses[0]
    assert loss.formula == "-C6H10O5"
