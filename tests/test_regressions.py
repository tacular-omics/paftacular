import pickle

import pytest
from tacular import ELEMENT_LOOKUP

import paftacular as p


def composition_mass(annotation):
    return sum(element.get_mass(monoisotopic=True) * count for element, count in annotation.comp().items()) - annotation.charge * 0.000548579909


@pytest.mark.parametrize("text", ["y2{PE}", "b2{PE}", "m2:3{PE}", "y2{PE}[M+H]", "y2{PE}[M+H+Na]^2"])
def test_sequence_composition_mass(text):
    pytest.importorskip("peptacular")
    annotation = p.parse(text)
    assert annotation.get_mass() == pytest.approx(composition_mass(annotation), rel=0, abs=1e-6)


def test_equivalent_protonation():
    assert p.parse("y2").get_mass() == pytest.approx(p.parse("y2[M+H]").get_mass(), rel=0, abs=1e-7)


def test_same_isotope_has_zero_shift():
    assert p.IsotopeSpecification(1, element="12C").get_mass() == 0
    assert not +p.IsotopeSpecification(1, element="12C").composition


def test_isotope_substitution_formula():
    annotation = p.parse("f{C2H4}+i")
    assert annotation.formula() == "C[13C]H4"
    assert annotation.comp()[ELEMENT_LOOKUP["C"]] == 1
    assert annotation.get_mass() == pytest.approx(composition_mass(annotation), rel=0, abs=1e-7)


def test_smiles_retains_isotope():
    pytest.importorskip("pysmiles")
    difference = p.parse("s{[13CH4]}").get_mass() - p.parse("s{C}").get_mass()
    assert difference == pytest.approx(1.00335483507, rel=0, abs=1e-9)


@pytest.mark.parametrize("text", ["y2{PE},b2{PT}", "m2:3{PE},y2{PT}", "y2{{Glycan:Hex}PE},b2{PT}"])
def test_sequence_annotation_boundaries(text):
    annotations = p.parse_multi(text)
    assert len(annotations) == 2
    assert annotations[1].sequence == "PT"


@pytest.mark.parametrize("value", [0.000001, 0.123456789, 1e-100, 0.9999999999999999])
def test_numeric_roundtrip(value):
    annotation = p.PafAnnotation.make_peptide("y", 2, confidence=value, mass_error=value)
    assert p.parse(annotation.serialize()) == annotation


@pytest.mark.parametrize("text", ["y1^0", "y0", "m4:2", "y1,", "y1,,b2"])
def test_reject_invalid_annotations(text):
    with pytest.raises(ValueError):
        p.parse(text)


@pytest.mark.parametrize("component,text", [(p.PeptideIon, "b2!"), (p.Adduct, "+H!"), (p.IsotopeSpecification, "+iN")])
def test_component_consumes_input(component, text):
    with pytest.raises(ValueError):
        component.parse(text)


@pytest.mark.parametrize("text", ["IK", "r[TMT126]", "y2-H2O", "y2[M+H]", "y2+i", "_{sample}", "?4", "p"])
def test_pickle_roundtrip(text):
    annotation = p.parse(text)
    assert pickle.loads(pickle.dumps(annotation)) == annotation


@pytest.mark.parametrize("series", ["ax", "bx", "cx", "ay", "by", "cy", "az", "bz", "cz"])
def test_internal_conversion_mass(series):
    pt = pytest.importorskip("peptacular")
    fragment = pt.parse("PEPTIDE").frag(ion_type=series, charge=2, position=(2, 4))
    annotation = p.to_mzpaf(fragment)
    assert annotation.get_mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)
    assert p.parse(annotation.serialize()).get_mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)


def test_explicit_internal_serialization_mass():
    pytest.importorskip("peptacular")
    ion = p.InternalFragment(2, 3, sequence="PE", nterm_ion_type=p.IonSeries.A, cterm_ion_type=p.IonSeries.Y)
    annotation = p.PafAnnotation(ion)
    assert p.parse(annotation.serialize()).get_mass() == pytest.approx(annotation.get_mass(), rel=0, abs=1e-6)


def test_reference_cache_bound():
    from paftacular import parser
    from paftacular.constants import MAX_CACHE_SIZE

    for index in range(MAX_CACHE_SIZE + 1):
        p.parse(f"r[regression-{index}]")
    assert len(parser._ION_CACHE) <= MAX_CACHE_SIZE


@pytest.mark.parametrize("series", ["ax", "bx", "cx", "ay", "by", "cy", "az", "bz", "cz"])
def test_internal_component_roundtrip(series):
    ion = p.InternalFragment(2, 3, nterm_ion_type=p.IonSeries(series[0]), cterm_ion_type=p.IonSeries(series[1]))
    restored = p.InternalFragment.parse(ion.serialize())
    assert restored.composition == ion.composition
    assert restored.get_mass() == ion.get_mass()


def test_invalid_construction_raises_paftacular_error():
    with pytest.raises(p.PaftacularError):
        p.Adduct(True, "H")
    assert type(p.Adduct(1, "H").count) is int


@pytest.mark.parametrize("value", [1e-100, 0.123456789123, 100000000.12345679])
def test_neutral_mass_precision(value):
    annotation = p.PafAnnotation.make_peptide("y", 2, neutral_losses=[p.NeutralLoss(-1, base_mass=value)])
    restored = p.parse(annotation.serialize())
    assert restored.neutral_losses[0].base_mass == value


def test_smiles_charge_must_be_separate():
    pytest.importorskip("pysmiles")
    with pytest.raises(ValueError, match="neutral molecule"):
        p.parse("s{[NH4+]}").get_mass()


def test_generated_annotation_roundtrips():
    from itertools import product

    for series, position, charge, confidence in product(["a", "b", "c", "x", "y", "z"], [1, 20], [1, 3], [None, 1e-7, 0.123456789]):
        annotation = p.PafAnnotation.make_peptide(series, position, charge=charge, confidence=confidence, neutral_losses=["-H2O"], isotopes=["+i13C"])
        assert p.parse(annotation.serialize()) == annotation


def test_component_cannot_swallow_two_sequences():
    with pytest.raises(ValueError):
        p.PeptideIon.parse("y2{PE},b2{PT}")


@pytest.mark.parametrize("text", ["y3{PEPTIDE}", "y7{PEP}", "b2{PEP}", "m2:5{PEPTIDE}", "m3:7{PEP}"])
def test_embedded_sequence_length_mismatch_warns(text):
    pytest.importorskip("peptacular")
    annotation = p.parse(text)
    with pytest.warns(UserWarning, match="embedded sequence"):
        annotation.get_mass()
    with pytest.warns(UserWarning, match="embedded sequence"):
        annotation.comp()


@pytest.mark.parametrize("text", ["y3{IDE}", "b2{PE}", "m2:5{EPTI}", "y4{M[Oxidation]ACK}", "b1{[Acetyl]-M}", "y5", "m2:5"])
def test_embedded_sequence_length_match_is_silent(text):
    pytest.importorskip("peptacular")
    import warnings

    annotation = p.parse(text)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        annotation.get_mass()
        annotation.comp()


@pytest.mark.parametrize("text", ["b3", "y3", "z3", "w3", "m2:4", "m2:4-CO", "p", "r[TMT126]"])
def test_ion_composition_is_a_copy(text):
    # PeptideIon, InternalFragment, PrecursorIon and ReferenceIon used to hand out
    # tacular's cached Counter itself, so mutating it changed every later calculation.
    before = p.PeptideIon("y", 3).get_mass()
    ion = p.parse(text).ion_type
    expected = dict(ion.composition)
    ion.composition.clear()
    ion.composition[ELEMENT_LOOKUP["C"]] += 99
    assert dict(ion.composition) == expected
    assert dict(p.parse(text).ion_type.composition) == expected
    assert p.PeptideIon("y", 3).get_mass() == before


@pytest.mark.parametrize("name", ["TMT126-ETD", "TMT131C-ETD", "TMTpro_zero", "sidechain_A"])
@pytest.mark.parametrize("template", ["a1+[{}]", "p-[{}]", "y2-2[{}]^2", "r[{}]-[{}]"])
def test_reference_name_loss_with_hyphen_or_underscore(name, template):
    # mzPAF 1.0.1 section 4.5 allows any reference molecule name as a bracketed loss or gain,
    # and Appendix B names contain "-" (TMT126-ETD) and "_" (TMTpro_zero, sidechain_A). The
    # section 6.2 grammar allows both characters; the loss regex used to reject them.
    text = template.format(name, name)
    annotation = p.parse(text)
    loss = annotation.neutral_losses[-1]
    assert name in str(loss)
    assert annotation.serialize() == text
    reference = p.NeutralLoss.parse(f"+[{name}]").get_mass()
    assert reference == pytest.approx(p.parse(f"r[{name}]").get_mass() - 1.007276466812, rel=0, abs=1e-6)


@pytest.mark.parametrize("name", ["HexNAc(2)", "Hex(1)HexNAc(2)", "dHex(1)Hex(1)"])
@pytest.mark.parametrize("template", ["y2-[{}]", "p+[{}]", "y2-2[{}]^2", "r[{}]-[{}]", "y2-H2O-[{}]+i"])
def test_reference_name_loss_with_parentheses(name, template):
    # mzPAF 1.0.1 section 4.4.7 names Unimod entries such as HexNAc(2), and the section 6.2
    # grammar allows "(" and ")" in bracketed content. The loss regex used to reject them.
    text = template.format(name, name)
    annotation = p.parse(text)
    assert annotation.neutral_losses[-1].base_reference == name
    assert annotation.serialize() == text
    assert p.parse(annotation.serialize()) == annotation
    reference = p.NeutralLoss.parse(f"+[{name}]").get_mass()
    assert reference == pytest.approx(p.parse(f"r[{name}]").get_mass() - 1.007276466812, rel=0, abs=1e-6)


@pytest.mark.parametrize("text", ["y2-[HexNAc(2]", "y2-[HexNAc2)]", "y2-[HexNAc)(2]", "y2-[Hex((2))]"])
def test_reference_name_loss_rejects_unbalanced_parentheses(text):
    with pytest.raises(p.PafParseError):
        p.parse(text)


@pytest.mark.parametrize(
    ("text", "formula", "count"),
    [("y2-H2O", "H2O", -1), ("y2-2H2O", "H2O", -2), ("y2+H2[18O1]", "H2[18O1]", 1), ("y2-NH3-H2O", "H2O", -1)],
)
def test_formula_losses_unchanged_by_reference_name_grammar(text, formula, count):
    loss = p.parse(text).neutral_losses[-1]
    assert (loss.base_formula, loss.count) == (formula, count)


@pytest.mark.parametrize("text", ["r[NotAMolecule]", "p-[NotAMolecule]", "y5-[NotAMolecule]", "y2-[Hex(9)NotAMolecule(1)]"])
def test_unknown_reference_error_is_also_key_error(text):
    # 1.3.2 raised KeyError for an unknown r[...] name. The error is now a ValueError, and it
    # subclasses KeyError too so existing "except KeyError" handlers keep working.
    annotation = p.parse(text)
    with pytest.raises(p.PafUnknownReferenceError, match="NotAMolecule") as info:
        annotation.get_mass()
    assert isinstance(info.value, KeyError)
    assert isinstance(info.value, ValueError)
    assert str(info.value).startswith("Unknown reference molecule '")
    assert "NotAMolecule" in info.value.name
    with pytest.raises(KeyError):
        _ = p.ReferenceIon("NotAMolecule").composition


def test_unknown_reference_error_pickles():
    error = p.PafUnknownReferenceError("NotAMolecule")
    restored = pickle.loads(pickle.dumps(error))
    assert type(restored) is p.PafUnknownReferenceError
    assert (restored.name, str(restored)) == (error.name, str(error))
