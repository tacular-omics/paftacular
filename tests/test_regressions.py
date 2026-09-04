import pickle

import pytest
from tacular import ELEMENT_LOOKUP

import paftacular as p


def composition_mass(annotation):
    return sum(element.get_mass(True) * count for element, count in annotation.comp().items()) - annotation.charge * 0.000548579909


@pytest.mark.parametrize("text", ["y2{PE}", "b2{PE}", "m2:3{PE}", "y2{PE}[M+H]", "y2{PE}[M+H+Na]^2"])
def test_sequence_composition_mass(text):
    pytest.importorskip("peptacular")
    annotation = p.parse_single(text)
    assert annotation.mass() == pytest.approx(composition_mass(annotation), rel=0, abs=1e-6)


def test_equivalent_protonation():
    assert p.parse_single("y2").mass() == pytest.approx(p.parse_single("y2[M+H]").mass(), rel=0, abs=1e-7)


def test_same_isotope_has_zero_shift():
    assert p.IsotopeSpecification(1, "12C").mass() == 0
    assert not +p.IsotopeSpecification(1, "12C").composition


def test_isotope_substitution_formula():
    annotation = p.parse_single("f{C2H4}+i")
    assert annotation.formula() == "C[13C]H4"
    assert annotation.comp()[ELEMENT_LOOKUP["C"]] == 1
    assert annotation.mass() == pytest.approx(composition_mass(annotation), rel=0, abs=1e-7)


def test_smiles_retains_isotope():
    pytest.importorskip("pysmiles")
    difference = p.parse_single("s{[13CH4]}").mass() - p.parse_single("s{C}").mass()
    assert difference == pytest.approx(1.00335483507, rel=0, abs=1e-9)


@pytest.mark.parametrize("text", ["y2{PE},b2{PT}", "m2:3{PE},y2{PT}", "y2{{Glycan:Hex}PE},b2{PT}"])
def test_sequence_annotation_boundaries(text):
    annotations = p.parse_multi(text)
    assert len(annotations) == 2
    assert annotations[1].sequence == "PT"


@pytest.mark.parametrize("value", [0.000001, 0.123456789, 1e-100, 0.9999999999999999])
def test_numeric_roundtrip(value):
    annotation = p.PafAnnotation.make_peptide("y", 2, confidence=value, mass_error=value)
    assert p.parse_single(annotation.serialize()) == annotation


@pytest.mark.parametrize("text", ["y1^0", "y0", "m4:2", "y1,", "y1,,b2"])
def test_reject_invalid_annotations(text):
    with pytest.raises(ValueError):
        p.parse_single(text)


@pytest.mark.parametrize("component,text", [(p.PeptideIon, "b2!"), (p.Adduct, "+H!"), (p.IsotopeSpecification, "+iN")])
def test_component_consumes_input(component, text):
    with pytest.raises(ValueError):
        component.parse(text)


@pytest.mark.parametrize("text", ["IK", "r[TMT126]", "y2-H2O", "y2[M+H]", "y2+i", "_{sample}", "?4", "p"])
def test_pickle_roundtrip(text):
    annotation = p.parse_single(text)
    assert pickle.loads(pickle.dumps(annotation)) == annotation


@pytest.mark.parametrize("series", ["ax", "bx", "cx", "ay", "by", "cy", "az", "bz", "cz"])
def test_internal_conversion_mass(series):
    pt = pytest.importorskip("peptacular")
    fragment = pt.parse("PEPTIDE").frag(ion_type=series, charge=2, position=(2, 4))
    annotation = p.to_mzpaf(fragment)
    assert annotation.mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)
    assert p.parse_single(annotation.serialize()).mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)


def test_explicit_internal_serialization_mass():
    pytest.importorskip("peptacular")
    ion = p.InternalFragment(2, 3, "PE", p.IonSeries.A, p.IonSeries.Y)
    annotation = p.PafAnnotation(ion)
    assert p.parse_single(annotation.serialize()).mass() == pytest.approx(annotation.mass(), rel=0, abs=1e-6)


def test_reference_cache_bound():
    from paftacular.constants import MAX_CACHE_SIZE

    for index in range(MAX_CACHE_SIZE + 1):
        p.ReferenceIon(f"regression-{index}")
    assert len(p.ReferenceIon._cache) <= MAX_CACHE_SIZE


@pytest.mark.parametrize("series", ["ax", "bx", "cx", "ay", "by", "cy", "az", "bz", "cz"])
def test_internal_component_roundtrip(series):
    ion = p.InternalFragment(2, 3, None, p.IonSeries(series[0]), p.IonSeries(series[1]))
    restored = p.InternalFragment.parse(ion.serialize())
    assert restored.composition == ion.composition
    assert restored.mass() == ion.mass()


def test_invalid_construction_cannot_mutate_cached_component():
    original = p.Adduct(1, "H")
    with pytest.raises(ValueError):
        p.Adduct(True, "H")
    assert type(original.count) is int
    assert p.Adduct(1, "H") is original


@pytest.mark.parametrize("value", [1e-100, 0.123456789123, 100000000.12345679])
def test_neutral_mass_precision(value):
    annotation = p.PafAnnotation.make_peptide("y", 2, neutral_losses=[p.NeutralLoss(-1, base_mass=value)])
    restored = p.parse_single(annotation.serialize())
    assert restored.neutral_losses[0].base_mass == value


def test_smiles_charge_must_be_separate():
    pytest.importorskip("pysmiles")
    with pytest.raises(ValueError, match="neutral molecule"):
        p.parse_single("s{[NH4+]}").mass()


def test_generated_annotation_roundtrips():
    from itertools import product

    for series, position, charge, confidence in product(["a", "b", "c", "x", "y", "z"], [1, 20], [1, 3], [None, 1e-7, 0.123456789]):
        annotation = p.PafAnnotation.make_peptide(series, position, charge=charge, confidence=confidence, neutral_losses=["-H2O"], isotopes=["+i13C"])
        assert p.parse_single(annotation.serialize()) == annotation


def test_component_cannot_swallow_two_sequences():
    with pytest.raises(ValueError):
        p.PeptideIon.parse("y2{PE},b2{PT}")
