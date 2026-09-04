import json

import pytest

pt = pytest.importorskip("peptacular")

import paftacular as p  # noqa: E402


@pytest.mark.parametrize("text,sequence", [("b2", "PE"), ("y2", "DE"), ("m2:4", "EPT"), ("p^2", "PEPTIDE")])
def test_resolve_sequence_and_mass(text, sequence):
    original = p.parse_single(text)
    resolved = original.resolve("PEPTIDE/3")
    assert resolved.sequence == sequence
    assert original.sequence is None
    assert resolved.serialize() == text
    expected = pt.parse(sequence).mass(ion_type=resolved.peptacular_ion_type, charge=resolved.charge)
    assert resolved.mass() == pytest.approx(expected, rel=0, abs=1e-6)
    assert resolved.mz() == pytest.approx(expected / resolved.charge, rel=0, abs=1e-6)
    assert p.PafAnnotation.from_dict(json.loads(json.dumps(resolved.to_dict()))) == resolved


def test_resolve_retains_modifications_and_isotopes():
    resolved = p.parse_single("b2+i").resolve("PE[Oxidation]PTIDE")
    expected = pt.parse("PE[Oxidation]").mass(ion_type="b", charge=1) + 1.00335483507
    assert resolved.mass() == pytest.approx(expected, rel=0, abs=1e-6)
    assert resolved.mass() == pytest.approx(sum(e.get_mass(True) * n for e, n in resolved.comp().items()) - 0.000548579909, rel=0, abs=1e-6)
    assert resolved.formula()


def test_resolve_mapping_and_embedded_match():
    assert p.resolve(p.parse_single("2@y2{DE}"), {1: "AAAA", 2: "PEPTIDE"}).sequence == "DE"
    assert p.parse_single("b2").resolve({1: "PEPTIDE"}).sequence == "PE"


@pytest.mark.parametrize(
    "text,context",
    [
        ("2@y2", {1: "PEPTIDE"}),
        ("y8", "PEPTIDE"),
        ("y2{PE}", "PEPTIDE"),
        ("m1:3", "PEPTIDE"),
        ("m2:7", "PEPTIDE"),
        ("?", "PEPTIDE"),
        ("d2", "PEPTIDE"),
        ("p", ""),
        ("b2{PE/2}", "PEPTIDE"),
    ],
)
def test_resolution_errors(text, context):
    with pytest.raises(ValueError):
        p.parse_single(text).resolve(context)
