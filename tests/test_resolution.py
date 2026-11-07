import json
import warnings

import pytest

pt = pytest.importorskip("peptacular")

import paftacular as p  # noqa: E402


@pytest.mark.parametrize("text,sequence", [("b2", "PE"), ("y2", "DE"), ("z2", "DE"), ("m2:4", "EPT"), ("p^2", "PEPTIDE")])
def test_resolve_sequence_and_mass(text, sequence):
    original = p.parse(text)
    resolved = original.resolve("PEPTIDE/3")
    assert resolved.sequence == sequence
    assert original.sequence is None
    assert resolved.serialize() == text
    expected = pt.parse(sequence).mass(ion_type=resolved.peptacular_ion_type, charge=resolved.charge)
    assert resolved.get_mass() == pytest.approx(expected, rel=0, abs=1e-6)
    assert resolved.mz() == pytest.approx(expected / resolved.charge, rel=0, abs=1e-6)
    assert p.PafAnnotation.from_dict(json.loads(json.dumps(resolved.to_dict()))) == resolved


def test_resolve_retains_modifications_and_isotopes():
    resolved = p.parse("b2+i").resolve("PE[Oxidation]PTIDE")
    expected = pt.parse("PE[Oxidation]").mass(ion_type="b", charge=1) + 1.00335483507
    assert resolved.get_mass() == pytest.approx(expected, rel=0, abs=1e-6)
    assert resolved.get_mass() == pytest.approx(sum(e.get_mass(monoisotopic=True) * n for e, n in resolved.comp().items()) - 0.000548579909, rel=0, abs=1e-6)
    assert resolved.formula()


@pytest.mark.parametrize(("text", "sequence"), [("d2", "PE"), ("v3", "IDE"), ("w2", "DE"), ("wa3", "IDE")])
def test_resolve_side_chain_ions(text, sequence):
    resolved = p.parse(text).resolve("PEPTIDE")
    assert resolved.sequence == sequence
    embedded = p.parse(f"{text}{{{sequence}}}")
    assert resolved.get_mass() == pytest.approx(embedded.get_mass(), rel=0, abs=1e-9)


def test_resolve_mapping_and_embedded_match():
    assert p.resolve(p.parse("2@y2{DE}"), {1: "AAAA", 2: "PEPTIDE"}).sequence == "DE"
    assert p.parse("b2").resolve({1: "PEPTIDE"}).sequence == "PE"


@pytest.mark.parametrize(
    "text,context",
    [
        ("2@y2", {1: "PEPTIDE"}),
        ("y8", "PEPTIDE"),
        ("y2{PE}", "PEPTIDE"),
        ("m1:3", "PEPTIDE"),
        ("m2:7", "PEPTIDE"),
        ("?", "PEPTIDE"),
        ("p", ""),
        ("b2{PE/2}", "PEPTIDE"),
    ],
)
def test_resolution_errors(text, context):
    with pytest.raises(ValueError):
        p.parse(text).resolve(context)


# (annotation, the fragment the analyte's residues select, the m/z paftacular 2.0.2 published)
UNLOCALIZED_CASES = [
    ("b2^1", "PE", 307.068964),
    ("y2^1", "DE", 343.053708),
    ("m2:4^1", "EPT", 408.116643),
]


@pytest.mark.parametrize(("text", "sequence", "published_mz"), UNLOCALIZED_CASES)
def test_resolve_excludes_an_unlocalized_modification(text, sequence, published_mz):
    """An unlocalized modification is assigned to no residue, so no fragment may claim it."""
    with pytest.warns(UserWarning, match="unlocalized modification"):
        resolved = p.parse(text).resolve("[Phospho]?PEPTIDE/2")
    assert resolved.resolved_sequence == sequence
    assert resolved.mz() == pytest.approx(p.parse(text).resolve("PEPTIDE/2").mz(), rel=0, abs=1e-9)
    # the published value was the fragment plus the whole phospho mass, on every fragment
    assert published_mz - resolved.mz() == pytest.approx(79.966331, rel=0, abs=1e-6)


def test_resolve_drops_a_labile_modification_from_the_sequence():
    """A fragment loses a labile modification, so the resolved text must not show one."""
    resolved = p.parse("b2^1").resolve("{Glycan:HexNAc}PEPTIDE/2")
    assert resolved.resolved_sequence == "PE"
    assert resolved.mz() == pytest.approx(p.parse("b2^1").resolve("PEPTIDE/2").mz(), rel=0, abs=1e-9)


def test_resolve_precursor_keeps_both():
    resolved = p.resolve(p.parse("p^2"), "[Phospho]?{Glycan:HexNAc}PEPTIDE")
    assert resolved.resolved_sequence == "{Glycan:HexNAc}[Phospho]?PEPTIDE"
    assert resolved.mz() == pytest.approx(pt.parse("[Phospho]?{Glycan:HexNAc}PEPTIDE").mz(charge=2), rel=0, abs=1e-9)


def test_an_embedded_unlocalized_modification_is_kept_verbatim():
    """mzPAF 4.4.3 uses the embedded sequence as written, so its own modification counts."""
    assert p.parse("b2{[Phospho]?PE}^1").mz() == pytest.approx(307.068964, rel=0, abs=1e-6)


def test_resolving_an_embedded_sequence_keeps_its_own_modification_silently():
    """The embedded sequence is the fragment's own claim, so the analyte's slice only checks agreement."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        resolved = p.parse("b2{[Phospho]?PE}^1").resolve("[Phospho]?PEPTIDE")
        labile = p.parse("b2{{Glycan:HexNAc}PE}^1").resolve("{Glycan:HexNAc}PEPTIDE")
    assert resolved.resolved_sequence == "[Phospho]?PE"
    assert resolved.mz() == pytest.approx(307.068964, rel=0, abs=1e-6)
    assert labile.resolved_sequence == "{Glycan:HexNAc}PE"
    assert labile.mz() == pytest.approx(227.102633, rel=0, abs=1e-6)
