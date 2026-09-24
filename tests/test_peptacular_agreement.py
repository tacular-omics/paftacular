"""paftacular and peptacular agree on ion m/z for sequences with named modifications.

Both add a named modification by its listed database mass (Unimod or PSI-MOD, 6 decimals),
not by its composition, so plain fragment and precursor ions match to 1e-9 Da. Labile
modifications count for the precursor only: fragments lose them, in both packages.

Ions with neutral losses or isotope peaks match too: both packages add the loss or isotope
shift to the plain ion, which keeps the listed modification masses.
"""

import pytest

pt = pytest.importorskip("peptacular")

import paftacular as pft  # noqa: E402

PEPTIDES = [
    "PEM[Oxidation]TIDEK",
    "[Acetyl]-S[Phospho]EQUENC[Carbamidomethyl]EK",
    "AC[Carbamidomethyl]DEFGHIK[Label:13C(6)15N(2)]",
    "PEPS[MOD:00046]TIDEK",
    "M[MOD:00719]AGIC[UNIMOD:4]PEPTIDER",
    "N[Deamidated]Q[Gln->pyro-Glu]TY[Nitro]K[GG]R",
    "<[Carbamidomethyl]@C>PEC[Oxidation]TIM[Dioxidation]DER",
    "{Glycan:Hex}PEPTIDEK",
    "{Glycan:Hex}PEM[Oxidation]TIDEK",
]

# mzPAF z is the z-dot radical, which peptacular calls "z.".
SERIES = {"a": "a", "b": "b", "c": "c", "x": "x", "y": "y", "z": "z."}


def _cases():
    for peptide in PEPTIDES:
        length = len(pt.parse(peptide))
        for series in SERIES:
            for position in (1, 3, length - 1):
                for charge in (1, 2, 3):
                    yield peptide, series, position, charge


@pytest.mark.parametrize(("peptide", "series", "position", "charge"), list(_cases()))
def test_mz_matches_peptacular(peptide, series, position, charge):
    expected = pt.parse(peptide).frag(ion_type=SERIES[series], charge=charge, position=position)
    resolved = pft.parse(f"{series}{position}^{charge}").resolve(peptide)
    assert resolved.mz() == pytest.approx(expected.mz, rel=0, abs=1e-9)
    assert resolved.get_mass() == pytest.approx(expected.mass, rel=0, abs=1e-9)
    assert pft.to_mzpaf(expected).mz() == pytest.approx(expected.mz, rel=0, abs=1e-9)


@pytest.mark.parametrize("peptide", PEPTIDES)
@pytest.mark.parametrize("charge", [1, 2])
def test_precursor_mz_matches_peptacular(peptide, charge):
    expected = pt.parse(peptide).mz(charge=charge)
    assert pft.parse(f"p^{charge}").resolve(peptide).mz() == pytest.approx(expected, rel=0, abs=1e-9)


@pytest.mark.parametrize("peptide", PEPTIDES)
def test_average_mass_matches_peptacular(peptide):
    expected = pt.parse(peptide).frag(ion_type="y", charge=2, position=4, monoisotopic=False)
    resolved = pft.parse("y4^2").resolve(peptide)
    assert resolved.get_mass(monoisotopic=False) == pytest.approx(expected.mass, rel=0, abs=1e-9)


# (peptide, mzPAF annotation, the same ion without the delta, peptacular frag() arguments)
DELTA_CASES = [
    ("PEM[Oxidation]TIDEK", "y6-H2O", "y6", {"ion_type": "y", "position": 6, "charge": 1, "deltas": {"H2O": 1}}),
    ("PEM[Oxidation]TIDEK", "b3-NH3^2", "b3^2", {"ion_type": "b", "position": 3, "charge": 2, "deltas": {"NH3": 1}}),
    ("PEPS[Phospho]TIDEK", "y6-H3PO4", "y6", {"ion_type": "y", "position": 6, "charge": 1, "deltas": {"H3PO4": 1}}),
    ("PEM[Oxidation]TIDEK", "y6+i", "y6", {"ion_type": "y", "position": 6, "charge": 1, "isotopes": 1}),
    ("PEPS[Phospho]TIDEK", "b4+2i^2", "b4^2", {"ion_type": "b", "position": 4, "charge": 2, "isotopes": 2}),
    ("PEM[Oxidation]TIDEK", "p-H2O^2", "p^2", {"ion_type": "p", "charge": 2, "deltas": {"H2O": 1}}),
]


@pytest.mark.parametrize(("peptide", "annotation", "plain", "kwargs"), DELTA_CASES)
def test_delta_and_isotope_mz_matches_peptacular(peptide, annotation, plain, kwargs):
    expected = pt.parse(peptide).frag(**kwargs)
    assert pft.parse(annotation).resolve(peptide).mz() == pytest.approx(expected.mz, rel=0, abs=1e-9)


@pytest.mark.parametrize(("peptide", "annotation", "plain", "kwargs"), DELTA_CASES)
def test_delta_and_isotope_ion_is_plain_ion_plus_delta(peptide, annotation, plain, kwargs):
    # paftacular is self-consistent: the delta ion is the plain ion (which matches peptacular)
    # plus the exact delta.
    plain_kwargs = {key: value for key, value in kwargs.items() if key not in ("deltas", "isotopes")}
    plain_ion = pft.parse(plain).resolve(peptide)
    assert plain_ion.get_mass() == pytest.approx(pt.parse(peptide).frag(**plain_kwargs).mass, rel=0, abs=1e-9)
    delta = pft.parse(annotation).get_mass(calculate_sequence=False) - pft.parse(plain).get_mass(calculate_sequence=False)
    assert pft.parse(annotation).resolve(peptide).get_mass() == pytest.approx(plain_ion.get_mass() + delta, rel=0, abs=1e-9)
