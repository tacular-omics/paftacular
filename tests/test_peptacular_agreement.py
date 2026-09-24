"""paftacular and peptacular agree on every ion m/z for sequences with named modifications.

Both add a named modification by its listed database mass (Unimod or PSI-MOD, 6 decimals),
not by its composition, so the two packages match to 1e-9 Da.
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
