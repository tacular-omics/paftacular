import pytest

pt = pytest.importorskip("peptacular")

import paftacular as paf  # noqa: E402


class TestConversion:
    """Test converting peptacular objects to mzPAF format"""

    def test_conversion_b(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.B, charge=2)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.PeptideIon)
        assert "b7{PEPTIDE}^2" == str(paf_annot)
        assert paf_annot.mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_b_pos(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.B, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.PeptideIon)
        assert "b3{PEP}^2" == str(paf_annot)
        assert paf_annot.mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_y(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.Y, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.PeptideIon)
        assert "y3{IDE}^2" == str(paf_annot)
        assert paf_annot.mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_i(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.IMMONIUM, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.ImmoniumIon)
        assert "IP^2" == str(paf_annot)
        assert paf_annot.mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_i_mod(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEP[+10]TIDE/2")
        frag = annot.frag(ion_type=pt.IonType.IMMONIUM, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.ImmoniumIon)
        assert "IP[+10]^2" == str(paf_annot)
        assert paf_annot.mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_i_modx2(self):
        annot = pt.parse("PEP[+10][Oxidation]TIDE/2")
        frag = annot.frag(ion_type=pt.IonType.IMMONIUM, charge=2, position=3)
        with pytest.raises(ValueError):
            paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
            assert isinstance(paf_annot.ion_type, paf.ImmoniumIon)

    def test_conversion_internal_by(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.BY, charge=2, position=(3, 5))
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.InternalFragment)
        assert "m3:5{PTI}^2" == str(paf_annot)
        assert paf_annot.mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_internal_ay(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.AY, charge=2, position=(3, 5))
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.InternalFragment)
        assert "m3:5{PTI}-CO^2" == str(paf_annot)
        assert paf_annot.mass() == pytest.approx(frag.mass, rel=1e-6)

    @pytest.mark.parametrize(
        "options",
        [
            {"charge": "Na:z+1"},
            {"charge": "Na:z+1^2"},
            {"charge": 2, "isotopes": 1},
            {"charge": 2, "isotopes": {"13C": 2}},
            {"charge": 2, "deltas": -18.010565},
            {"charge": 2, "deltas": "H-2O-1"},
        ],
    )
    def test_conversion_preserves_modifiers(self, options):
        fragment = pt.parse("PEPTIDE").frag(ion_type="y", position=3, **options)
        annotation = paf.to_mzpaf(fragment, confidence=0.123456789, mass_error=1e-7)
        assert annotation.mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)
        assert paf.parse_single(annotation.serialize()) == annotation

    def test_precursor_conversion_retains_context(self):
        fragment = pt.parse("PEPTIDE").frag(ion_type="p", charge=2)
        annotation = paf.to_mzpaf(fragment)
        assert annotation.serialize() == "p^2"
        assert annotation.sequence == "PEPTIDE"
        assert annotation.mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)
        assert paf.PafAnnotation.from_dict(annotation.to_dict()) == annotation
        assert paf.to_mzpaf(fragment, include_annotation=False).sequence is None


if __name__ == "__main__":
    import pytest

    pytest.main([__file__])


@pytest.mark.parametrize(
    ("ion_type", "expected"),
    [
        (pt.IonType.Z, "z3{IDE}-H^2"),
        (pt.IonType.Z_RADICAL, "z3{IDE}^2"),
        (pt.IonType.Z_PLUS_H, "z3{IDE}+H^2"),
        (pt.IonType.C_MINUS_H, "c3{PEP}-H^2"),
    ],
)
def test_conversion_z_and_c_variants(ion_type, expected):
    # mzPAF 1.0.1 section 4.4.3: z is the z-dot radical (sum + H2O - NH2).
    frag = pt.parse("PEPTIDE/2").frag(ion_type=ion_type, charge=2, position=3)
    annotation = paf.to_mzpaf(frag)
    assert str(annotation) == expected
    assert annotation.mass() == pytest.approx(frag.mass, rel=0, abs=1e-6)


PROTON = 1.007276466812
SIDE_CHAIN_PEPTIDE = "PVTIKDEWTLR"


def _side_chain_fragment(ion_type: str, position: int, charge: int):
    # Built directly so the test does not depend on how the installed peptacular
    # validates or weighs these series (4.0.0 and earlier used other formulas).
    return pt.Fragment(
        ion_type=pt.IonType(ion_type),
        position=position,
        mass=0.0,
        monoisotopic=True,
        charge_state=charge,
        parent_sequence=f"{SIDE_CHAIN_PEPTIDE}/{charge}",
        parent_sequence_length=len(SIDE_CHAIN_PEPTIDE),
    )


# peptacular flags d, v, da, db, wa and wb as FORWARD|AA_SPECIFIC_FWD or
# BACKWARD|AA_SPECIFIC_BWD, so to_mzpaf must test flag membership, not equality.
# Expected m/z values follow mzPAF 1.0.1 section 4.4.3 (they match peptacular main).
@pytest.mark.parametrize(
    ("ion_type", "position", "expected", "mz"),
    [
        ("d-valine", 2, "d2{PV}", 154.110064),
        ("d", 5, "d5{PVTIK}", 453.294571),
        ("v", 2, "v2{LR}", 230.124766),
        ("w", 2, "w2{LR}", 229.129517),
        ("w-valine", 10, "w10{VTIKDEWTLR}", 1229.652464),
        ("da-threonine", 3, "da3{PVT}", 255.157743),
        ("db-threonine", 3, "db3{PVT}", 253.178478),
        ("da-isoleucine", 4, "da4{PVTI}", 368.241807),
        ("db-isoleucine", 4, "db4{PVTI}", 354.226157),
        ("wa-threonine", 3, "wa3{TLR}", 358.208495),
        ("wb-threonine", 3, "wb3{TLR}", 356.229231),
        ("wa-isoleucine", 8, "wa8{IKDEWTLR}", 1029.536372),
        ("wb-isoleucine", 8, "wb8{IKDEWTLR}", 1015.520721),
    ],
)
@pytest.mark.parametrize("charge", [1, 2])
def test_conversion_side_chain_ions(ion_type, position, expected, mz, charge):
    annotation = paf.to_mzpaf(_side_chain_fragment(ion_type, position, charge))
    assert str(annotation) == (expected if charge == 1 else f"{expected}^{charge}")
    expected_mz = (mz + (charge - 1) * PROTON) / charge
    assert annotation.mz() == pytest.approx(expected_mz, rel=0, abs=1e-5)
    parsed = paf.parse_single(annotation.serialize())
    assert parsed == annotation
    assert parsed.mz() == pytest.approx(expected_mz, rel=0, abs=1e-5)


@pytest.mark.parametrize(("ion_type", "position"), [("d", 5), ("v", 2), ("da-threonine", 3), ("wb-isoleucine", 8)])
def test_conversion_side_chain_ions_negative_charge(ion_type, position):
    # mzPAF 1.0.1 charges are positive integers, as for every other series.
    with pytest.raises(ValueError, match="Charge"):
        paf.to_mzpaf(_side_chain_fragment(ion_type, position, -1))
