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
        assert paf_annot.get_mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_b_pos(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.B, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.PeptideIon)
        assert "b3{PEP}^2" == str(paf_annot)
        assert paf_annot.get_mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_y(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.Y, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.PeptideIon)
        assert "y3{IDE}^2" == str(paf_annot)
        assert paf_annot.get_mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_i(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.IMMONIUM, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.ImmoniumIon)
        assert "IP^2" == str(paf_annot)
        assert paf_annot.get_mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_i_mod(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEP[+10]TIDE/2")
        frag = annot.frag(ion_type=pt.IonType.IMMONIUM, charge=2, position=3)
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.ImmoniumIon)
        assert "IP[+10]^2" == str(paf_annot)
        assert paf_annot.get_mass() == pytest.approx(frag.mass, rel=1e-6)

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
        assert paf_annot.get_mass() == pytest.approx(frag.mass, rel=1e-6)

    def test_conversion_internal_ay(self):
        """Test converting a peptacular Fragment to mzPAF format string"""
        annot = pt.parse("PEPTIDE/2")
        frag = annot.frag(ion_type=pt.IonType.AY, charge=2, position=(3, 5))
        paf_annot: paf.PafAnnotation = paf.to_mzpaf(frag)
        assert isinstance(paf_annot.ion_type, paf.InternalFragment)
        assert "m3:5{PTI}-CO^2" == str(paf_annot)
        assert paf_annot.get_mass() == pytest.approx(frag.mass, rel=1e-6)

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
        assert annotation.get_mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)
        assert paf.parse(annotation.serialize()) == annotation

    def test_precursor_conversion_retains_context(self):
        fragment = pt.parse("PEPTIDE").frag(ion_type="p", charge=2)
        annotation = paf.to_mzpaf(fragment)
        assert annotation.serialize() == "p^2"
        assert annotation.sequence == "PEPTIDE"
        assert annotation.get_mass() == pytest.approx(fragment.mass, rel=0, abs=1e-6)
        assert paf.PafAnnotation.from_dict(annotation.to_dict()) == annotation
        assert paf.to_mzpaf(fragment, include_sequence=False).sequence is None


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
    assert annotation.get_mass() == pytest.approx(frag.mass, rel=0, abs=1e-6)


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
    parsed = paf.parse(annotation.serialize())
    assert parsed == annotation
    assert parsed.mz() == pytest.approx(expected_mz, rel=0, abs=1e-5)


@pytest.mark.parametrize(("ion_type", "position"), [("d", 5), ("v", 2), ("da-threonine", 3), ("wb-isoleucine", 8)])
def test_conversion_side_chain_ions_negative_charge(ion_type, position):
    # Negative charge is kept and written as ^-1, so the round trip is lossless.
    fragment = _side_chain_fragment(ion_type, position, -1)
    annotation = paf.to_mzpaf(fragment)
    assert annotation.charge == -1
    assert annotation.serialize().endswith("^-1")
    assert paf.parse(annotation.serialize()) == annotation
    positive = paf.to_mzpaf(_side_chain_fragment(ion_type, position, 1))
    assert positive.get_mass() - annotation.get_mass() == pytest.approx(2 * PROTON, rel=0, abs=1e-6)


_PARITY_CASES = [
    (
        "[Acetyl]-PEM[Oxidation]TIDEK/2",
        {
            "ion_types": ("b", "y", "a", "c", "x", "z", "z+H", "c-H", "p", "i", "by", "ax", "cz"),
            "charges": [1, 2],
            "neutral_deltas": ("H2O", "NH3", "HCOOH", "HCONH2"),
            "max_ndeltas": 2,
        },
    ),
    ("PEPTIDEK", {"ion_types": ("b", "y"), "charges": [1], "isotopes": (0, 1, {"13C": 2}), "deltas": (None, 15.9949, "H2O", {"Na": 1})}),
    ("PEPTIDEK", {"ion_types": ("b", "y", "d", "w", "v", "p"), "charges": [-1, -2, 3]}),
    ("PEPTIDEK/[Na:z+1,H:z+1]", {"ion_types": ("b", "y", "p")}),
]


@pytest.mark.parametrize(("sequence", "options"), _PARITY_CASES)
def test_to_mzpaf_parity_with_peptacular(sequence, options):
    """to_mzpaf(f) serializes to text that parses back to the same annotation and m/z as f."""
    for fragment in pt.fragment(sequence, **options):
        annotation = paf.to_mzpaf(fragment)
        text = annotation.serialize()
        parsed = paf.parse(text)
        if isinstance(annotation.ion_type, paf.PrecursorIon):
            # p carries its sequence as resolved context, not in mzPAF text.
            assert parsed.serialize() == text
            assert annotation.mz() == pytest.approx(fragment.mz, rel=0, abs=1e-5), text
            continue
        assert parsed == annotation, text
        assert parsed.mz() == pytest.approx(fragment.mz, rel=0, abs=1e-5), text
        if fragment.charge_state > 0 and not isinstance(annotation.ion_type, paf.ImmoniumIon):
            # peptacular's own label names the same ion. Its immonium labels drop terminal
            # modifications (IP for [Acetyl]-P), so those are checked against fragment.mz only.
            assert paf.parse(fragment.to_mzpaf()).mz() == pytest.approx(parsed.mz(), rel=0, abs=1e-5), text


def test_terminal_modified_immonium_keeps_its_mass():
    fragment = pt.fragment("[Acetyl]-PEK", ion_types=("i",), charges=[1])[0]
    annotation = paf.to_mzpaf(fragment)
    assert annotation.serialize() == "IP[Acetyl]"
    assert annotation.mz() == pytest.approx(fragment.mz, rel=0, abs=1e-5)


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        ("<[Oxidation]@P>PEK", "IP[Oxidation]"),
        ("<13C>PEK", "IP+4i13C"),
        ("<13C><15N>PEK", "IP+4i13C+i15N"),
        ("<13C>P[Acetyl]EK", "IP[Acetyl]+6i13C"),
        ("<[Oxidation]@P><13C>PEK", "IP[Oxidation]+4i13C"),
    ],
)
def test_immonium_keeps_global_mods_and_labels(sequence, expected):
    fragment = pt.fragment(sequence, ion_types=("i",), charges=[1])[0]
    annotation = paf.to_mzpaf(fragment)
    assert annotation.serialize() == expected
    assert annotation.serialize() == fragment.to_mzpaf()
    assert paf.parse(expected).mz() == pytest.approx(fragment.mz, rel=0, abs=1e-6)


def test_immonium_mass_uses_exact_constants():
    fragment = pt.fragment("PEK", ion_types=("i",), charges=[1])[0]
    # Only the proton may differ: peptacular before 5.0 charged with H minus an electron.
    assert paf.parse("IP").mz() == pytest.approx(fragment.mz, rel=0, abs=5e-8)


@pytest.mark.parametrize(
    ("ion_type", "loss"),
    [("cy", "+NH3"), ("bz", "-NH3"), ("az", "-HCONH2"), ("ax", "-H2"), ("bx", "+CO-H2"), ("cx", "+CHNO"), ("ay", "-CO")],
)
def test_internal_offsets_use_canonical_names(ion_type, loss):
    fragment = pt.fragment("PEPTIDEK", ion_types=(ion_type,), charges=[1])[0]
    text = paf.to_mzpaf(fragment).serialize()
    assert text == fragment.to_mzpaf()
    assert text.endswith(loss)


@pytest.mark.parametrize(
    ("sequence", "ion_type", "position", "deltas"),
    [
        ("<13C>RYK", "a", 2, None),
        ("<13C>QDK", "x", 2, None),
        ("<15N>PLK", "c", 2, None),
        ("<15N>KSL", "z", 2, None),
        ("<13C>KMIH", "v", 3, None),
        ("<13C>PEPTIDE", "az", (2, 4), None),
        ("<13C>PEPTIDE", "p", None, None),
        ("<15N>KEK", "b", 2, "N-1H-3"),
        ("<13C>PEK", "b", 2, "C-1O-2"),
        ("<15N>KEK", "b", 2, {-17.026549: 1}),
    ],
)
def test_global_label_covers_offset_and_formula_deltas(sequence, ion_type, position, deltas):
    # The label replaces its element in the final ion: residues, ion offset and formula
    # deltas. Mass-only deltas and the charge carrier stay unlabelled.
    fragment = pt.parse(sequence).frag(ion_type=ion_type, position=position, charge=1, deltas=deltas)
    annotation = paf.to_mzpaf(fragment)
    if ion_type != "p":  # a precursor keeps its sequence outside the text
        assert paf.parse(annotation.serialize()).get_mass() == pytest.approx(annotation.get_mass(), rel=0, abs=1e-9)
    assert annotation.mz() == pytest.approx(fragment.mz, rel=0, abs=1e-7)


def test_immonium_label_counts_atoms_after_deltas():
    # <15N>K immonium has two N. Losing NH3 removes one of them, so one 15N label remains.
    fragment = pt.parse("<15N>K").frag(ion_type="i", position=1, charge=1, deltas="N-1H-3")
    annotation = paf.to_mzpaf(fragment)
    assert annotation.serialize() == "IK-NH3+i15N"
    assert annotation.mz() == pytest.approx(fragment.mz, rel=0, abs=1e-7)
