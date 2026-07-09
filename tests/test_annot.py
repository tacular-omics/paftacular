from collections import Counter

import pytest
from tacular import ELEMENT_LOOKUP

from paftacular import PafAnnotation, parse, parse_multi
from paftacular.comps import (
    Adduct,
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    IsotopeSpecification,
    MassError,
    NamedCompound,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)
from paftacular.constants import AminoAcids, IonSeries


class TestPafAnnotationBasics:
    """Test basic PafAnnotation functionality"""

    def test_simple_peptide_ion(self):
        """Test basic peptide ion annotation"""
        annotation = PafAnnotation(ion_type=PeptideIon(series=IonSeries.B, position=2))
        assert annotation.charge == 1
        assert annotation.ion_type.series == IonSeries.B
        assert annotation.serialize() == "b2"

    def test_peptide_ion_with_sequence(self):
        """Test peptide ion with sequence specification"""
        annotation = PafAnnotation(ion_type=PeptideIon(series=IonSeries.Y, position=3, sequence="PEP"))
        assert annotation.serialize() == "y3{PEP}"

    def test_charge_validation(self):
        """Test that charge must be >= 1"""
        with pytest.raises(ValueError, match="Charge must be >= 1"):
            PafAnnotation(
                ion_type=PeptideIon(series=IonSeries.B, position=2),
                charge=0,
            )

    def test_confidence_validation(self):
        """Test confidence must be between 0 and 1"""
        with pytest.raises(ValueError, match="Confidence must be between"):
            PafAnnotation(
                ion_type=PeptideIon(series=IonSeries.B, position=2),
                confidence=1.5,
            )


class TestIonTypes:
    """Test all ion type variations"""

    def test_peptide_ions(self):
        """Test various peptide ion series"""
        series_tests = [
            (IonSeries.A, 1, "a1"),
            (IonSeries.B, 2, "b2"),
            (IonSeries.C, 3, "c3"),
            (IonSeries.X, 4, "x4"),
            (IonSeries.Y, 5, "y5"),
            (IonSeries.Z, 6, "z6"),
        ]
        for series, pos, expected in series_tests:
            ann = PafAnnotation(ion_type=PeptideIon(series=series, position=pos))
            assert ann.serialize() == expected

    def test_internal_fragment(self):
        """Test internal fragment notation"""
        annotation = PafAnnotation(ion_type=InternalFragment(start_position=2, end_position=5))
        assert annotation.serialize() == "m2:5"

        # With sequence
        annotation = PafAnnotation(ion_type=InternalFragment(start_position=2, end_position=5, sequence="PEP"))
        assert annotation.serialize() == "m2:5{PEP}"

    def test_precursor_ion(self):
        """Test precursor ion"""
        annotation = PafAnnotation(ion_type=PrecursorIon())
        assert annotation.serialize() == "p"

    def test_immonium_ion(self):
        """Test immonium ions"""
        # Simple immonium
        annotation = PafAnnotation(ion_type=ImmoniumIon(amino_acid=AminoAcids.A))
        assert annotation.serialize() == "IA"

        # With modification
        annotation = PafAnnotation(ion_type=ImmoniumIon(amino_acid=AminoAcids.K, modification="Acetyl"))
        assert annotation.serialize() == "IK[Acetyl]"

    def test_reference_ion(self):
        """Test reference ion"""
        annotation = PafAnnotation(ion_type=ReferenceIon(name="TMT126"))
        assert annotation.serialize() == "r[TMT126]"

    def test_chemical_formula(self):
        """Test chemical formula ion"""
        annotation = PafAnnotation(ion_type=ChemicalFormula(formula="C6H12O6"))
        assert annotation.serialize() == "f{C6H12O6}"

    def test_named_compound(self):
        """Test named compound"""
        annotation = PafAnnotation(ion_type=NamedCompound(name="Glucose"))
        assert annotation.serialize() == "_{Glucose}"

    def test_smiles_compound(self):
        """Test SMILES compound"""
        annotation = PafAnnotation(ion_type=SMILESCompound(smiles="CC(C)O"))
        assert annotation.serialize() == "s{CC(C)O}"

    def test_unknown_ion(self):
        """Test unknown/unannotated ion"""
        # Without label
        annotation = PafAnnotation(ion_type=UnknownIon())
        assert annotation.serialize() == "?"

        # With label
        annotation = PafAnnotation(ion_type=UnknownIon(label=42))
        assert annotation.serialize() == "?42"


class TestModifiers:
    """Test annotation modifiers"""

    def test_neutral_loss_water(self):
        """Test water loss"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            neutral_losses=(NeutralLoss(count=-1, base_formula="H2O"),),
        )
        assert annotation.serialize() == "b2-H2O"

    def test_neutral_loss_ammonia(self):
        """Test ammonia loss"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.Y, position=3),
            neutral_losses=(NeutralLoss(count=-1, base_formula="NH3"),),
        )
        assert annotation.serialize() == "y3-NH3"

    def test_multiple_neutral_losses(self):
        """Test multiple neutral losses"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            neutral_losses=(
                NeutralLoss(count=-1, base_formula="H2O"),
                NeutralLoss(count=-1, base_formula="NH3"),
            ),
        )
        assert annotation.serialize() == "b2-H2O-NH3"

    def test_neutral_gain(self):
        """Test neutral gain (addition)"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            neutral_losses=(NeutralLoss(count=1, base_formula="H2O"),),
        )
        assert annotation.serialize() == "b2+H2O"

    def test_isotope_specification(self):
        """Test isotope notation"""
        # Single isotope
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            isotopes=(IsotopeSpecification(count=1),),
        )
        assert annotation.serialize() == "b2+i"

        # Multiple isotopes
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            isotopes=(IsotopeSpecification(count=2),),
        )
        assert annotation.serialize() == "b2+2i"

    def test_adducts_single(self):
        """Test single adduct"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            adducts=(Adduct(count=1, base_formula="Na"),),
        )
        assert annotation.serialize() == "b2[M+Na]"

    def test_adducts_multiple(self):
        """Test multiple adducts"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            adducts=(
                Adduct(count=2, base_formula="H"),
                Adduct(count=1, base_formula="Na"),
            ),
        )
        assert annotation.serialize() == "b2[M+2H+Na]"

    def test_charge_state(self):
        """Test charge state annotation"""
        # Default charge (not shown)
        annotation = PafAnnotation(ion_type=PeptideIon(series=IonSeries.B, position=2))
        assert annotation.serialize() == "b2"

        # Explicit charge
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            charge=2,
        )
        assert annotation.serialize() == "b2^2"

    def test_mass_error_da(self):
        """Test mass error in Daltons"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            mass_error=MassError(value=0.5, unit="da"),
        )
        assert annotation.serialize() == "b2/0.5"

    def test_mass_error_ppm(self):
        """Test mass error in ppm"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            mass_error=MassError(value=10.0, unit="ppm"),
        )
        assert annotation.serialize() == "b2/10ppm"

    def test_confidence(self):
        """Test confidence score"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            confidence=0.95,
        )
        assert annotation.serialize() == "b2*0.95"


class TestComplexAnnotations:
    """Test complex combinations of modifiers"""

    def test_full_annotation(self):
        """Test annotation with all modifiers"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2, sequence="PEP"),
            analyte_reference=1,
            is_auxiliary=True,
            neutral_losses=(NeutralLoss(count=-1, base_formula="H2O"),),
            isotopes=(IsotopeSpecification(count=1),),
            adducts=(Adduct(count=1, base_formula="Na"),),
            charge=2,
            mass_error=MassError(value=5.0, unit="ppm"),
            confidence=0.9,
        )
        serialized = annotation.serialize()
        assert "&1@" in serialized
        assert "b2{PEP}" in serialized
        assert "-H2O" in serialized
        assert "+i" in serialized
        assert "[M+Na]" in serialized
        assert "^2" in serialized
        assert "/5ppm" in serialized
        assert "*0.9" in serialized


class TestMassCalculations:
    """Test mass calculations"""

    def test_basic_mass(self):
        """Test basic mass calculation"""
        annotation = PafAnnotation(ion_type=PeptideIon(series=IonSeries.B, position=2))
        # Should have some positive mass
        assert annotation.mass() > 0

    def test_mass_with_neutral_loss(self):
        """Test mass with neutral loss"""
        base = PafAnnotation(ion_type=PeptideIon(series=IonSeries.B, position=2))
        with_loss = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            neutral_losses=(NeutralLoss(count=-1, base_formula="H2O"),),
        )
        # Mass should be reduced by water
        assert with_loss.mass() < base.mass()

    def test_mass_with_charge(self):
        """Test mass calculation includes protonation"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.B, position=2),
            charge=2,
        )
        # Should add 2 protons
        assert annotation.mass() > 0

    @pytest.mark.parametrize(
        "amino_acid,literature_mz",
        [
            (AminoAcids.G, 30.034),
            (AminoAcids.A, 44.050),
            (AminoAcids.P, 70.065),
            (AminoAcids.V, 72.081),
            (AminoAcids.L, 86.096),
            (AminoAcids.F, 120.081),
            (AminoAcids.Y, 136.076),
            (AminoAcids.W, 159.092),
        ],
    )
    def test_immonium_mass_matches_literature(self, amino_acid, literature_mz):
        """ImmoniumIon mass must equal published immonium ion m/z values.

        Regression test: ImmoniumIon.mass() previously (accidentally) used the
        internal by-fragment shift instead of the immonium-specific (-CO) shift.
        """
        annotation = PafAnnotation(ion_type=ImmoniumIon(amino_acid=amino_acid, modification=None), charge=1)
        assert annotation.mass() == pytest.approx(literature_mz, abs=1e-3)

    def test_internal_by_fragment_mass_equals_residue_sum(self):
        """A 'by' internal fragment's mass shift must be zero (default, per mzPAF)."""
        from tacular import AA_LOOKUP

        sequence = "PTI"
        proton_mass = 1.007276466812
        residue_sum = sum(AA_LOOKUP[AminoAcids(c)].get_mass(True) for c in sequence) + proton_mass
        annotation = PafAnnotation(
            ion_type=InternalFragment(start_position=3, end_position=5, sequence=sequence),
            charge=1,
        )
        assert annotation.mass() == pytest.approx(residue_sum, rel=1e-9)

    def test_internal_fragment_respects_backbone_cleavage_type(self):
        """A non-default (e.g. 'ax') internal fragment must not silently reuse the 'by' shift.

        Regression test: InternalFragment previously ignored its own nterm_ion_type/
        cterm_ion_type fields and always used the 'by' (0 shift) lookup.
        """
        from tacular import FRAGMENT_ION_LOOKUP

        by_frag = InternalFragment(start_position=3, end_position=5, sequence="PTI")
        ax_frag = InternalFragment(start_position=3, end_position=5, sequence="PTI", nterm_ion_type=IonSeries.A, cterm_ion_type=IonSeries.X)
        assert by_frag.mass() != ax_frag.mass()
        assert ax_frag.mass() == pytest.approx(FRAGMENT_ION_LOOKUP["ax"].get_mass(True) + by_frag.mass(), rel=1e-9)

    def test_internal_fragment_requires_both_backbone_cleavage_types(self):
        """Setting only one of nterm_ion_type/cterm_ion_type must raise, not silently default."""
        with pytest.raises(ValueError):
            InternalFragment(start_position=1, end_position=2, nterm_ion_type=IonSeries.A)

    def test_immonium_composition_matches_mass_with_modification(self):
        """A modified ImmoniumIon's composition-implied mass must match .mass().

        Regression test: ImmoniumIon.composition accumulated the modification's
        composition into an empty Counter before the amino acid's own composition;
        Counter's `+=` drops any element whose running total is <= 0 at that step,
        silently losing an atom-removing modification (e.g. Deamidated).
        """
        pytest.importorskip("peptacular")

        ion = ImmoniumIon(amino_acid=AminoAcids.A, modification="Deamidated")
        mass_from_composition = sum(elem.mass for elem, count in ion.composition.items() for _ in range(count))
        assert mass_from_composition == pytest.approx(ion.mass(), abs=1e-6)


class TestComposition:
    """Test elemental composition calculations"""

    def test_basic_composition(self):
        """Test basic composition calculation"""
        annotation = PafAnnotation(ion_type=ChemicalFormula(formula="C6H12O6"))
        comp = annotation.comp()
        assert isinstance(comp, Counter)
        # Should contain C, H, O elements
        assert any(elem.symbol == "C" for elem in comp)
        assert any(elem.symbol == "H" for elem in comp)
        assert any(elem.symbol == "O" for elem in comp)

    def test_composition_with_charge(self):
        """Test composition includes protonation"""
        annotation = PafAnnotation(
            ion_type=ChemicalFormula(formula="C6H12O6"),
            charge=2,
        )
        comp = annotation.comp()
        # Should have protons added
        h_element = ELEMENT_LOOKUP["H"]
        assert h_element in comp

    def test_composition_with_neutral_loss(self):
        """Test composition with neutral loss"""
        annotation = PafAnnotation(
            ion_type=ChemicalFormula(formula="C6H12O6"),
            neutral_losses=(NeutralLoss(count=-1, base_formula="H2O"),),
        )
        comp = annotation.comp()
        # Composition should be affected by loss
        assert isinstance(comp, Counter)


class TestParsing:
    """Test parsing mzPAF strings"""

    def test_parse_simple(self):
        """Test parsing simple annotation"""
        annotation = parse("b2")
        assert isinstance(annotation.ion_type, PeptideIon)
        assert annotation.ion_type.series == IonSeries.B
        assert annotation.ion_type.position == 2

    def test_parse_with_sequence(self):
        """Test parsing with sequence"""
        annotation = parse("y3{PEP}")
        assert isinstance(annotation.ion_type, PeptideIon)
        assert annotation.ion_type.series == IonSeries.Y
        assert annotation.ion_type.position == 3
        assert annotation.ion_type.sequence == "PEP"

    def test_parse_with_neutral_loss(self):
        """Test parsing with neutral loss"""
        annotation = parse("b2-H2O")
        assert len(annotation.neutral_losses) == 1
        assert annotation.neutral_losses[0].count == -1

    def test_parse_with_adduct(self):
        """Test parsing with adduct"""
        annotation = parse("b2[M+Na]")
        assert len(annotation.adducts) == 1
        assert annotation.adducts[0].count == 1

    def test_parse_with_charge(self):
        """Test parsing with charge"""
        annotation = parse("b2^2")
        assert annotation.charge == 2

    def test_parse_multiple_annotations(self):
        """Test parsing comma-separated annotations"""
        annotations = parse_multi("b2, y3, a1")
        assert len(annotations) == 3
        assert isinstance(annotations[0].ion_type, PeptideIon)
        assert isinstance(annotations[1].ion_type, PeptideIon)
        assert isinstance(annotations[2].ion_type, PeptideIon)
        assert annotations[0].ion_type.series == IonSeries.B
        assert annotations[1].ion_type.series == IonSeries.Y
        assert annotations[2].ion_type.series == IonSeries.A

    def test_parse_roundtrip(self):
        """Test that parse -> serialize -> parse gives same result"""
        test_strings = [
            "b2",
            "y3{PEP}",
            "b2-H2O",
            "b2[M+Na]^2",
            "m2:5",
            "IA",
            "p",
            "r[TMT126]",
            "f{C6H12O6}",
        ]
        for test_str in test_strings:
            annotation = parse(test_str)
            serialized = annotation.serialize()
            re_parsed = parse(serialized)
            assert annotation.ion_type == re_parsed.ion_type

    def test_invalid_annotation(self):
        """Test that invalid annotations raise errors"""
        with pytest.raises(ValueError, match="Invalid mzPAF annotation"):
            parse("invalid_annotation!")


class TestFormulaProperties:
    """Test formula string properties"""

    def test_formula_property(self):
        """Test formula property returns string"""
        annotation = PafAnnotation(ion_type=ChemicalFormula(formula="C6H12O6"))
        formula = annotation.formula()
        assert isinstance(formula, str)
        assert len(formula) > 0

    def test_proforma_formula_property(self):
        """Test ProForma formula property"""
        annotation = PafAnnotation(ion_type=ChemicalFormula(formula="C6H12O6"))
        formula = annotation.proforma_formula()
        assert isinstance(formula, str)
        assert len(formula) > 0


class TestSequenceProperty:
    """Test sequence property"""

    def test_sequence_for_peptide_ion(self):
        """Test sequence property for peptide ions"""
        annotation = PafAnnotation(ion_type=PeptideIon(series=IonSeries.B, position=2, sequence="PEP"))
        assert annotation.sequence == "PEP"

    def test_sequence_for_internal_fragment(self):
        """Test sequence property for internal fragments"""
        annotation = PafAnnotation(ion_type=InternalFragment(start_position=2, end_position=5, sequence="PEP"))
        assert annotation.sequence == "PEP"

    def test_sequence_none_for_other_ions(self):
        """Test sequence property is None for non-peptide ions"""
        annotation = PafAnnotation(ion_type=PrecursorIon())
        assert annotation.sequence is None


class TestAsDictMethod:
    """Test as_dict method"""

    def test_as_dict_basic(self):
        """Test as_dict returns proper dictionary"""
        annotation = PafAnnotation(ion_type=PeptideIon(series=IonSeries.B, position=2))
        result = annotation.as_dict()
        assert isinstance(result, dict)
        assert "ion" in result
        assert "charge" in result
        assert result["charge"] == 1

    def test_as_dict_with_all_fields(self):
        """Test as_dict with all fields populated"""
        annotation = PafAnnotation(
            ion_type=PeptideIon(series=IonSeries.Y, position=5),
            analyte_reference=1,
            is_auxiliary=True,
            neutral_losses=(NeutralLoss(count=-1, base_formula="H2O"),),
            isotopes=(IsotopeSpecification(count=1),),
            adducts=(Adduct(count=1, base_formula="Na"),),
            charge=2,
            mass_error=MassError(0.5, "ppm"),
            confidence=0.95,
        )
        result = annotation.as_dict()
        assert result["analyte_reference"] == 1
        assert result["is_auxiliary"] is True
        assert len(result["neutral_losses"]) == 1
        assert len(result["isotopes"]) == 1
        assert len(result["adducts"]) == 1
        assert result["charge"] == 2
        assert result["mass_error"] is not None
        assert result["confidence"] == 0.95


class TestReprMethod:
    """Test __repr__ method"""

    def test_repr_basic(self):
        """Test __repr__ returns string"""
        annotation = PafAnnotation(ion_type=PeptideIon(series=IonSeries.B, position=2))
        repr_str = repr(annotation)
        assert isinstance(repr_str, str)
        assert "PafAnnotation" in repr_str


class TestParseEmptyAndMultiple:
    """Test parse with empty and multiple annotations"""

    def test_parse_empty_string(self):
        """Test parsing empty string returns empty list"""
        result = parse_multi("")
        assert result == []

    def test_parse_multiple_annotations(self):
        """Test parsing multiple comma-separated annotations"""
        result = parse_multi("b5, y10")
        assert len(result) == 2
        assert isinstance(result[0].ion_type, PeptideIon)
        assert isinstance(result[1].ion_type, PeptideIon)

    def test_parse_multiple_with_whitespace(self):
        """Test parsing handles whitespace correctly"""
        result = parse_multi("b5,  y10  , c3")
        assert len(result) == 3


class TestCompositionMethods:
    """Test composition-related methods"""

    def test_composition_basic(self):
        """Test composition method returns Counter"""
        annotation = PafAnnotation(ion_type=PrecursorIon())
        comp = annotation.comp()
        assert isinstance(comp, Counter)

    def test_dict_composition(self):
        """Test dict_composition returns dict with string keys"""
        annotation = PafAnnotation(ion_type=PrecursorIon())
        comp = annotation.dict_composition()
        assert isinstance(comp, dict)
        for key in comp.keys():
            assert isinstance(key, str)

    def test_composition_with_neutral_loss(self):
        """Test composition includes neutral loss"""
        annotation = PafAnnotation(
            ion_type=PrecursorIon(),
            neutral_losses=(NeutralLoss(count=-1, base_formula="H2O"),),
        )
        comp = annotation.comp()
        # Should have adjusted composition
        assert isinstance(comp, Counter)

    def test_composition_with_adduct(self):
        """Test composition includes adduct"""
        annotation = PafAnnotation(
            ion_type=PrecursorIon(),
            adducts=(Adduct(count=1, base_formula="Na"),),
        )
        comp = annotation.comp()
        # Should include Na
        assert any("Na" in str(elem) for elem in comp.keys())

    def test_composition_with_charge_no_adducts(self):
        """Test composition includes protons for charge when no adducts"""
        annotation = PafAnnotation(ion_type=PrecursorIon(), charge=2)
        comp = annotation.comp()
        # Should have H added for charge
        h_elem = ELEMENT_LOOKUP["H"]
        assert h_elem in comp


class TestMzMethod:
    """Test m/z calculation"""

    def test_mz_basic(self):
        """Test m/z calculation"""
        annotation = PafAnnotation(ion_type=PrecursorIon(), charge=2)
        mz = annotation.mz()
        assert isinstance(mz, float)
        assert mz > 0

    def test_mz_equals_mass_when_charge_one(self):
        """Test m/z equals mass when charge is 1"""
        annotation = PafAnnotation(ion_type=PrecursorIon(), charge=1)
        mass = annotation.mass()
        mz = annotation.mz()
        assert abs(mass - mz) < 0.01


class TestInvalidAnnotations:
    """Test error handling for invalid annotations"""

    def test_invalid_annotation_string(self):
        """Test parsing invalid annotation raises ValueError"""
        with pytest.raises(ValueError, match="Invalid mzPAF annotation"):
            parse("not_valid_annotation")

    def test_negative_confidence(self):
        """Test negative confidence raises ValueError"""
        with pytest.raises(ValueError, match="Confidence must be between"):
            PafAnnotation(
                ion_type=PrecursorIon(),
                confidence=-0.1,
            )

    def test_confidence_greater_than_one(self):
        """Test confidence > 1 raises ValueError"""
        with pytest.raises(ValueError, match="Confidence must be between"):
            PafAnnotation(
                ion_type=PrecursorIon(),
                confidence=1.1,
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
