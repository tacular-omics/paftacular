"""Tests for PafAnnotation static make_* factory methods"""

from paftacular import PafAnnotation
from paftacular.comps import (
    Adduct,
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    IsotopeSpecification,
    NamedCompound,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)
from paftacular.constants import AminoAcids, IonSeries


class TestMakePrecursor:
    """Test make_precursor factory method"""

    def test_basic_precursor(self):
        """Test creating a basic precursor annotation"""
        annot = PafAnnotation.make_precursor()
        assert isinstance(annot.ion_type, PrecursorIon)
        assert annot.charge == 1
        assert annot.serialize() == "p"

    def test_precursor_with_charge(self):
        """Test precursor with charge state"""
        annot = PafAnnotation.make_precursor(charge=2)
        assert annot.charge == 2
        assert annot.serialize() == "p^2"

    def test_precursor_with_neutral_loss(self):
        """Test precursor with neutral loss"""
        annot = PafAnnotation.make_precursor(neutral_losses=["-H2O"])
        assert len(annot.neutral_losses) == 1
        assert annot.neutral_losses[0].base_formula == "H2O"
        assert "H2O" in annot.serialize()

    def test_precursor_with_adduct(self):
        """Test precursor with adduct"""
        annot = PafAnnotation.make_precursor(adducts=["+Na"])
        assert len(annot.adducts) == 1
        assert annot.adducts[0].base_formula == "Na"
        assert "[M+Na]" in annot.serialize()

    def test_precursor_with_isotope(self):
        """Test precursor with isotope specification"""
        annot = PafAnnotation.make_precursor(isotopes=["+i"])
        assert len(annot.isotopes) == 1
        assert annot.isotopes[0].count == 1
        assert "+i" in annot.serialize()

    def test_precursor_with_mass_error(self):
        """Test precursor with mass error"""
        annot = PafAnnotation.make_precursor(mass_error=0.5, mass_error_unit="ppm")
        assert annot.mass_error is not None
        assert annot.mass_error.value == 0.5
        assert annot.mass_error.unit == "ppm"
        assert "/0.5ppm" in annot.serialize()

    def test_precursor_with_confidence(self):
        """Test precursor with confidence score"""
        annot = PafAnnotation.make_precursor(confidence=0.95)
        assert annot.confidence == 0.95
        assert "*0.95" in annot.serialize()


class TestMakePeptide:
    """Test make_peptide factory method"""

    def test_basic_peptide_b_ion(self):
        """Test creating a basic b-ion"""
        annot = PafAnnotation.make_peptide("b", 5)
        assert isinstance(annot.ion_type, PeptideIon)
        assert annot.ion_type.series == IonSeries.B
        assert annot.ion_type.position == 5
        assert annot.serialize() == "b5"

    def test_peptide_y_ion(self):
        """Test creating a y-ion"""
        annot = PafAnnotation.make_peptide("y", 10)
        assert isinstance(annot.ion_type, PeptideIon)
        assert annot.ion_type.series == IonSeries.Y
        assert annot.ion_type.position == 10
        assert annot.serialize() == "y10"

    def test_peptide_with_sequence(self):
        """Test peptide ion with sequence"""
        annot = PafAnnotation.make_peptide("b", 3, sequence="PEP")
        assert isinstance(annot.ion_type, PeptideIon)
        assert annot.ion_type.sequence == "PEP"
        assert annot.serialize() == "b3{PEP}"

    def test_peptide_with_charge(self):
        """Test peptide ion with charge state"""
        annot = PafAnnotation.make_peptide("y", 5, charge=2)
        assert annot.charge == 2
        assert annot.serialize() == "y5^2"

    def test_peptide_with_neutral_loss(self):
        """Test peptide ion with neutral loss"""
        annot = PafAnnotation.make_peptide("b", 7, neutral_losses=["-H2O", "-NH3"])
        assert len(annot.neutral_losses) == 2
        assert "-H2O" in annot.serialize()
        assert "-NH3" in annot.serialize()

    def test_peptide_with_ion_series_enum(self):
        """Test using IonSeries enum directly"""
        annot = PafAnnotation.make_peptide(IonSeries.C, 4)
        assert isinstance(annot.ion_type, PeptideIon)
        assert annot.ion_type.series == IonSeries.C
        assert annot.serialize() == "c4"


class TestMakeInternal:
    """Test make_internal factory method"""

    def test_basic_internal_fragment(self):
        """Test creating a basic internal fragment"""
        annot = PafAnnotation.make_internal(5, 10)
        assert isinstance(annot.ion_type, InternalFragment)
        assert annot.ion_type.start_position == 5
        assert annot.ion_type.end_position == 10
        assert annot.serialize().startswith("m5:10")

    def test_internal_with_sequence(self):
        """Test internal fragment with sequence"""
        annot = PafAnnotation.make_internal(2, 5, sequence="TIDE")
        assert isinstance(annot.ion_type, InternalFragment)
        assert annot.ion_type.sequence == "TIDE"
        assert "m2:5{TIDE}" in annot.serialize()

    def test_internal_with_ion_type(self):
        """Test internal fragment with specific ion type"""
        annot = PafAnnotation.make_internal(3, 7, ion_type="ay")
        # Should add the appropriate neutral loss for ay type
        assert len(annot.neutral_losses) > 0

    def test_internal_with_charge(self):
        """Test internal fragment with charge"""
        annot = PafAnnotation.make_internal(1, 4, charge=2)
        assert annot.charge == 2
        assert "^2" in annot.serialize()


class TestMakeImmonium:
    """Test make_immonium factory method"""

    def test_basic_immonium(self):
        """Test creating a basic immonium ion"""
        annot = PafAnnotation.make_immonium("K")
        assert isinstance(annot.ion_type, ImmoniumIon)
        assert annot.ion_type.amino_acid == AminoAcids.K
        assert annot.serialize() == "IK"

    def test_immonium_with_amino_acid_enum(self):
        """Test using AminoAcids enum directly"""
        annot = PafAnnotation.make_immonium(AminoAcids.W)
        assert isinstance(annot.ion_type, ImmoniumIon)
        assert annot.ion_type.amino_acid == AminoAcids.W
        assert annot.serialize() == "IW"

    def test_immonium_with_modification(self):
        """Test immonium ion with modification"""
        annot = PafAnnotation.make_immonium("M", modification="Oxidation")
        assert isinstance(annot.ion_type, ImmoniumIon)
        assert annot.ion_type.modification == "Oxidation"
        assert annot.serialize() == "IM[Oxidation]"

    def test_immonium_with_charge(self):
        """Test immonium ion with charge state"""
        annot = PafAnnotation.make_immonium("H", charge=2)
        assert annot.charge == 2
        assert "^2" in annot.serialize()


class TestMakeReference:
    """Test make_reference factory method"""

    def test_basic_reference_ion(self):
        """Test creating a basic reference ion"""
        annot = PafAnnotation.make_reference("Phospho")
        assert isinstance(annot.ion_type, ReferenceIon)
        assert annot.ion_type.name == "Phospho"
        assert annot.serialize() == "r[Phospho]"

    def test_reference_with_charge(self):
        """Test reference ion with charge"""
        annot = PafAnnotation.make_reference("iTRAQ4plex", charge=2)
        assert annot.charge == 2
        assert "^2" in annot.serialize()

    def test_reference_with_neutral_loss(self):
        """Test reference ion with neutral loss"""
        annot = PafAnnotation.make_reference("TMT6plex", neutral_losses=["-H2O"])
        assert len(annot.neutral_losses) == 1


class TestMakeNamedCompound:
    """Test make_named_compound factory method"""

    def test_basic_named_compound(self):
        """Test creating a basic named compound"""
        annot = PafAnnotation.make_named_compound("Urocanic Acid")
        assert isinstance(annot.ion_type, NamedCompound)
        assert annot.ion_type.name == "Urocanic Acid"
        assert annot.serialize() == "_{Urocanic Acid}"

    def test_named_compound_with_charge(self):
        """Test named compound with charge"""
        annot = PafAnnotation.make_named_compound("Some Compound", charge=3)
        assert annot.charge == 3
        assert "^3" in annot.serialize()


class TestMakeFormula:
    """Test make_formula factory method"""

    def test_basic_formula(self):
        """Test creating a basic formula annotation"""
        annot = PafAnnotation.make_formula("C13H9")
        assert isinstance(annot.ion_type, ChemicalFormula)
        assert annot.ion_type.formula == "C13H9"
        assert annot.serialize() == "f{C13H9}"

    def test_formula_with_mass_error(self):
        """Test formula with mass error"""
        annot = PafAnnotation.make_formula("C12H9N", mass_error=-0.55, mass_error_unit="ppm")
        assert annot.mass_error is not None
        assert annot.mass_error.value == -0.55
        assert annot.mass_error.unit == "ppm"
        assert "/-0.55ppm" in annot.serialize()

    def test_formula_with_adduct(self):
        """Test formula with adduct"""
        annot = PafAnnotation.make_formula("C14H10N", adducts=["+H"])
        assert len(annot.adducts) == 1
        assert "[M+H]" in annot.serialize()


class TestMakeSmiles:
    """Test make_smiles factory method"""

    def test_basic_smiles(self):
        """Test creating a basic SMILES annotation"""
        annot = PafAnnotation.make_smiles("CN=C=O")
        assert isinstance(annot.ion_type, SMILESCompound)
        assert annot.ion_type.smiles == "CN=C=O"
        assert annot.serialize() == "s{CN=C=O}"

    def test_smiles_with_adduct(self):
        """Test SMILES with adduct"""
        annot = PafAnnotation.make_smiles("COc(c1)cccc1C#N", adducts=["+H", "+Na"], charge=2)
        assert len(annot.adducts) == 2
        assert annot.charge == 2
        assert "[M+H+Na]^2" in annot.serialize()


class TestMakeUnknown:
    """Test make_unknown factory method"""

    def test_basic_unknown(self):
        """Test creating a basic unknown ion"""
        annot = PafAnnotation.make_unknown()
        assert isinstance(annot.ion_type, UnknownIon)
        assert annot.ion_type.label is None
        assert annot.serialize() == "?"

    def test_unknown_with_label(self):
        """Test unknown ion with label"""
        annot = PafAnnotation.make_unknown(label=5)
        assert isinstance(annot.ion_type, UnknownIon)
        assert annot.ion_type.label == 5
        assert annot.serialize() == "?5"

    def test_unknown_with_charge(self):
        """Test unknown ion with charge"""
        annot = PafAnnotation.make_unknown(label=3, charge=2)
        assert annot.charge == 2
        assert "^2" in annot.serialize()


class TestMakeMethodsWithComplexModifiers:
    """Test make methods with multiple complex modifiers"""

    def test_peptide_fully_annotated(self):
        """Test peptide ion with all modifiers"""
        annot = PafAnnotation.make_peptide(
            "y",
            7,
            sequence="PEPTIDE",
            neutral_losses=["-H2O", "-NH3"],
            isotopes=["+i13C", 3],
            charge=2,
            mass_error=0.05,
            mass_error_unit="da",
            confidence=0.99,
            is_auxiliary=True,
            analyte_reference=1,
        )
        serialized = annot.serialize()
        assert serialized.startswith("&1@y7{PEPTIDE}")
        assert "-H2O" in serialized
        assert "-NH3" in serialized
        assert "+i13C" in serialized
        assert "+3i" in serialized
        assert "^2" in serialized
        assert "/0.05" in serialized
        assert "*0.99" in serialized

    def test_precursor_with_multiple_adducts(self):
        """Test precursor with multiple adducts"""
        annot = PafAnnotation.make_precursor(adducts=["+2H", "+Na"], charge=3)
        assert len(annot.adducts) == 2
        assert annot.charge == 3
        assert "[M+2H+Na]^3" in annot.serialize()

    def test_string_to_object_conversion(self):
        """Test that strings are properly converted to objects"""
        annot = PafAnnotation.make_peptide(
            "b",
            5,
            neutral_losses=["-H2O"],  # String
            isotopes=["+i"],  # String
            adducts=["+Na"],  # String
        )
        # Verify they were converted to proper objects
        assert isinstance(annot.neutral_losses[0], NeutralLoss)
        assert isinstance(annot.isotopes[0], IsotopeSpecification)
        assert isinstance(annot.adducts[0], Adduct)

    def test_object_direct_usage(self):
        """Test using objects directly instead of strings"""
        annot = PafAnnotation.make_peptide(
            "b",
            5,
            neutral_losses=[NeutralLoss(count=-1, base_formula="H2O")],
            isotopes=[IsotopeSpecification(count=1, element="13C")],
            adducts=[Adduct(count=1, base_formula="Na")],
        )
        assert len(annot.neutral_losses) == 1
        assert len(annot.isotopes) == 1
        assert len(annot.adducts) == 1


class TestMakeMethodsRoundTrip:
    """Test that make methods produce parseable annotations"""

    def test_precursor_roundtrip(self):
        """Test precursor can be serialized and parsed back"""
        annot = PafAnnotation.make_precursor(charge=2, neutral_losses=["-H2O"])
        serialized = annot.serialize()
        parsed = PafAnnotation.parse(serialized)
        assert isinstance(parsed.ion_type, PrecursorIon)
        assert parsed.charge == 2
        assert len(parsed.neutral_losses) == 1

    def test_peptide_roundtrip(self):
        """Test peptide can be serialized and parsed back"""
        annot = PafAnnotation.make_peptide("y", 10, sequence="PEPTIDE", charge=2)
        serialized = annot.serialize()
        parsed = PafAnnotation.parse(serialized)
        assert isinstance(parsed.ion_type, PeptideIon)
        assert parsed.ion_type.series == IonSeries.Y
        assert parsed.ion_type.position == 10
        assert parsed.ion_type.sequence == "PEPTIDE"
        assert parsed.charge == 2

    def test_immonium_roundtrip(self):
        """Test immonium can be serialized and parsed back"""
        annot = PafAnnotation.make_immonium("K", modification="Acetyl")
        serialized = annot.serialize()
        parsed = PafAnnotation.parse(serialized)
        assert isinstance(parsed.ion_type, ImmoniumIon)
        assert parsed.ion_type.amino_acid == AminoAcids.K
        assert parsed.ion_type.modification == "Acetyl"

    def test_formula_roundtrip(self):
        """Test formula can be serialized and parsed back"""
        annot = PafAnnotation.make_formula("C13H9", mass_error=0.55, mass_error_unit="ppm")
        serialized = annot.serialize()
        parsed = PafAnnotation.parse(serialized)
        assert isinstance(parsed.ion_type, ChemicalFormula)
        assert parsed.ion_type.formula == "C13H9"
        assert parsed.mass_error is not None
        assert parsed.mass_error.value == 0.55
        assert parsed.mass_error.unit == "ppm"
