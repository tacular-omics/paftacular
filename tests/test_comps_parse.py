"""Tests for parse methods in comps.py classes"""

import pytest

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


class TestMassErrorParse:
    """Test MassError.parse method"""

    def test_parse_ppm(self):
        """Test parsing ppm mass error"""
        me = MassError.parse("0.55ppm")
        assert me.value == 0.55
        assert me.unit == "ppm"

    def test_parse_da(self):
        """Test parsing Da mass error"""
        me = MassError.parse("0.06")
        assert me.value == 0.06
        assert me.unit == "da"

    def test_parse_negative_ppm(self):
        """Test parsing negative ppm"""
        me = MassError.parse("-2.5ppm")
        assert me.value == -2.5
        assert me.unit == "ppm"

    def test_parse_with_whitespace(self):
        """Test parsing with whitespace"""
        me = MassError.parse("  1.23ppm  ")
        assert me.value == 1.23
        assert me.unit == "ppm"


class TestIsotopeSpecificationParse:
    """Test IsotopeSpecification.parse method"""

    def test_parse_basic_positive(self):
        """Test parsing basic positive isotope"""
        iso = IsotopeSpecification.parse("+i")
        assert iso.count == 1
        assert iso.element is None
        assert iso.is_average is False

    def test_parse_basic_negative(self):
        """Test parsing basic negative isotope"""
        iso = IsotopeSpecification.parse("-i")
        assert iso.count == -1
        assert iso.element is None

    def test_parse_with_count(self):
        """Test parsing isotope with count"""
        iso = IsotopeSpecification.parse("+2i")
        assert iso.count == 2

    def test_parse_with_element(self):
        """Test parsing isotope with element"""
        iso = IsotopeSpecification.parse("+i13C")
        assert iso.count == 1
        assert iso.element == "13C"

    def test_parse_negative_with_element(self):
        """Test parsing negative isotope with element"""
        iso = IsotopeSpecification.parse("-2i13C")
        assert iso.count == -2
        assert iso.element == "13C"

    def test_parse_average(self):
        """Test parsing average isotopomer"""
        iso = IsotopeSpecification.parse("+iA")
        assert iso.count == 1
        assert iso.is_average is True

    def test_parse_negative_average(self):
        """Test parsing negative average"""
        iso = IsotopeSpecification.parse("-2iA")
        assert iso.count == -2
        assert iso.is_average is True

    def test_parse_no_sign(self):
        """Test parsing without explicit sign (defaults to positive)"""
        iso = IsotopeSpecification.parse("i")
        assert iso.count == 1

    def test_parse_without_element(self):
        """Test parsing isotope without element specified"""
        iso = IsotopeSpecification.parse("i")
        assert iso.count == 1
        assert iso.element is None


class TestPeptideIonParse:
    """Test PeptideIon.parse method"""

    def test_parse_basic_b_ion(self):
        """Test parsing basic b-ion"""
        ion = PeptideIon.parse("b5")
        assert ion.series == IonSeries.B
        assert ion.position == 5
        assert ion.sequence is None

    def test_parse_y_ion(self):
        """Test parsing y-ion"""
        ion = PeptideIon.parse("y10")
        assert ion.series == IonSeries.Y
        assert ion.position == 10

    def test_parse_with_sequence(self):
        """Test parsing ion with sequence"""
        ion = PeptideIon.parse("b3{PEPTIDE}")
        assert ion.series == IonSeries.B
        assert ion.position == 3
        assert ion.sequence == "PEPTIDE"

    def test_parse_c_ion(self):
        """Test parsing c-ion"""
        ion = PeptideIon.parse("c7")
        assert ion.series == IonSeries.C
        assert ion.position == 7

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid peptide ion"):
            PeptideIon.parse("invalid")


class TestInternalFragmentParse:
    """Test InternalFragment.parse method"""

    def test_parse_basic(self):
        """Test parsing basic internal fragment"""
        frag = InternalFragment.parse("m5:10")
        assert frag.start_position == 5
        assert frag.end_position == 10
        assert frag.sequence is None

    def test_parse_with_sequence(self):
        """Test parsing internal fragment with sequence"""
        frag = InternalFragment.parse("m2:5{TIDE}")
        assert frag.start_position == 2
        assert frag.end_position == 5
        assert frag.sequence == "TIDE"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid internal fragment"):
            InternalFragment.parse("invalid")


class TestImmoniumIonParse:
    """Test ImmoniumIon.parse method"""

    def test_parse_basic(self):
        """Test parsing basic immonium ion"""
        ion = ImmoniumIon.parse("IK")
        assert ion.amino_acid == AminoAcids.K
        assert ion.modification is None

    def test_parse_with_modification(self):
        """Test parsing immonium ion with modification"""
        ion = ImmoniumIon.parse("IM[Oxidation]")
        assert ion.amino_acid == AminoAcids.M
        assert ion.modification == "Oxidation"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid immonium ion"):
            ImmoniumIon.parse("invalid")


class TestReferenceIonParse:
    """Test ReferenceIon.parse method"""

    def test_parse_basic(self):
        """Test parsing basic reference ion"""
        ion = ReferenceIon.parse("r[Phospho]")
        assert ion.name == "Phospho"

    def test_parse_complex_name(self):
        """Test parsing reference with complex name"""
        ion = ReferenceIon.parse("r[iTRAQ4plex]")
        assert ion.name == "iTRAQ4plex"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid reference ion"):
            ReferenceIon.parse("invalid")


class TestNamedCompoundParse:
    """Test NamedCompound.parse method"""

    def test_parse_basic(self):
        """Test parsing basic named compound"""
        comp = NamedCompound.parse("_{Urocanic Acid}")
        assert comp.name == "Urocanic Acid"

    def test_parse_simple_name(self):
        """Test parsing compound with simple name"""
        comp = NamedCompound.parse("_{Water}")
        assert comp.name == "Water"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid named compound"):
            NamedCompound.parse("invalid")


class TestChemicalFormulaParse:
    """Test ChemicalFormula.parse method"""

    def test_parse_basic(self):
        """Test parsing basic formula"""
        formula = ChemicalFormula.parse("f{C13H9}")
        assert formula.formula == "C13H9"

    def test_parse_complex(self):
        """Test parsing complex formula"""
        formula = ChemicalFormula.parse("f{C12H9N}")
        assert formula.formula == "C12H9N"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid chemical formula"):
            ChemicalFormula.parse("invalid")


class TestSMILESCompoundParse:
    """Test SMILESCompound.parse method"""

    def test_parse_basic(self):
        """Test parsing basic SMILES"""
        smiles = SMILESCompound.parse("s{CN=C=O}")
        assert smiles.smiles == "CN=C=O"

    def test_parse_complex(self):
        """Test parsing complex SMILES"""
        smiles = SMILESCompound.parse("s{COc(c1)cccc1C#N}")
        assert smiles.smiles == "COc(c1)cccc1C#N"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid SMILES compound"):
            SMILESCompound.parse("invalid")


class TestUnknownIonParse:
    """Test UnknownIon.parse method"""

    def test_parse_without_label(self):
        """Test parsing unknown ion without label"""
        ion = UnknownIon.parse("?")
        assert ion.label is None

    def test_parse_with_label(self):
        """Test parsing unknown ion with label"""
        ion = UnknownIon.parse("?5")
        assert ion.label == 5

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid unknown ion"):
            UnknownIon.parse("invalid")


class TestPrecursorIonParse:
    """Test PrecursorIon.parse method"""

    def test_parse_basic(self):
        """Test parsing precursor ion"""
        ion = PrecursorIon.parse("p")
        assert isinstance(ion, PrecursorIon)

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid precursor ion"):
            PrecursorIon.parse("invalid")

    def test_parse_with_extra(self):
        """Test parsing with extra characters fails"""
        with pytest.raises(ValueError, match="Invalid precursor ion"):
            PrecursorIon.parse("p2")


class TestAdductParse:
    """Test Adduct.parse method"""

    def test_parse_positive_single(self):
        """Test parsing positive single adduct"""
        adduct = Adduct.parse("+H")
        assert adduct.count == 1
        assert adduct.base_formula == "H"

    def test_parse_positive_multiple(self):
        """Test parsing positive adduct with count"""
        adduct = Adduct.parse("+2Na")
        assert adduct.count == 2
        assert adduct.base_formula == "Na"

    def test_parse_negative(self):
        """Test parsing negative adduct"""
        adduct = Adduct.parse("-H")
        assert adduct.count == -1
        assert adduct.base_formula == "H"

    def test_parse_negative_multiple(self):
        """Test parsing negative adduct with count"""
        adduct = Adduct.parse("-2H")
        assert adduct.count == -2
        assert adduct.base_formula == "H"

    def test_parse_complex_formula(self):
        """Test parsing adduct with complex formula"""
        adduct = Adduct.parse("+NH4")
        assert adduct.count == 1
        assert adduct.base_formula == "NH4"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError, match="Invalid adduct"):
            Adduct.parse("invalid")


class TestNeutralLossParse:
    """Test NeutralLoss.parse method"""

    def test_parse_formula_loss(self):
        """Test parsing formula-based loss"""
        loss = NeutralLoss.parse("-H2O")
        assert loss.count == -1
        assert loss.base_formula == "H2O"

    def test_parse_formula_gain(self):
        """Test parsing formula-based gain"""
        loss = NeutralLoss.parse("+NH3")
        assert loss.count == 1
        assert loss.base_formula == "NH3"

    def test_parse_multiple_loss(self):
        """Test parsing multiple neutral loss"""
        loss = NeutralLoss.parse("-2H2O")
        assert loss.count == -2
        assert loss.base_formula == "H2O"

    def test_parse_mass_loss(self):
        """Test parsing mass-based loss"""
        loss = NeutralLoss.parse("-98.5")
        assert loss.count == -1
        assert loss.base_mass == 98.5

    def test_parse_reference_loss(self):
        """Test parsing reference-based loss"""
        loss = NeutralLoss.parse("-[Phospho]")
        assert loss.count == -1
        assert loss.base_reference == "Phospho"

    def test_parse_reference_with_count(self):
        """Test parsing reference loss with count"""
        loss = NeutralLoss.parse("-2[Phospho]")
        assert loss.count == -2
        assert loss.base_reference == "Phospho"

    def test_parse_invalid(self):
        """Test parsing invalid string"""
        with pytest.raises(ValueError):
            NeutralLoss.parse("invalid")


class TestParseRoundtrip:
    """Test that parse methods are inverse of serialize"""

    def test_mass_error_roundtrip(self):
        """Test MassError parse/serialize roundtrip"""
        original = MassError(0.55, "ppm")
        parsed = MassError.parse(original.serialize())
        assert parsed.value == original.value
        assert parsed.unit == original.unit

    def test_isotope_roundtrip(self):
        """Test IsotopeSpecification parse/serialize roundtrip"""
        original = IsotopeSpecification(count=2, element="13C")
        parsed = IsotopeSpecification.parse(original.serialize())
        assert parsed.count == original.count
        assert parsed.element == original.element

    def test_peptide_ion_roundtrip(self):
        """Test PeptideIon parse/serialize roundtrip"""
        original = PeptideIon(series=IonSeries.B, position=5, sequence="PEPTIDE")
        parsed = PeptideIon.parse(original.serialize())
        assert parsed.series == original.series
        assert parsed.position == original.position
        assert parsed.sequence == original.sequence

    def test_adduct_roundtrip(self):
        """Test Adduct parse/serialize roundtrip"""
        original = Adduct(count=2, base_formula="Na")
        parsed = Adduct.parse(original.serialize())
        assert parsed.count == original.count
        assert parsed.base_formula == original.base_formula

    def test_neutral_loss_roundtrip(self):
        """Test NeutralLoss parse/serialize roundtrip"""
        original = NeutralLoss(count=-1, base_formula="H2O")
        parsed = NeutralLoss.parse(original.serialize())
        assert parsed.count == original.count
        assert parsed.base_formula == original.base_formula
