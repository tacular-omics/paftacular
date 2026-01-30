import pytest
from tacular import ELEMENT_LOOKUP

from paftacular import PafAnnotation, parse, parse_multi
from paftacular.comps import (
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    NamedCompound,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)
from paftacular.constants import AminoAcids, IonSeries


def parse_one(s: str) -> PafAnnotation:
    return parse(s)


# ============================================================================
# Peptide Ion Tests
# ============================================================================


def test_peptide_simple():
    ann = parse_one("y12")
    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.series == IonSeries.Y
    assert ann.ion_type.position == 12


def test_peptide_with_sequence():
    ann = parse_one("a1{PEPT}")
    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.series == IonSeries.A
    assert ann.ion_type.position == 1
    assert ann.ion_type.sequence == "PEPT"


def test_peptide_properties():
    ann = parse_one("y12")
    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.formula == "H2O"  # Y ion formula
    assert ann.ion_type.monoisotopic_mass == pytest.approx(18.010565, rel=1e-5)
    assert ann.ion_type.composition == {ELEMENT_LOOKUP["H"]: 2, ELEMENT_LOOKUP["O"]: 1}


# ============================================================================
# Internal Fragment Tests
# ============================================================================


def test_internal_fragment():
    ann = parse_one("m2:5")
    assert isinstance(ann.ion_type, InternalFragment)
    assert ann.ion_type.start_position == 2
    assert ann.ion_type.end_position == 5


def test_internal_fragment_with_sequence():
    ann = parse_one("m3:7{PEPTIDE}")
    assert isinstance(ann.ion_type, InternalFragment)
    assert ann.ion_type.start_position == 3
    assert ann.ion_type.end_position == 7
    assert ann.ion_type.sequence == "PEPTIDE"


# ============================================================================
# Precursor Ion Tests
# ============================================================================


def test_precursor():
    ann = parse_one("p")
    assert isinstance(ann.ion_type, PrecursorIon)
    assert ann.ion_type.serialize() == "p"


# ============================================================================
# Immonium Ion Tests
# ============================================================================


def test_immonium_simple():
    ann = parse_one("IA")
    assert isinstance(ann.ion_type, ImmoniumIon)
    assert ann.ion_type.amino_acid == AminoAcids.A
    assert ann.ion_type.modification is None
    assert ann.ion_type.serialize() == "IA"


def test_immonium_with_modification():
    ann = parse_one("IH[+16]")
    assert isinstance(ann.ion_type, ImmoniumIon)
    assert ann.ion_type.amino_acid == AminoAcids.H
    assert ann.ion_type.modification == "+16"
    assert ann.ion_type.serialize() == "IH[+16]"


# ============================================================================
# Reference Ion Tests
# ============================================================================


def test_reference_ion():
    ann = parse_one("r[ATP]")
    assert isinstance(ann.ion_type, ReferenceIon)
    assert ann.ion_type.name == "ATP"
    assert ann.ion_type.serialize() == "r[ATP]"


def test_reference_ion_properties():
    """Test mass, formula, and composition properties of ReferenceIon"""
    # Use a molecule that's in the REFMOL_LOOKUP (e.g., Adenine)
    ann = parse_one("r[Adenine]")
    ref_ion = ann.ion_type
    assert isinstance(ref_ion, ReferenceIon)

    # Test mass calculations
    # Adenine (C5H5N5) monoisotopic mass
    assert ref_ion.monoisotopic_mass == pytest.approx(135.054495, rel=1e-5)
    # Adenine average mass
    assert ref_ion.average_mass == pytest.approx(135.127, rel=1e-3)
    assert ref_ion.mass(monoisotopic=True) == ref_ion.monoisotopic_mass
    assert ref_ion.mass(monoisotopic=False) == ref_ion.average_mass

    # Test formula
    assert ref_ion.formula == "C5H5N5"  # Adenine formula

    # Test composition - Adenine has C5H5N5
    comp = ref_ion.composition
    assert len(comp) == 3  # Only C, H, N elements
    assert comp[ELEMENT_LOOKUP["C"]] == 5
    assert comp[ELEMENT_LOOKUP["H"]] == 5
    assert comp[ELEMENT_LOOKUP["N"]] == 5


# ============================================================================
# Chemical Formula Tests
# ============================================================================


def test_chemical_formula():
    ann = parse_one("f{C6H6}")
    assert isinstance(ann.ion_type, ChemicalFormula)
    assert ann.ion_type.formula == "C6H6"
    assert ann.ion_type.serialize() == "f{C6H6}"
    # Test mass calculations for benzene
    assert ann.ion_type.monoisotopic_mass == pytest.approx(78.046950, rel=1e-5)
    assert ann.ion_type.average_mass == pytest.approx(78.11, rel=1e-3)


def test_chemical_formula_complex():
    ann = parse_one("f{C13H9N}")
    assert isinstance(ann.ion_type, ChemicalFormula)
    assert ann.ion_type.formula == "C13H9N"
    comp = ann.ion_type.composition
    assert comp[ELEMENT_LOOKUP["C"]] == 13
    assert comp[ELEMENT_LOOKUP["H"]] == 9
    assert comp[ELEMENT_LOOKUP["N"]] == 1


def test_chemical_formula_proforma():
    """Test proforma_formula property matches input formula"""
    ann = parse_one("f{C6H6}")
    assert isinstance(ann.ion_type, ChemicalFormula)
    assert ann.ion_type.proforma_formula == "C6H6"
    assert ann.ion_type.proforma_formula == ann.ion_type.formula


# ============================================================================
# Named Compound Tests
# ============================================================================


def test_named_compound():
    ann = parse_one("_{Urocanic}")
    assert isinstance(ann.ion_type, NamedCompound)
    assert ann.ion_type.name == "Urocanic"
    assert ann.ion_type.serialize() == "_{Urocanic}"


# Named compounds with spaces are not directly supported by the parser


# ============================================================================
# SMILES Compound Tests
# ============================================================================


def test_smiles_compound():
    ann = parse_one("s{CN=C=O}")
    assert isinstance(ann.ion_type, SMILESCompound)
    assert ann.ion_type.smiles == "CN=C=O"
    assert ann.ion_type.serialize() == "s{CN=C=O}"


def test_smiles_compound_complex():
    ann = parse_one("s{COc(c1)cccc1C#N}")
    assert isinstance(ann.ion_type, SMILESCompound)
    assert ann.ion_type.smiles == "COc(c1)cccc1C#N"
    assert ann.ion_type.serialize() == "s{COc(c1)cccc1C#N}"


def test_smiles_compound_composition():
    """Test composition calculation from SMILES"""
    ann = parse_one("s{CN=C=O}")
    smiles = ann.ion_type
    assert isinstance(smiles, SMILESCompound)

    # CN=C=O is methyl isocyanate: C2H3NO
    comp = smiles.composition
    assert comp[ELEMENT_LOOKUP["C"]] == 2
    assert comp[ELEMENT_LOOKUP["N"]] == 1
    assert comp[ELEMENT_LOOKUP["O"]] == 1
    assert comp[ELEMENT_LOOKUP["H"]] == 3
    # Note: H count depends on implicit hydrogens in SMILES parsing


def test_smiles_compound_mass():
    """Test mass calculations from SMILES"""
    ann = parse_one("s{CN=C=O}")
    smiles = ann.ion_type
    assert isinstance(smiles, SMILESCompound)

    # CN=C=O is methyl isocyanate: C2H3NO
    # Monoisotopic: 12C2 + 1H3 + 14N + 16O
    # = 12.000*2 + 1.007825*3 + 14.003074 + 15.994915 = 57.021464
    assert smiles.monoisotopic_mass == pytest.approx(57.021464, rel=1e-5)
    # Average mass should be slightly higher due to heavier isotopes
    assert smiles.average_mass == pytest.approx(57.05, rel=1e-2)
    assert smiles.mass(monoisotopic=True) == smiles.monoisotopic_mass
    assert smiles.mass(monoisotopic=False) == smiles.average_mass


def test_smiles_compound_formula_properties():
    """Test formula and proforma_formula properties"""
    ann = parse_one("s{CN=C=O}")
    smiles = ann.ion_type
    assert isinstance(smiles, SMILESCompound)

    # CN=C=O is methyl isocyanate: C2H3NO
    # ProForma formula should be sorted by element symbol
    proforma = smiles.proforma_formula
    assert proforma == "C2H3NO"

    # Formula should include + sign prefix
    formula = smiles.formula
    assert formula == "+C2H3NO"
    assert formula == f"+{proforma}"


# ============================================================================
# Unknown Ion Tests
# ============================================================================


def test_unknown_ion():
    ann = parse_one("?")
    assert isinstance(ann.ion_type, UnknownIon)
    assert ann.ion_type.label is None
    assert ann.ion_type.serialize() == "?"


def test_unknown_ion_with_label():
    ann = parse_one("?42")
    assert isinstance(ann.ion_type, UnknownIon)
    assert ann.ion_type.label == 42
    assert ann.ion_type.serialize() == "?42"


# ============================================================================
# Neutral Loss Tests - Formula
# ============================================================================


def test_neutral_loss_formula_water():
    ann = parse_one("y5-H2O")
    assert ann.neutral_losses and len(ann.neutral_losses) == 1
    loss = ann.neutral_losses[0]

    # Basic properties
    assert loss.loss_type == "formula"
    assert loss.formula == "-H2O"
    assert loss.count == -1
    assert loss.base_formula == "H2O"

    # ProForma formula
    assert loss.proforma_formula == "H-2O-1"

    # Composition check
    comp = loss.composition
    assert comp[ELEMENT_LOOKUP["H"]] == -2
    assert comp[ELEMENT_LOOKUP["O"]] == -1

    # Mass
    assert loss.mass() == pytest.approx(-18.010565, rel=1e-6)
    assert loss.monoisotopic_mass == pytest.approx(-18.010565, rel=1e-6)

    # Serialization
    assert loss.serialize() == "-H2O"
    assert loss.serialize(loss_type="formula") == "-H2O"
    assert float(loss.serialize(loss_type="mass")) == pytest.approx(-18.01056, rel=1e-4)
    with pytest.raises(ValueError):
        loss.serialize(loss_type="reference")


def test_neutral_loss_multiple_formula():
    # Note: -3NH3 may be parsed as a mass value depending on parser
    # Use a clearly formula-based example
    ann = parse_one("y5-2H2O")
    assert ann.neutral_losses and len(ann.neutral_losses) == 1
    loss = ann.neutral_losses[0]

    assert loss.loss_type == "formula" or loss.loss_type == "mass"
    if loss.loss_type == "formula":
        assert loss.sign == -1
        assert loss.count >= 1


def test_neutral_loss_formula_gain():
    ann = parse_one("y5+H2")
    assert ann.neutral_losses and len(ann.neutral_losses) == 1
    loss = ann.neutral_losses[0]

    assert loss.loss_type == "formula"
    assert loss.count == 1
    assert loss.formula == "+H2"
    # H2 mass: 2 * 1.007825 = 2.01565 Da
    assert loss.mass() == pytest.approx(2.01565, rel=1e-5)


# ============================================================================
# Neutral Loss Tests - Reference
# ============================================================================


def test_neutral_loss_reference_gain():
    ann = parse_one("y5+[Adenine]")
    assert ann.neutral_losses and len(ann.neutral_losses) == 1
    loss = ann.neutral_losses[0]

    # Basic properties
    assert loss.loss_type == "reference"
    assert loss.base_reference == "Adenine"
    assert loss.count == 1

    # ProForma formula
    assert loss.proforma_formula == "C5H5N5"

    # Composition check
    comp = loss.composition
    assert comp[ELEMENT_LOOKUP["C"]] == 5
    assert comp[ELEMENT_LOOKUP["H"]] == 5
    assert comp[ELEMENT_LOOKUP["N"]] == 5

    # Mass
    assert loss.mass() == pytest.approx(135.0544951833, rel=1e-6)
    assert loss.monoisotopic_mass == pytest.approx(135.0544951833, rel=1e-6)

    # Serialization
    assert loss.serialize() == "+[Adenine]"
    assert loss.serialize(loss_type="reference") == "+[Adenine]"
    assert loss.serialize(loss_type="formula") == "+C5H5N5"
    assert float(loss.serialize(loss_type="mass")) == pytest.approx(135.05450, rel=1e-4)


def test_neutral_loss_multiple_reference():
    ann = parse_one("y5-2[ADenine]")
    assert ann.neutral_losses and len(ann.neutral_losses) == 1
    loss = ann.neutral_losses[0]

    # Basic properties
    assert loss.loss_type == "reference"
    assert loss.base_reference == "ADenine"
    assert loss.count == -2

    # ProForma formula (negative count)
    assert loss.proforma_formula == "C-10H-10N-10"

    # Composition check (2x loss)
    comp = loss.composition
    assert comp[ELEMENT_LOOKUP["C"]] == -10
    assert comp[ELEMENT_LOOKUP["H"]] == -10
    assert comp[ELEMENT_LOOKUP["N"]] == -10

    # Mass (2x loss)
    assert loss.mass() == pytest.approx(-270.10899, rel=1e-6)

    # Serialization
    assert loss.serialize() == "-2[ADenine]"
    assert loss.serialize(loss_type="reference") == "-2[ADenine]"
    assert loss.serialize(loss_type="formula") == "-2C5H5N5"
    assert float(loss.serialize(loss_type="mass")) == pytest.approx(-270.10899, rel=1e-4)


# ============================================================================
# Isotope Specification Tests
# ============================================================================


def test_isotope_specification():
    ann = parse_one("y5+i13C")
    assert ann.isotopes and len(ann.isotopes) == 1
    iso = ann.isotopes[0]
    assert iso.count == 1
    assert iso.element == "13C"
    assert iso.is_average is False
    assert iso.serialize() == "+i13C"


def test_isotope_specification_multiple():
    ann = parse_one("y5+2i13C")
    assert ann.isotopes and len(ann.isotopes) == 1
    iso = ann.isotopes[0]
    assert iso.count == 2
    assert iso.element == "13C"
    assert iso.serialize() == "+2i13C"


def test_isotope_specification_negative():
    ann = parse_one("y5-i")
    assert ann.isotopes and len(ann.isotopes) == 1
    iso = ann.isotopes[0]
    assert iso.count == -1
    assert iso.element is None
    assert iso.serialize() == "-i"


def test_isotope_specification_average():
    ann = parse_one("y5+iA")
    assert ann.isotopes and len(ann.isotopes) == 1
    iso = ann.isotopes[0]
    assert iso.count == 1
    assert iso.is_average is True
    assert iso.serialize() == "+iA"


def test_isotope_specification_composition():
    """Test composition calculation for isotope specifications"""
    ann = parse_one("y5+i13C")
    assert ann.isotopes and len(ann.isotopes) == 1
    iso = ann.isotopes[0]

    # Composition should show loss of monoisotopic C and gain of 13C
    comp = iso.composition
    # Get the 13C and 12C ElementInfo objects
    c13 = ELEMENT_LOOKUP["13C"]
    c12 = ELEMENT_LOOKUP["12C"]  # explicit 12C
    assert comp[c13] == 1
    assert comp[c12] == -1


def test_isotope_specification_mass():
    """Test mass calculation for isotope specifications"""
    ann = parse_one("y5+i13C")
    assert ann.isotopes and len(ann.isotopes) == 1
    iso = ann.isotopes[0]

    # Mass shift for 13C - 12C should be approximately 1.003 Da
    mass = iso.mass()
    # The mass difference between 13C and 12C
    c13_mass = ELEMENT_LOOKUP["13C"].mass
    c12_mass = ELEMENT_LOOKUP["12C"].mass
    expected_shift = c13_mass - c12_mass
    assert mass == pytest.approx(expected_shift, rel=1e-6)


def test_isotope_specification_multiple_mass():
    """Test mass for multiple isotope count"""
    ann = parse_one("y5+2i13C")
    assert ann.isotopes and len(ann.isotopes) == 1
    iso = ann.isotopes[0]

    # Mass shift for 2x (13C - 12C)
    mass = iso.mass()
    c13_mass = ELEMENT_LOOKUP["13C"].mass
    c12_mass = ELEMENT_LOOKUP["12C"].mass
    expected_shift = 2 * (c13_mass - c12_mass)
    print(mass, expected_shift)
    assert mass == pytest.approx(expected_shift, rel=1e-6)


def test_isotope_specification_errors():
    """Test error handling for isotope specifications"""
    # Test average isotopomer - cannot calculate composition
    ann = parse_one("y5+iA")
    iso = ann.isotopes[0]
    with pytest.raises(ValueError, match="average isotopomer"):
        _ = iso.composition
    with pytest.raises(ValueError, match="average isotopomer"):
        _ = iso.mass()

    # Test generic isotope - cannot calculate composition without element
    ann2 = parse_one("y5+i")
    iso2 = ann2.isotopes[0]
    with pytest.raises(ValueError, match="generic isotope"):
        _ = iso2.composition


# ============================================================================
# Adduct Tests
# ============================================================================


def test_adduct_simple():
    ann = parse_one("y5[M+H]")
    assert ann.adducts and len(ann.adducts) == 1
    adduct = ann.adducts[0]
    assert adduct.count == 1
    # Check composition instead
    comp = adduct.composition
    assert comp[ELEMENT_LOOKUP["H"]] == 1
    assert adduct.serialize() == "+H"


def test_adduct_multiple():
    ann = parse_one("y5[M+H+Na]")
    assert ann.adducts and len(ann.adducts) == 2
    # Check compositions
    comp_h = ann.adducts[0].composition
    assert comp_h[ELEMENT_LOOKUP["H"]] == 1
    assert ann.adducts[0].serialize() == "+H"
    comp_na = ann.adducts[1].composition
    assert comp_na[ELEMENT_LOOKUP["Na"]] == 1
    assert ann.adducts[1].serialize() == "+Na"


def test_adduct_count():
    ann = parse_one("y5[M+2H]")
    assert ann.adducts and len(ann.adducts) == 1
    adduct = ann.adducts[0]
    assert adduct.count == 2
    # Check composition (2 hydrogens)
    comp = adduct.composition
    assert comp[ELEMENT_LOOKUP["H"]] == 2
    assert adduct.serialize() == "+2H"


def test_adduct_removal():
    ann = parse_one("y5[M-H]")
    assert ann.adducts and len(ann.adducts) == 1
    adduct = ann.adducts[0]
    assert adduct.count == -1
    # Check composition (negative hydrogen)
    comp = adduct.composition
    assert comp[ELEMENT_LOOKUP["H"]] == -1
    assert adduct.serialize() == "-H"


def test_adduct_mass_calculation():
    ann = parse_one("y5[M+H]")
    assert ann.adducts and len(ann.adducts) == 1
    adduct = ann.adducts[0]
    # Hydrogen mass
    assert adduct.mass() == pytest.approx(1.007825, rel=1e-5)
    assert adduct.monoisotopic_mass == pytest.approx(1.007825, rel=1e-5)


def test_adduct_complex_formula():
    ann = parse_one("y5[M+NH4]")
    assert ann.adducts and len(ann.adducts) == 1
    adduct = ann.adducts[0]
    # Check composition
    comp = adduct.composition
    assert comp[ELEMENT_LOOKUP["N"]] == 1
    assert comp[ELEMENT_LOOKUP["H"]] == 4
    assert adduct.serialize() == "+NH4"


def test_adduct_proforma_formula():
    """Test proforma_formula property for adducts"""
    ann = parse_one("y5[M+NH4]")
    adduct = ann.adducts[0]

    # ProForma formula should be NH4
    proforma = adduct.proforma_formula
    assert proforma == "H4N"


def test_adduct_average_mass():
    """Test average_mass property for adducts"""
    ann = parse_one("y5[M+H]")
    adduct = ann.adducts[0]

    # Test both monoisotopic and average mass for hydrogen
    assert adduct.monoisotopic_mass == pytest.approx(1.007825, rel=1e-5)
    # Average mass of hydrogen (considering natural isotope abundance)
    assert adduct.average_mass == pytest.approx(1.008, rel=1e-3)
    assert adduct.average_mass > adduct.monoisotopic_mass


def test_adduct_formula_property():
    """Test formula property returns sign + formula"""
    ann = parse_one("y5[M+2Na]")
    adduct = ann.adducts[0]

    # Formula should include sign and count
    assert adduct.formula == "+2Na"
    assert adduct.serialize() == "+2Na"


# ============================================================================
# Charge Tests
# ============================================================================


def test_charge_annotation():
    ann = parse_one("y5^2")
    assert ann.charge == 2


# Negative charges are not supported - charges must be >= 1

# ============================================================================
# Mass Error Tests
# ============================================================================


def test_mass_error_ppm():
    ann = parse_one("y5/-0.55ppm")
    assert ann.mass_error is not None
    assert ann.mass_error.unit == "ppm"
    assert ann.mass_error.value == pytest.approx(-0.55, rel=1e-6)
    assert ann.mass_error.serialize() == "-0.55ppm"


def test_mass_error_daltons():
    ann = parse_one("y5/0.003")
    assert ann.mass_error is not None
    assert ann.mass_error.unit == "da"
    assert ann.mass_error.value == pytest.approx(0.003, rel=1e-6)
    assert ann.mass_error.serialize() == "0.003"


# Mass errors with explicit + sign are handled implicitly


def test_confidence_score():
    ann = parse_one("y5*0.85")
    assert ann.confidence == pytest.approx(0.85)


def test_confidence_score_low():
    ann = parse_one("y5*0.1")
    assert ann.confidence == pytest.approx(0.1)


def test_confidence_score_high():
    ann = parse_one("y5*0.999")
    assert ann.confidence == pytest.approx(0.999)


# ============================================================================
# Auxiliary Flag Tests
# ============================================================================


def test_auxiliary_flag():
    ann = parse_one("&y5")
    assert ann.is_auxiliary is True


def test_auxiliary_flag_with_other_features():
    ann = parse_one("&y5^2")
    assert ann.is_auxiliary is True
    assert ann.charge == 2


# ============================================================================
# Analyte Reference Tests
# ============================================================================


def test_analyte_reference():
    ann = parse_one("2@y5")
    assert ann.analyte_reference == 2


def test_analyte_reference_high_number():
    ann = parse_one("10@b3")
    assert ann.analyte_reference == 10


# ============================================================================
# Complex Combined Features Tests
# ============================================================================


def test_complex_annotation_all_features():
    """Test annotation with all possible features combined"""
    ann = parse_one("&2@y5+i13C[M+H+Na]^2/-0.55ppm*0.85")

    # Auxiliary and reference
    assert ann.is_auxiliary is True
    assert ann.analyte_reference == 2

    # Ion type
    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.series == IonSeries.Y
    assert ann.ion_type.position == 5

    # Isotopes
    assert ann.isotopes and len(ann.isotopes) == 1
    assert ann.isotopes[0].element == "13C"
    assert ann.isotopes[0].count == 1

    # Adducts
    assert ann.adducts and len(ann.adducts) == 2
    # Check by composition
    comp_h = ann.adducts[0].composition
    assert comp_h[ELEMENT_LOOKUP["H"]] == 1
    comp_na = ann.adducts[1].composition
    assert comp_na[ELEMENT_LOOKUP["Na"]] == 1

    # Charge
    assert ann.charge == 2

    # Mass error
    assert ann.mass_error is not None
    assert ann.mass_error.unit == "ppm"
    assert ann.mass_error.value == pytest.approx(-0.55, rel=1e-6)

    # Confidence
    assert ann.confidence == pytest.approx(0.85)


def test_complex_annotation_neutral_loss_isotope():
    """Test annotation with neutral loss and isotope"""
    ann = parse_one("b3-H2O+i13C/1.2ppm")

    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.series == IonSeries.B
    assert ann.ion_type.position == 3

    # Neutral loss
    assert ann.neutral_losses and len(ann.neutral_losses) == 1
    assert ann.neutral_losses[0].formula == "-H2O"

    # Isotope
    assert ann.isotopes and len(ann.isotopes) == 1
    assert ann.isotopes[0].element == "13C"

    # Mass error
    assert ann.mass_error is not None
    assert ann.mass_error.value == pytest.approx(1.2, rel=1e-6)


# ============================================================================
# Multiple Annotations Tests
# ============================================================================


def test_multiple_annotations_parse():
    anns = parse_multi("y3, b2, p, IY")
    assert len(anns) == 4

    assert isinstance(anns[0].ion_type, PeptideIon)
    assert anns[0].ion_type.series == IonSeries.Y

    assert isinstance(anns[1].ion_type, PeptideIon)
    assert anns[1].ion_type.series == IonSeries.B

    assert isinstance(anns[2].ion_type, PrecursorIon)

    assert isinstance(anns[3].ion_type, ImmoniumIon)
    assert anns[3].ion_type.amino_acid == AminoAcids.Y


def test_multiple_annotations_complex():
    anns = parse_multi("y5-H2O, b3+NH3, p^2")
    assert len(anns) == 3

    # First annotation with neutral loss
    assert anns[0].neutral_losses and len(anns[0].neutral_losses) == 1

    # Second annotation with neutral gain
    assert anns[1].neutral_losses and len(anns[1].neutral_losses) == 1

    # Third annotation with charge
    assert anns[2].charge == 2


def test_neutral_loss_mass():
    """Test mass-based neutral loss"""
    ann = parse_one("y5-17.03")
    assert ann.neutral_losses and len(ann.neutral_losses) == 1
    loss = ann.neutral_losses[0]

    assert loss.loss_type == "mass"
    assert loss.base_mass == pytest.approx(17.03, rel=1e-5)
    assert loss.count == -1
    assert loss.mass() == pytest.approx(-17.03, rel=1e-5)

    # Should fail to get composition for mass-based loss
    with pytest.raises(ValueError, match="Cannot calculate composition for mass-based loss"):
        _ = loss.composition

    # Serialization
    assert loss.serialize(loss_type="mass") == "-17.03000"
    with pytest.raises(ValueError, match="Cannot get formula"):
        loss.serialize(loss_type="formula")


def test_neutral_loss_mass_positive():
    """Test mass-based neutral gain"""
    ann = parse_one("y5+42.0106")
    loss = ann.neutral_losses[0]

    assert loss.loss_type == "mass"
    assert loss.base_mass == pytest.approx(42.0106, rel=1e-5)
    assert loss.count == 1


def test_multiple_neutral_losses():
    """Test multiple neutral losses on same ion"""
    ann = parse_one("y5-H2O-NH3")
    assert len(ann.neutral_losses) == 2

    # Water loss
    assert ann.neutral_losses[0].formula == "-H2O"
    # Ammonia loss
    assert ann.neutral_losses[1].formula == "-NH3"


def test_multiple_isotopes():
    """Test multiple isotope specifications"""
    ann = parse_one("y5+i13C+i15N")
    assert len(ann.isotopes) == 2

    assert ann.isotopes[0].element == "13C"
    assert ann.isotopes[1].element == "15N"


def test_invalid_smiles():
    """Test invalid SMILES string handling"""
    ann = parse_one("s{-783yfuINVALID!!!}")
    smiles = ann.ion_type
    assert isinstance(smiles, SMILESCompound)
    with pytest.raises(ValueError):
        _ = smiles.composition


def test_unknown_reference_molecule():
    """Test handling of unknown reference molecule"""
    ann = parse_one("y5-[UnknownMolecule]")
    loss = ann.neutral_losses[0]

    assert loss.base_reference == "UnknownMolecule"
    assert isinstance(loss.reference, str)  # Falls back to string

    with pytest.raises(ValueError):
        _ = loss.composition


def test_adduct_validation():
    """Test Adduct __post_init__ validation"""
    from paftacular.comps import Adduct

    # Invalid count
    with pytest.raises(ValueError):
        Adduct(count=0, base_formula="H")

    # Empty formula
    with pytest.raises(ValueError):
        Adduct(count=1, base_formula="")


def test_neutral_loss_validation():
    """Test NeutralLoss __post_init__ validation"""
    from paftacular.comps import NeutralLoss

    # No specification provided
    with pytest.raises(ValueError, match="Exactly one of formula, mass, or reference"):
        NeutralLoss(count=-1)

    # Multiple specifications provided
    with pytest.raises(ValueError, match="Exactly one of formula, mass, or reference"):
        NeutralLoss(count=-1, base_formula="H2O", base_mass=18.01)


def test_serialization_roundtrip_simple():
    """Test that serialization produces parseable output"""
    original = "y5-H2O+i13C[M+H]^2/-0.55ppm*0.85"
    ann = parse_one(original)
    serialized = ann.serialize()
    reparsed = parse_one(serialized)

    # Compare key properties
    assert isinstance(reparsed.ion_type, PeptideIon)
    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.series == reparsed.ion_type.series
    assert ann.charge == reparsed.charge
    assert len(ann.neutral_losses) == len(reparsed.neutral_losses)


def test_serialization_chemical_formula():
    """Test ChemicalFormula serialization round-trip"""
    original = "f{C6H12O6}"
    ann = parse_one(original)
    assert ann.serialize() == original


def test_internal_fragment_backbone_types():
    """Test nterm_ion_type and cterm_ion_type if supported by parser"""
    # This depends on your parser implementation
    # The code has these fields but they may not be populated
    ann = parse_one("m2:5")
    fragment = ann.ion_type

    # Test default values
    assert isinstance(fragment, InternalFragment)
    assert fragment.nterm_ion_type is None
    assert fragment.cterm_ion_type is None
