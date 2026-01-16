import pytest

from paftacular import PafAnnotation, mzPAFParser
from paftacular.comps import (
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    NamedCompound,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)
from paftacular.constants import AminoAcids, IonSeries


def parse_one(s: str) -> PafAnnotation:
    return mzPAFParser().parse_single(s)


def test_peptide_simple():
    ann = parse_one("y12")
    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.series == IonSeries.Y
    assert ann.ion_type.position == 12


def test_peptide_with_sequence():
    ann = parse_one("a1{PEPT}")
    assert isinstance(ann.ion_type, PeptideIon)
    assert ann.ion_type.series == IonSeries.A
    assert ann.ion_type.sequence == "PEPT"


def test_internal_fragment():
    ann = parse_one("m2:5")
    assert isinstance(ann.ion_type, InternalFragment)
    assert ann.ion_type.start_position == 2
    assert ann.ion_type.end_position == 5


def test_precursor():
    ann = parse_one("p")
    assert isinstance(ann.ion_type, PrecursorIon)


def test_immonium_simple_and_mod():
    ann = parse_one("IA")
    assert isinstance(ann.ion_type, ImmoniumIon)
    assert ann.ion_type.amino_acid == AminoAcids.A

    ann2 = parse_one("IH[+16]")
    assert isinstance(ann2.ion_type, ImmoniumIon)
    assert ann2.ion_type.amino_acid == AminoAcids.H
    assert ann2.ion_type.modification == "+16"


def test_reference_named_formula_smiles_unknown():
    ann_r = parse_one("r[ATP]")
    assert isinstance(ann_r.ion_type, ReferenceIon)
    assert ann_r.ion_type.name == "ATP"

    ann_f = parse_one("f{C6H6}")
    assert isinstance(ann_f.ion_type, ChemicalFormula)
    assert ann_f.ion_type.formula == "C6H6"

    ann_n = parse_one("_{Urocanic}")
    assert isinstance(ann_n.ion_type, NamedCompound)
    assert ann_n.ion_type.name == "Urocanic"

    ann_s = parse_one("s{CN=C=O}")
    assert isinstance(ann_s.ion_type, SMILESCompound)
    assert ann_s.ion_type.smiles == "CN=C=O"

    ann_u = parse_one("?")
    assert isinstance(ann_u.ion_type, UnknownIon)


def test_neutral_losses_and_reference_loss():
    ann = parse_one("y5-H2O")
    assert ann.neutral_losses and isinstance(ann.neutral_losses[0], NeutralLoss)
    assert ann.neutral_losses[0].loss_type == "formula"
    assert ann.neutral_losses[0].formula == "H2O"

    ann2 = parse_one("y5+[Phospho]")
    assert ann2.neutral_losses and ann2.neutral_losses[0].loss_type == "reference"
    assert ann2.neutral_losses[0]._reference == "Phospho"


def test_isotopes_adducts_charge_mass_confidence_ref_aux():
    ann = parse_one("&2@y5+i13C[M+H+Na]^2/-0.55ppm*0.85")
    assert ann.is_auxiliary is True
    assert ann.analyte_reference == 2
    assert ann.charge == 2
    assert ann.mass_error is not None and ann.mass_error.unit == "ppm"
    assert pytest.approx(ann.mass_error.value, rel=1e-6) == 0.55
    assert ann.confidence == pytest.approx(0.85)
    # isotopes
    assert ann.isotopes and ann.isotopes[0].element == "13C"
    # adducts
    assert ann.adducts and len(ann.adducts) >= 1


def test_multiple_annotations_parse():
    anns = mzPAFParser().parse("y3, b2, p, IY")
    assert len(anns) == 4
