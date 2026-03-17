import peptacular as pt
import pytest

import paftacular as paf


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


if __name__ == "__main__":
    import pytest

    pytest.main([__file__])
