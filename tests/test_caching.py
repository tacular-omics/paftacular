"""Tests for object caching in frozen dataclasses"""

from paftacular.comps import (
    Adduct,
    ImmoniumIon,
    IsotopeSpecification,
    NamedCompound,
    NeutralLoss,
    PrecursorIon,
    ReferenceIon,
    UnknownIon,
)
from paftacular.constants import MAX_CACHE_SIZE, AminoAcids


class TestIsotopeSpecificationCaching:
    """Test that IsotopeSpecification instances are cached"""

    def test_same_isotope_returns_same_object(self):
        """Test that creating identical isotopes returns the same object"""
        iso1 = IsotopeSpecification(count=1, element="13C")
        iso2 = IsotopeSpecification(count=1, element="13C")
        assert iso1 is iso2

    def test_different_count_returns_different_object(self):
        """Test that different counts return different objects"""
        iso1 = IsotopeSpecification(count=1, element="13C")
        iso2 = IsotopeSpecification(count=2, element="13C")
        assert iso1 is not iso2

    def test_different_element_returns_different_object(self):
        """Test that different elements return different objects"""
        iso1 = IsotopeSpecification(count=1, element="13C")
        iso2 = IsotopeSpecification(count=1, element="15N")
        assert iso1 is not iso2

    def test_average_vs_not_average_returns_different_object(self):
        """Test that average isotopes are cached separately"""
        iso1 = IsotopeSpecification(count=1, is_average=True)
        iso2 = IsotopeSpecification(count=1, is_average=False)
        assert iso1 is not iso2

    def test_parsed_isotope_uses_cache(self):
        """Test that parsed isotopes use the cache"""
        iso1 = IsotopeSpecification(count=1, element="13C")
        iso2 = IsotopeSpecification.parse("+i13C")
        assert iso1 is iso2


class TestNeutralLossCaching:
    """Test that NeutralLoss instances are cached"""

    def test_same_loss_returns_same_object(self):
        """Test that creating identical losses returns the same object"""
        loss1 = NeutralLoss(count=-1, base_formula="H2O")
        loss2 = NeutralLoss(count=-1, base_formula="H2O")
        assert loss1 is loss2

    def test_different_count_returns_different_object(self):
        """Test that different counts return different objects"""
        loss1 = NeutralLoss(count=-1, base_formula="H2O")
        loss2 = NeutralLoss(count=-2, base_formula="H2O")
        assert loss1 is not loss2

    def test_different_formula_returns_different_object(self):
        """Test that different formulas return different objects"""
        loss1 = NeutralLoss(count=-1, base_formula="H2O")
        loss2 = NeutralLoss(count=-1, base_formula="NH3")
        assert loss1 is not loss2

    def test_mass_based_loss_cached(self):
        """Test that mass-based losses are cached"""
        loss1 = NeutralLoss(count=-1, base_mass=18.01)
        loss2 = NeutralLoss(count=-1, base_mass=18.01)
        assert loss1 is loss2

    def test_reference_based_loss_cached(self):
        """Test that reference-based losses are cached"""
        loss1 = NeutralLoss(count=-1, base_reference="Phospho")
        loss2 = NeutralLoss(count=-1, base_reference="Phospho")
        assert loss1 is loss2

    def test_parsed_loss_uses_cache(self):
        """Test that parsed losses use the cache"""
        loss1 = NeutralLoss(count=-1, base_formula="H2O")
        loss2 = NeutralLoss.parse("-H2O")
        assert loss1 is loss2


class TestAdductCaching:
    """Test that Adduct instances are cached"""

    def test_same_adduct_returns_same_object(self):
        """Test that creating identical adducts returns the same object"""
        adduct1 = Adduct(count=1, base_formula="H")
        adduct2 = Adduct(count=1, base_formula="H")
        assert adduct1 is adduct2

    def test_different_count_returns_different_object(self):
        """Test that different counts return different objects"""
        adduct1 = Adduct(count=1, base_formula="H")
        adduct2 = Adduct(count=2, base_formula="H")
        assert adduct1 is not adduct2

    def test_different_formula_returns_different_object(self):
        """Test that different formulas return different objects"""
        adduct1 = Adduct(count=1, base_formula="H")
        adduct2 = Adduct(count=1, base_formula="Na")
        assert adduct1 is not adduct2

    def test_negative_count_cached(self):
        """Test that negative adducts are cached"""
        adduct1 = Adduct(count=-1, base_formula="H")
        adduct2 = Adduct(count=-1, base_formula="H")
        assert adduct1 is adduct2

    def test_parsed_adduct_uses_cache(self):
        """Test that parsed adducts use the cache"""
        adduct1 = Adduct(count=1, base_formula="Na")
        adduct2 = Adduct.parse("+Na")
        assert adduct1 is adduct2


class TestCachingInParsing:
    """Test that caching works correctly when parsing annotations"""

    def test_multiple_annotations_share_neutral_loss(self):
        """Test that parsing multiple annotations with same loss shares objects"""
        from paftacular import parse_multi

        annots = parse_multi("b5-H2O, y10-H2O")
        assert len(annots) == 2
        assert len(annots[0].neutral_losses) == 1
        assert len(annots[1].neutral_losses) == 1
        # Should be the same cached object
        assert annots[0].neutral_losses[0] is annots[1].neutral_losses[0]

    def test_multiple_annotations_share_isotope(self):
        """Test that parsing multiple annotations with same isotope shares objects"""
        from paftacular import parse_multi

        annots = parse_multi("b5+i13C, y10+i13C")
        assert len(annots) == 2
        assert len(annots[0].isotopes) == 1
        assert len(annots[1].isotopes) == 1
        # Should be the same cached object
        assert annots[0].isotopes[0] is annots[1].isotopes[0]

    def test_multiple_annotations_share_adduct(self):
        """Test that parsing multiple annotations with same adduct shares objects"""
        from paftacular import parse_multi

        annots = parse_multi("p[M+Na], p[M+Na]")
        assert len(annots) == 2
        assert len(annots[0].adducts) == 1
        assert len(annots[1].adducts) == 1
        # Should be the same cached object
        assert annots[0].adducts[0] is annots[1].adducts[0]


class TestCacheMemoryBehavior:
    """Test memory behavior of caching"""

    def test_cache_reduces_memory_for_repeated_objects(self):
        """Test that cache reduces memory usage for repeated objects"""
        # Create 1000 identical isotopes
        isotopes = [IsotopeSpecification(count=1, element="13C") for _ in range(1000)]

        # All should be the same object
        first = isotopes[0]
        assert all(iso is first for iso in isotopes)

    def test_cache_size_grows_with_unique_objects(self):
        """Test that cache grows with unique objects"""
        initial_size = len(IsotopeSpecification._cache)

        # Create 10 unique isotopes with high counts to avoid collision with existing cache
        for i in range(100, 110):
            IsotopeSpecification(count=i, element="13C")

        # Cache should have grown by at least 10
        assert len(IsotopeSpecification._cache) >= initial_size + 10

    def test_cache_eviction_when_full(self):
        """Test that cache evicts oldest entry when MAX_CACHE_SIZE is reached"""
        # Clear the cache for this test
        UnknownIon._cache.clear()

        # Fill cache to MAX_CACHE_SIZE
        for i in range(MAX_CACHE_SIZE):
            UnknownIon(label=i)

        # Verify cache is at capacity
        assert len(UnknownIon._cache) == MAX_CACHE_SIZE

        # First entry should still be in cache
        first_key = (0,)
        assert first_key in UnknownIon._cache
        first_ion = UnknownIon._cache[first_key]

        # Add one more entry, which should trigger eviction
        new_ion = UnknownIon(label=MAX_CACHE_SIZE)

        # Cache should not exceed MAX_CACHE_SIZE
        assert len(UnknownIon._cache) == MAX_CACHE_SIZE

        # New entry should be in cache
        new_key = (MAX_CACHE_SIZE,)
        assert new_key in UnknownIon._cache

        # First entry should have been evicted
        assert first_key not in UnknownIon._cache

        # Creating the first entry again should create a new object (not the same as before)
        recreated_first = UnknownIon(label=0)
        assert recreated_first is not first_ion


class TestImmoniumIonCaching:
    """Test that ImmoniumIon instances are cached"""

    def test_same_immonium_returns_same_object(self):
        """Test that creating identical immonium ions returns the same object"""
        ion1 = ImmoniumIon(amino_acid=AminoAcids.A)
        ion2 = ImmoniumIon(amino_acid=AminoAcids.A)
        assert ion1 is ion2

    def test_different_amino_acid_returns_different_object(self):
        """Test that different amino acids return different objects"""
        ion1 = ImmoniumIon(amino_acid=AminoAcids.A)
        ion2 = ImmoniumIon(amino_acid=AminoAcids.G)
        assert ion1 is not ion2

    def test_with_modification_returns_different_object(self):
        """Test that modified immonium ions are cached separately"""
        ion1 = ImmoniumIon(amino_acid=AminoAcids.A)
        ion2 = ImmoniumIon(amino_acid=AminoAcids.A, modification="Oxidation")
        assert ion1 is not ion2

    def test_same_modification_returns_same_object(self):
        """Test that identical modified ions return same object"""
        ion1 = ImmoniumIon(amino_acid=AminoAcids.A, modification="Oxidation")
        ion2 = ImmoniumIon(amino_acid=AminoAcids.A, modification="Oxidation")
        assert ion1 is ion2


class TestReferenceIonCaching:
    """Test that ReferenceIon instances are cached"""

    def test_same_reference_returns_same_object(self):
        """Test that creating identical reference ions returns the same object"""
        ion1 = ReferenceIon(name="w")
        ion2 = ReferenceIon(name="w")
        assert ion1 is ion2

    def test_different_name_returns_different_object(self):
        """Test that different names return different objects"""
        ion1 = ReferenceIon(name="w")
        ion2 = ReferenceIon(name="d")
        assert ion1 is not ion2


class TestNamedCompoundCaching:
    """Test that NamedCompound instances are cached"""

    def test_same_compound_returns_same_object(self):
        """Test that creating identical named compounds returns the same object"""
        comp1 = NamedCompound(name="H2O")
        comp2 = NamedCompound(name="H2O")
        assert comp1 is comp2

    def test_different_name_returns_different_object(self):
        """Test that different names return different objects"""
        comp1 = NamedCompound(name="H2O")
        comp2 = NamedCompound(name="NH3")
        assert comp1 is not comp2


class TestUnknownIonCaching:
    """Test that UnknownIon instances are cached"""

    def test_same_unknown_returns_same_object(self):
        """Test that creating identical unknown ions returns the same object"""
        ion1 = UnknownIon(label=1)
        ion2 = UnknownIon(label=1)
        assert ion1 is ion2

    def test_different_label_returns_different_object(self):
        """Test that different labels return different objects"""
        ion1 = UnknownIon(label=1)
        ion2 = UnknownIon(label=2)
        assert ion1 is not ion2


class TestPrecursorIonCaching:
    """Test that PrecursorIon is a singleton"""

    def test_precursor_always_returns_same_object(self):
        """Test that PrecursorIon is effectively a singleton"""
        ion1 = PrecursorIon()
        ion2 = PrecursorIon()
        assert ion1 is ion2

    def test_multiple_precursors_are_same(self):
        """Test that creating many precursors all return the same object"""
        ions = [PrecursorIon() for _ in range(100)]
        first = ions[0]
        assert all(ion is first for ion in ions)
