from paftacular.comps import IsotopeSpecification


def test_isotope_serialize_positive_generic():
    iso = IsotopeSpecification(count=1, element=None, is_average=False)
    assert str(iso) == "+i"


def test_isotope_serialize_negative_element():
    iso = IsotopeSpecification(count=-2, element="13C", is_average=False)
    assert str(iso) == "-2i13C"


def test_isotope_serialize_average():
    iso = IsotopeSpecification(count=1, element=None, is_average=True)
    assert str(iso) == "+iA"
