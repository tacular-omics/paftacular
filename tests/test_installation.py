import importlib.util
import os
from pathlib import Path

import pytest

import paftacular as p


def test_installed_artifact_location():
    if os.environ.get("PAFTACULAR_TEST_WHEEL") == "1":
        assert not Path(p.__file__).resolve().is_relative_to(Path(__file__).resolve().parents[1] / "src")


def test_optional_installation():
    expected = os.environ.get("PAFTACULAR_EXPECT_EXTRAS")
    for package, extra in (("peptacular", "peptacular"), ("pysmiles", "smiles")):
        available = importlib.util.find_spec(package) is not None
        if expected is not None:
            assert available == (expected in (extra, "all"))
        text = "y2{PE}" if extra == "peptacular" else "s{C}"
        annotation = p.parse_single(text)
        assert annotation.serialize() == text
        if available:
            assert annotation.mass() > 0
        else:
            with pytest.raises(ImportError, match=package):
                annotation.mass()


def test_core_operations_need_no_extras():
    annotation = p.parse_single("f{C2H4}+i")
    assert annotation.formula() == "C[13C]H4"
    assert p.PafAnnotation.from_dict(annotation.to_dict()) == annotation
    assert [result.ok for result in p.parse_batch(["y2", "bad"])] == [True, False]
