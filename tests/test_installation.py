import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

import paftacular as p


def test_installed_artifact_location():
    if os.environ.get("PAFTACULAR_TEST_WHEEL") == "1":
        assert not Path(p.__file__).resolve().is_relative_to(Path(__file__).resolve().parents[1] / "src")


def test_optional_installation():
    expected = os.environ.get("PAFTACULAR_EXPECT_EXTRAS")
    selected = set(expected.split(",")) if expected else set()
    if "all" in selected:
        selected.update(("mcp", "peptacular", "smiles"))
    if "mcp" in selected:
        selected.add("peptacular")
    for package, extra in (("peptacular", "peptacular"), ("pysmiles", "smiles")):
        available = importlib.util.find_spec(package) is not None
        if expected is not None:
            assert available == (extra in selected)
        text = "y2{PE}" if extra == "peptacular" else "s{C}"
        annotation = p.parse_single(text)
        assert annotation.serialize() == text
        if available:
            assert annotation.mass() > 0
        else:
            with pytest.raises(ImportError, match=package):
                annotation.mass()

    if expected is not None:
        assert (importlib.util.find_spec("mcp") is not None) == ("mcp" in selected)


def test_mcp_import_boundary_and_cli():
    code = "import sys\nimport paftacular\nimport paftacular.mcp\nassert 'mcp' not in sys.modules\nassert 'pydantic' not in sys.modules\n"
    subprocess.run([sys.executable, "-I", "-c", code], check=True, capture_output=True, timeout=15)
    result = subprocess.run([sys.executable, "-I", "-m", "paftacular.mcp", "--version"], capture_output=True, text=True, timeout=15)
    assert result.returncode == 0
    assert result.stdout.strip() == f"paftacular-mcp {p.__version__}"
    if importlib.util.find_spec("mcp") is None:
        result = subprocess.run([sys.executable, "-I", "-m", "paftacular.mcp"], input="", capture_output=True, text=True, timeout=15)
        assert result.returncode == 2
        assert "paftacular[mcp]" in result.stderr
        assert not result.stdout


def test_core_operations_need_no_extras():
    annotation = p.parse_single("f{C2H4}+i")
    assert annotation.formula() == "C[13C]H4"
    assert p.PafAnnotation.from_dict(annotation.to_dict()) == annotation
    assert [result.ok for result in p.parse_batch(["y2", "bad"])] == [True, False]
