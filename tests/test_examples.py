import re
from pathlib import Path

import pytest


def test_readme_python_examples():
    pytest.importorskip("peptacular")
    root = Path(__file__).resolve().parents[1]
    namespace = {"__name__": "__example__"}
    for code in re.findall(r"```python\n(.*?)```", (root / "README.md").read_text(), flags=re.DOTALL):
        exec(compile(code, "README.md", "exec"), namespace)


def test_usage_code_blocks():
    pytest.importorskip("peptacular")
    pytest.importorskip("pysmiles")
    root = Path(__file__).resolve().parents[1]
    namespace = {"__name__": "__example__"}
    lines = (root / "docs/usage.rst").read_text().splitlines()
    index = 0
    while index < len(lines):
        if lines[index] not in (".. code-block:: python", ".. testcode::"):
            index += 1
            continue
        index += 1
        code = []
        while index < len(lines) and (not lines[index].strip() or lines[index].startswith("   ")):
            code.append(lines[index][3:])
            index += 1
        exec(compile("\n".join(code), "docs/usage.rst", "exec"), namespace)
