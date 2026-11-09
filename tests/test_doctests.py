"""Execute every docstring example in the package."""

import doctest
import importlib
import importlib.util
import pkgutil

import paftacular

# Measured 2026-11-08. Raise these when you add modules or examples; never lower them
# without saying why in the commit message.
MIN_MODULES = 22
MIN_EXAMPLES = 0
# The walk discovers all 22 names in every environment. Three of them import a third-party
# package that only the mcp extra installs, so without that extra they cannot be imported at
# all. They are skipped by name, and only when the package they need is genuinely absent, so
# a broken import with the extra installed still fails this test.
NEEDS_EXTRA = {
    "paftacular.mcp.handlers": "pydantic",
    "paftacular.mcp.models": "pydantic",
    "paftacular.mcp.server": "mcp",
}
# With every mcp module skipped the remaining 19 must still import.
MIN_IMPORTED = MIN_MODULES - len(NEEDS_EXTRA)


def _module_names():
    names = [paftacular.__name__]
    names.extend(name for _finder, name, _ispkg in pkgutil.walk_packages(paftacular.__path__, prefix=f"{paftacular.__name__}."))
    return names


def test_docstring_examples():
    names = _module_names()
    assert len(names) >= MIN_MODULES, f"module walk found only {len(names)}: {names}"
    assert not set(NEEDS_EXTRA) - set(names), f"NEEDS_EXTRA names a module the walk does not find: {sorted(set(NEEDS_EXTRA) - set(names))}"
    mods = []
    for name in names:
        required = NEEDS_EXTRA.get(name)
        if required is not None and importlib.util.find_spec(required) is None:
            continue
        mods.append(importlib.import_module(name))
    assert len(mods) >= MIN_IMPORTED, f"only {len(mods)} of {len(names)} modules were imported"
    attempted = failed = 0
    for mod in mods:
        result = doctest.testmod(mod, verbose=False)
        attempted += result.attempted
        failed += result.failed
    assert failed == 0, f"{failed} of {attempted} docstring examples failed"
    assert attempted >= MIN_EXAMPLES, f"only {attempted} examples ran, expected >= {MIN_EXAMPLES}"
