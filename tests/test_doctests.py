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


def _walk_error(name: str) -> None:
    # pkgutil.walk_packages catches ImportError and keeps going when no onerror is passed, so a
    # subpackage whose __init__ breaks is still yielded while its submodules silently vanish. The
    # MIN_MODULES floor below already catches that drop. This callback only adds the name of the
    # broken subpackage to the failure, instead of a bare count.
    raise AssertionError(f"walking into {name} failed, so its submodules were not seen")


def _module_names():
    names = [paftacular.__name__]
    names.extend(name for _finder, name, _ispkg in pkgutil.walk_packages(paftacular.__path__, prefix=f"{paftacular.__name__}.", onerror=_walk_error))
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


def test_needs_extra_gate_agrees_with_reality():
    """Check every NEEDS_EXTRA entry against what importing that module actually does.

    NEEDS_EXTRA is a declaration, and nothing else in the suite checks it against what happens. A
    key naming the wrong root makes the gate skip a module in an environment where it imports fine,
    or claim the extra is present where it is absent, and MIN_IMPORTED cannot see either because it
    is a subtraction. So compute the gate's verdict the way the runner does, then import the module
    for real and compare, per module rather than over the subtree.

    The per-module form rests on one invariant: each entry's root is a root that module genuinely
    needs, directly or through its own relative imports. That is what makes "root absent implies the
    import fails" true. It is not the same as the module's own module-scope import list, and here it
    is not equal to it: paftacular.mcp.server is keyed to mcp while its module-scope roots are mcp
    and pydantic. The entry is sound anyway, because the mcp distribution requires pydantic, so mcp
    present implies pydantic present and mcp absent makes the import fail.

    The form breaks for one entry shared across several modules whose roots differ, because a
    supported half-installed environment then makes that entry right for one module and wrong for
    another. Keep one entry per module.

    Do not assert that every gated module fails: whether a sibling fails depends on whether the
    gated subtree's __init__ imports the failing root, which varies by package and by how much of an
    extra is installed.
    """
    # Both halves of the mechanism can go inert without any assertion noticing. An empty mapping
    # gates nothing, and an entry with no root reads as present in every environment.
    assert NEEDS_EXTRA, "the gate maps no modules, so it can never skip or fail anything"
    for name, required in sorted(NEEDS_EXTRA.items()):
        assert required, f"{name} is listed with no root to look for, so the gate reads as present everywhere"
        predicted_absent = importlib.util.find_spec(required) is None
        failures = []
        try:
            importlib.import_module(name)
        except ImportError as error:
            failures.append(error)
        # Checked before the verdict so that a failure naming an unlisted root is reported as the
        # broken module it is, rather than as a disagreement about the extra.
        for error in failures:
            root = (getattr(error, "name", None) or "").split(".")[0]
            unnamed = root or "an import the error does not name"
            assert root == required, f"{name} failed on {unnamed}, which its key does not list, so it is broken rather than gated: {error}"
        # Ordered before the equality so that a verdict stuck at absent, which skips the subtree
        # while every floor still passes, is reported as that and not as a count.
        assert not (predicted_absent and not failures), f"the gate says {required} is absent, but {name} imported, so the key is wrong or the verdict is stuck"
        verdict = "absent" if predicted_absent else "present"
        assert predicted_absent == bool(failures), f"the gate says {required} is {verdict}, but importing {name} gave {failures}"
