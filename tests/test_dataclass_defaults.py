"""Guard against descriptor objects leaking into dataclass field defaults.

A read-only attribute declared on a mixin as a real ``property`` is picked up as the field
default of any dataclass subclass that redeclares it, because ``dataclasses`` resolves a
default with ``getattr(cls, name, MISSING)``. When a non-default field follows, the class
fails loudly at import. When it does not, as in the ``NeutralLoss`` shape where every later
field is keyword-only with a default, nothing raises at class creation or construction and
the ``property`` object itself becomes the value. Mixin declarations that exist only for the
type checker therefore live behind ``if TYPE_CHECKING`` (see ``comps/base.py``). This test
walks the package so the guard also covers dataclasses nobody has written yet.
"""

import dataclasses
import importlib
import inspect
import pkgutil
import types
from dataclasses import KW_ONLY, dataclass

import paftacular

# Slot descriptors are the one benign case: a slots dataclass inheriting from another can pick
# one up through the same getattr, and dataclasses treats it as "no default".
ALLOWED_DEFAULT_DESCRIPTORS = (types.MemberDescriptorType,)


def _descriptor_defaults(cls):
    """Report fields of one dataclass whose default is a descriptor object."""
    offenders = []
    for field in dataclasses.fields(cls):
        default = field.default
        if default is dataclasses.MISSING or isinstance(default, ALLOWED_DEFAULT_DESCRIPTORS):
            continue
        if hasattr(type(default), "__get__"):
            offenders.append(f"{cls.__module__}.{cls.__qualname__}.{field.name} defaults to a {type(default).__name__}")
    return offenders


def _package_modules():
    """Import every paftacular module, skipping those whose optional extra is absent."""
    modules = [paftacular]
    for info in pkgutil.walk_packages(paftacular.__path__, prefix=f"{paftacular.__name__}."):
        try:
            modules.append(importlib.import_module(info.name))
        except ImportError:
            continue
    return modules


def _package_dataclasses():
    """Collect every dataclass the package defines, keyed by identity to drop re-exports.

    Reading module attributes rather than ``__subclasses__()`` matters: ``slots=True`` rebuilds
    the class, so the discarded pre-slots original stays registered as a subclass while only the
    live class is bound to its module.
    """
    found = {}
    for module in _package_modules():
        for obj in vars(module).values():
            if inspect.isclass(obj) and dataclasses.is_dataclass(obj) and obj.__module__.startswith(paftacular.__name__):
                found[obj] = obj
    return list(found)


def test_no_descriptor_field_defaults():
    classes = _package_dataclasses()
    assert len(classes) >= 10, f"the module walk collected only {len(classes)} dataclasses, so it is not covering the package"

    offenders = [entry for cls in classes for entry in _descriptor_defaults(cls)]
    detail = ", ".join(offenders)
    assert not offenders, (
        "descriptor objects became dataclass field defaults, most likely a mixin "
        f"attribute declared as a real property instead of behind TYPE_CHECKING: {detail}"
    )


def test_guard_detects_a_property_default():
    """The silent shape the guard exists for, so a passing run above means something."""

    class Mixin:
        @property
        def count(self) -> int:
            raise NotImplementedError

    @dataclass(frozen=True, slots=True)
    class SilentShape(Mixin):
        count: int
        _: KW_ONLY
        base_formula: str | None = None

    # Nothing raised at class creation, and count is now optional with the property as its value.
    assert isinstance(SilentShape().count, property)
    expected = f"{__name__}.test_guard_detects_a_property_default.<locals>.SilentShape.count defaults to a property"
    assert _descriptor_defaults(SilentShape) == [expected]
