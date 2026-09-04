"""Version 1 interchange data, independent of the legacy as_dict() display."""

from collections.abc import Mapping
from dataclasses import asdict, fields
from enum import Enum
from types import UnionType
from typing import Any, Literal, cast, get_args, get_origin, get_type_hints

from .annotation import PafAnnotation
from .comps import (
    Adduct,
    ChemicalFormula,
    ImmoniumIon,
    InternalFragment,
    IsotopeSpecification,
    MassError,
    NamedCompound,
    NeutralLoss,
    PeptideIon,
    PrecursorIon,
    ReferenceIon,
    SMILESCompound,
    UnknownIon,
)

_IONS = {
    cls.__name__: cls
    for cls in (
        PeptideIon,
        InternalFragment,
        PrecursorIon,
        ImmoniumIon,
        ReferenceIon,
        ChemicalFormula,
        NamedCompound,
        SMILESCompound,
        UnknownIon,
    )
}


def _coerce(value: object, hint: Any) -> Any:
    if get_origin(hint) is UnionType:
        for member in get_args(hint):
            try:
                return _coerce(value, member)
            except ValueError:
                pass
        raise ValueError("Value does not match an allowed type")
    if get_origin(hint) is Literal:
        if any(type(value) is type(choice) and value == choice for choice in get_args(hint)):
            return value
    elif isinstance(hint, type) and issubclass(hint, Enum):
        if isinstance(value, str):
            return hint(value)
    elif hint is float:
        if type(value) in (float, int):
            return value
    elif type(value) is hint:
        return value
    raise ValueError(f"Unexpected value type for {hint}")


def _load(cls: Any, data: object) -> Any:
    if not isinstance(data, Mapping):
        raise ValueError(f"{cls.__name__} must be an object")
    data = cast(Mapping[str, object], data)
    expected = {field.name for field in fields(cls)}
    if set(data) != expected:
        raise ValueError(f"{cls.__name__} requires exactly these fields: {sorted(expected)}")
    hints = get_type_hints(cls)
    values = {}
    for name in expected:
        try:
            values[name] = _coerce(data[name], hints[name])
        except ValueError as error:
            raise ValueError(f"Invalid {cls.__name__}.{name}: {error}") from error
    return cls(**values)


def to_dict(annotation: PafAnnotation) -> dict:
    """Return JSON-compatible component data with an explicit schema version."""
    data = asdict(annotation)
    ion = data.pop("ion_type")
    ion["type"] = type(annotation.ion_type).__name__
    return {"schema_version": 1, "ion": ion, **data}


def from_dict(data: Mapping[str, object]) -> PafAnnotation:
    """Reject unknown fields and invalid types before constructing components."""
    expected = {field.name for field in fields(PafAnnotation)} - {"ion_type"}
    expected.update(("schema_version", "ion"))
    if not isinstance(data, Mapping) or set(data) != expected:
        raise ValueError(f"Annotation requires exactly these fields: {sorted(expected)}")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise ValueError("Unsupported annotation schema_version")
    raw_ion = data["ion"]
    if not isinstance(raw_ion, Mapping):
        raise ValueError("Ion must be an object")
    ion = cast(Mapping[str, object], raw_ion)
    kind = ion.get("type")
    if not isinstance(kind, str) or kind not in _IONS:
        raise ValueError("Unknown or missing ion type")
    values: dict[str, Any] = {"ion_type": _load(_IONS[kind], {key: value for key, value in ion.items() if key != "type"})}
    for name, cls in (("neutral_losses", NeutralLoss), ("isotopes", IsotopeSpecification), ("adducts", Adduct)):
        components = data[name]
        if not isinstance(components, list | tuple):
            raise ValueError(f"{name} must be an array")
        values[name] = tuple(_load(cls, item) for item in components)
    values["mass_error"] = None if data["mass_error"] is None else _load(MassError, data["mass_error"])
    hints = get_type_hints(PafAnnotation)
    for name in ("analyte_reference", "is_auxiliary", "charge", "confidence", "resolved_sequence"):
        values[name] = _coerce(data[name], hints[name])
    annotation = PafAnnotation(**values)
    # Validate textual fields without requiring optional chemistry dependencies.
    PafAnnotation.parse(annotation.serialize())
    return annotation
