"""Guard the MCP tool vocabulary against drift from the library and the shared tacular names."""

import asyncio
import re
from typing import get_args

import pytest

pytest.importorskip("mcp")
pytest.importorskip("peptacular")

from mcp.server.mcpserver.exceptions import ToolError
from tacular.types import ToleranceUnit

import paftacular as pft
from paftacular.mcp.server import create_server

# Old or non-shared spellings that must not reappear in any input or output schema.
BANNED = re.compile(
    r"^(unit|tolerance_type|.*_tolerance_type|retention_time.*|inverse_reduced.*|ion_mobility_.*|target_mz|"
    r"scan_start_time|ce|tic|TIC|time|one_over_k0.*|mz_begin|mz_end|window_group|monoisotopic_mz|ion_series|resolvable_series)$"
)

# Property names ending in "unit" that are not a Da/ppm tolerance switch.
# mass_error_unit is the mzPAF mass-error notation field (library MassError.unit), not a
# tolerance: the round-2 vocabulary spec keeps it as is. Its values must follow the library.
ALLOWED_UNIT_NAMES = {"mass_error_unit"}

# Keys of calculation records that describe the MCP calculation rather than a PafAnnotation attribute.
MCP_ONLY_CALCULATION_KEYS = {
    "status",
    "ion",
    "mode",
    "context_source",
    "mass_type",
    "mass_basis",
    "mass_da",
    "mz_th",
    "composition",
    "property_errors",
}
MCP_ONLY_VIEW_KEYS = {"canonical", "annotation", "notices"}

# Minimal valid arguments for every tool.
MINIMAL = {
    "get_capabilities": {},
    "parse_annotations": {"request": {"records": ["y3"]}},
    "resolve_annotation": {"request": {"annotation": "y3", "analyte": "PEPTIDE"}},
    "serialize_annotation": {"request": {"annotation": pft.parse("y3").to_dict()}},
    "build_annotation": {"request": {"ion": "y3"}},
    "calculate_ion": {"request": {"annotation": "y3", "analyte": "PEPTIDE"}},
    "calculate_ions": {"request": {"requests": [{"annotation": "y3", "analyte": "PEPTIDE"}]}},
    "generate_fragments": {"request": {"analyte": "PEPTIDE", "positions": [2]}},
    "match_mz": {"request": {"observed_mz": 100.0, "candidates": [{"annotation": "y3", "analyte": "PEPTIDE"}]}},
}


@pytest.fixture(scope="module")
def server():
    return create_server()


@pytest.fixture(scope="module")
def tools(server):
    return asyncio.run(server.list_tools())


def schema_props(schema):
    """Yield (path, name, subschema) for every property, following $ref into $defs once each."""
    defs = schema.get("$defs", {})
    seen = set()

    def walk(node, path):
        if not isinstance(node, dict):
            return
        ref = node.get("$ref")
        if ref:
            name = ref.rsplit("/", 1)[-1]
            if name not in seen:
                seen.add(name)
                yield from walk(defs.get(name, {}), f"{path}<{name}>")
        for key, value in (node.get("properties") or {}).items():
            yield f"{path}.{key}", key, value
            yield from walk(value, f"{path}.{key}")
        for key in ("items", "anyOf", "oneOf", "allOf", "additionalProperties"):
            value = node.get(key)
            for child in value if isinstance(value, list) else [value]:
                yield from walk(child, path)

    yield from walk(schema, "")


def enum_values(sub):
    values = set()
    for node in [sub, *sub.get("anyOf", [])]:
        values |= set(node.get("enum", []))
        if "const" in node:
            values.add(node["const"])
    return values


def all_props(tools):
    for tool in tools:
        for kind, schema in (("input", tool.input_schema), ("output", tool.output_schema or {})):
            for path, name, sub in schema_props(schema):
                yield f"{tool.name} {kind}{path}", name, sub


def test_tools_covered(tools):
    assert {tool.name for tool in tools} == set(MINIMAL)


def test_no_banned_names(tools):
    offenders = [path for path, name, _ in all_props(tools) if BANNED.match(name)]
    assert not offenders


def test_capability_names(tools):
    # Lists of ion kinds are plural ion_types, like peptacular. The singular is a per-record field there.
    capabilities = next(tool for tool in tools if tool.name == "get_capabilities")
    names = {name for path, name, _ in schema_props(capabilities.output_schema) if "<Capabilities>" in path}
    assert {"ion_types", "resolvable_ion_types"} <= names
    assert not names & {"ion_type", "ion_series", "resolvable_series"}


def test_tolerance_switches(tools):
    tolerance_values = set(get_args(ToleranceUnit))
    for path, name, sub in all_props(tools):
        if name.endswith("unit") and name not in ALLOWED_UNIT_NAMES:
            assert re.fullmatch(r"(\w+_)?tolerance_unit", name), path
            assert enum_values(sub) == tolerance_values, path
        if name in ALLOWED_UNIT_NAMES:
            # Not a tolerance: follows the library's mzPAF mass-error unit, whatever its values.
            library_units = set(get_args(pft.MassError.__dataclass_fields__["unit"].type))
            assert enum_values(sub) == library_units, path


def test_charges_are_signed(tools):
    # Library 2.0 accepts negative charges, so no charge input may have a positive-only floor.
    found = 0
    for path, name, sub in all_props(tools):
        if name in ("charge", "charges") and " input" in path:
            found += 1
            item = sub.get("items", sub)
            assert item.get("type") == "integer", path
            assert item.get("minimum", -1) < 0, path
    assert found == 2


def test_unknown_argument_rejected(server):
    async def run():
        for name, arguments in MINIMAL.items():
            result = await server.call_tool(name, arguments)
            assert not getattr(result, "is_error", False), name
            with pytest.raises(ToolError, match="Extra inputs are not permitted"):
                await server.call_tool(name, {**arguments, "__bogus__": 1})
            if "request" in arguments:
                with pytest.raises(ToolError, match="Extra inputs are not permitted"):
                    await server.call_tool(name, {"request": {**arguments["request"], "__bogus__": 1}})

    asyncio.run(run())


def test_record_keys_match_library(server):
    result = asyncio.run(
        server.call_tool("calculate_ion", {"request": {"annotation": "y3^-2", "analyte": "PEPTIDE", "properties": ["mass", "mz", "formula", "composition"]}})
    )
    data = result.structured_content["data"]
    library = pft.parse("y3^-2").resolve("PEPTIDE")
    attributes = {name for name in dir(library) if not name.startswith("_")}
    assert set(data) - attributes <= MCP_ONLY_CALCULATION_KEYS
    assert set(data["ion"]) - attributes <= MCP_ONLY_VIEW_KEYS
    assert set(data["ion"]["annotation"]) == set(library.to_dict())
    assert pft.PafAnnotation.from_dict(data["ion"]["annotation"]) == library
    assert data["charge"] == library.charge
    assert data["formula"] == library.formula()
