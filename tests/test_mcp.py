"""Exercise the optional server through the SDK and through an installed subprocess."""

import asyncio
import json
import os
import shutil
import socket
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp")
pytest.importorskip("peptacular")

import anyio
from jsonschema import validate
from mcp import Client
from mcp.client.stdio import StdioServerParameters

import paftacular as pft
from paftacular.mcp import create_server, handlers
from paftacular.mcp.guidance import EXAMPLES
from paftacular.mcp.models import MAX_REQUEST_BYTES, MAX_TEXT_BYTES


def call(name, request=None):
    async def run():
        async with Client(create_server()) as client:
            args = {"request": request} if request is not None else {}
            result = await client.call_tool(name, args)
            if result.structured_content is not None:
                tools = await client.list_tools()
                schema = next(tool.output_schema for tool in tools.tools if tool.name == name)
                validate(result.structured_content, schema)
                assert json.loads(result.content[0].text) == result.structured_content
            return result

    return asyncio.run(run())


def test_discovery_resources_and_prompts():
    async def run():
        async with Client(create_server()) as client:
            tools = await client.list_tools()
            assert {tool.name for tool in tools.tools} == {
                "get_capabilities",
                "parse_annotations",
                "calculate_ion",
                "calculate_ions",
                "resolve_annotation",
                "serialize_annotation",
                "build_annotation",
                "generate_fragments",
                "match_mz",
            }
            for tool in tools.tools:
                assert tool.input_schema and tool.output_schema
                assert tool.annotations.read_only_hint
                assert tool.annotations.idempotent_hint
                assert not tool.annotations.destructive_hint
                assert not tool.annotations.open_world_hint
            resources = await client.list_resources()
            assert len(resources.resources) == 4
            for resource in resources.resources:
                result = await client.read_resource(str(resource.uri))
                assert result.contents[0].text
            prompts = await client.list_prompts()
            assert {item.name for item in prompts.prompts} == {"analyze_fragment", "review_annotations"}
            analysis = await client.get_prompt("analyze_fragment", {"annotation": "y3^2", "analyte": "PEPTIDE"})
            assert "calculate_ion" in analysis.messages[0].content.text
            review = await client.get_prompt("review_annotations", {"records_json": '["y3", "bad"]'})
            assert "parse_annotations" in review.messages[0].content.text

    asyncio.run(run())


def test_capabilities():
    result = call("get_capabilities")
    assert not result.is_error
    data = result.structured_content["data"]
    assert data["dependencies"]["mcp"]
    assert data["integrations"]["peptacular"]
    assert data["limits"]["records"] == 100
    assert data["mass_types"] == ["monoisotopic"]


def test_parse_mixed_records():
    result = call("parse_annotations", {"records": ["y3{IDE},b2", "y3,garbage", "b2^0"]})
    assert not result.is_error
    rows = result.structured_content["data"]["records"]
    assert [item["canonical"] for item in rows[0]["annotations"]] == ["y3{IDE}", "b2"]
    assert rows[1]["annotations"] == []
    assert rows[1]["error"]["annotation_index"] == 1
    assert rows[1]["error"]["position"] == 3
    assert rows[2]["error"]["code"] == "parse_error"


def test_empty_parse_record():
    result = call("parse_annotations", {"records": ["", "y3"]})
    assert not result.is_error
    records = result.structured_content["data"]["records"]
    assert records[0]["annotations"] == []
    assert records[0]["error"] is None
    assert records[1]["annotations"][0]["canonical"] == "y3"


def test_bundled_tool_examples():
    name = None
    for line in EXAMPLES.splitlines():
        if line.endswith(":"):
            name = line[:-1]
        elif line.startswith("{"):
            assert name is not None
            assert not call(name, json.loads(line)["request"]).is_error


@pytest.mark.parametrize(
    "annotation,analyte",
    [
        ("y3^2", "PEPTIDE"),
        ("b3+i", "PEPTIDE"),
        ("m2:4", "PEPTIDE"),
        ("p^2", "[Acetyl]-PEPTIDE"),
        ("y3{IDE}", None),
        ("b2[M+2Na]^2", "PEPTIDE"),
        ("f{C2H4}+i", None),
        ("IM[Oxidation]", None),
        ("b3", "[Acetyl]-PEPTIDE"),
    ],
)
def test_calculation_matches_core(annotation, analyte):
    request = {"annotation": annotation, "properties": ["mass", "mz", "formula", "composition"]}
    expected = pft.parse_single(annotation)
    if analyte is not None:
        request["analyte"] = analyte
        expected = expected.resolve(analyte)
    result = call("calculate_ion", request)
    assert not result.is_error
    data = result.structured_content["data"]
    assert data["status"] == "success"
    assert data["mass_da"] == expected.mass()
    assert data["mz_th"] == expected.mz()
    assert data["formula"] == expected.formula()
    assert data["composition"] == expected.dict_composition()
    assert data["mass_basis"] == "charged_species"


def test_resolution_interchange_and_references():
    result = call("resolve_annotation", {"annotation": "2@y3^2", "analytes": [{"reference": 2, "sequence": "PEPTIDE"}]})
    assert not result.is_error
    data = result.structured_content["data"]
    assert data["sequence"] == "IDE"
    assert data["annotation"]["resolved_sequence"] == "IDE"
    result = call("serialize_annotation", {"annotation": data["annotation"]})
    assert result.structured_content["data"] == data
    assert data["notices"]
    calculated = call("calculate_ion", {"annotation": data["annotation"]})
    assert not calculated.is_error
    assert calculated.structured_content["data"]["context_source"] == "resolved"
    assert calculated.structured_content["data"]["mz_th"] == pft.parse_single("y3^2").resolve("PEPTIDE").mz()


def test_builder():
    result = call(
        "build_annotation",
        {"ion": "y3", "charge": 2, "neutral_losses": ["-H2O"], "isotopes": ["+i"], "mass_error": 1.2, "mass_error_unit": "ppm", "confidence": 0.9},
    )
    assert not result.is_error
    annotation = pft.parse_single(result.structured_content["data"]["canonical"])
    assert annotation.charge == 2
    assert annotation.neutral_losses[0] == pft.NeutralLoss.parse("-H2O")
    assert annotation.mass_error == pft.MassError(1.2, "ppm")
    assert annotation.confidence == 0.9
    assert call("build_annotation", {"ion": "y3^2"}).is_error


@pytest.mark.parametrize(
    "arguments",
    [
        {"annotation": "y3"},
        {"annotation": "p"},
        {"annotation": "m2:4"},
        {"annotation": "y30", "analyte": "PEPTIDE"},
        {"annotation": "m1:3", "analyte": "PEPTIDE"},
        {"annotation": "y3{PEP}", "analyte": "PEPTIDE"},
        {"annotation": "2@y3", "analytes": [{"reference": 1, "sequence": "PEPTIDE"}]},
        {"annotation": "d3", "analyte": "PEPTIDE"},
        {"annotation": "f{H2}", "analyte": "PEPTIDE"},
        {"annotation": "y3+iA", "analyte": "PEPTIDE"},
        {"annotation": "_{unknown}"},
        {"annotation": "?"},
        {"annotation": "r[definitely_missing]"},
    ],
)
def test_calculation_domain_errors(arguments):
    result = call("calculate_ion", arguments)
    assert result.is_error
    body = result.structured_content
    assert body["error"] is not None or body["data"]["status"] == "error"


def test_explicit_offsets():
    result = call("calculate_ion", {"annotation": "y3{IDE}", "mode": "offsets"})
    data = result.structured_content["data"]
    assert not result.is_error
    assert data["mass_da"] == pft.parse_single("y3").mass()
    assert data["mass_basis"] == "offsets_and_modifiers"
    assert data["context_source"] == "ignored"
    assert call("calculate_ion", {"annotation": "f{H2}", "mode": "offsets"}).is_error


def test_numeric_loss_partial_result():
    result = call("calculate_ion", {"annotation": "y3-18.01056", "analyte": "PEPTIDE", "properties": ["mass", "composition"]})
    assert not result.is_error
    data = result.structured_content["data"]
    assert data["status"] == "partial"
    assert data["mass_da"] > 0
    assert data["composition"] is None
    assert "composition" in data["property_errors"]


def test_batch_errors_preserve_order():
    result = call(
        "calculate_ions",
        {
            "requests": [
                {"annotation": "y3", "analyte": "PEPTIDE"},
                {"annotation": "bad"},
                {"annotation": "b2"},
                {"annotation": "y3+iA", "analyte": "PEPTIDE"},
            ]
        },
    )
    assert not result.is_error
    rows = result.structured_content["data"]["records"]
    assert [row["index"] for row in rows] == [0, 1, 2, 3]
    assert rows[0]["result"]["status"] == "success"
    assert rows[1]["error"]["code"] == "parse_error"
    assert rows[2]["error"]["code"] == "missing_context"
    assert rows[3]["result"]["status"] == "error"


def test_generate_fragments():
    result = call("generate_fragments", {"analyte": "[Acetyl]-PEPTIDE", "series": ["b", "y"], "charges": [1, 2], "positions": [2, 3]})
    assert not result.is_error
    rows = result.structured_content["data"]["records"]
    assert len(rows) == 8
    assert rows[0]["result"]["ion"]["canonical"] == "b2"
    assert rows[-1]["result"]["ion"]["canonical"] == "y3^2"
    for row in rows:
        data = row["result"]
        core = pft.parse_single(data["ion"]["canonical"]).resolve("[Acetyl]-PEPTIDE")
        assert data["mz_th"] == core.mz()
    assert call("generate_fragments", {"analyte": "A" * 60}).is_error
    assert call("generate_fragments", {"analyte": "PEPTIDE", "positions": [7]}).is_error


@pytest.mark.parametrize("unit,tolerance", [("ppm", 2.0), ("Th", 0.001)])
def test_matching(unit, tolerance):
    theoretical = pft.parse_single("y3^2").resolve("PEPTIDE").mz()
    result = call(
        "match_mz",
        {
            "observed_mz": theoretical + 0.0001,
            "tolerance": tolerance,
            "tolerance_unit": unit,
            "candidates": [{"annotation": "y3^2", "analyte": "PEPTIDE"}, {"annotation": "b2", "analyte": "PEPTIDE"}, {"annotation": "bad"}],
        },
    )
    assert not result.is_error
    data = result.structured_content["data"]
    assert data["matching_indices"] == [0]
    assert data["candidates"][0]["delta_th"] == pytest.approx(0.0001)
    assert data["candidates"][0]["delta_ppm"] == pytest.approx(0.0001 / theoretical * 1e6)
    assert data["candidates"][2]["error"]["code"] == "parse_error"


@pytest.mark.parametrize(
    "arguments",
    [
        {"annotation": "y3", "mode": "offsets", "analyte": "PEPTIDE"},
        {"annotation": "y3", "analyte": "PEPTIDE", "analytes": [{"reference": 1, "sequence": "PEPTIDE"}]},
        {"annotation": "y3", "analytes": [{"reference": 1, "sequence": "AA"}, {"reference": 1, "sequence": "BB"}]},
        {"annotation": "y3", "analytes": [{"reference": True, "sequence": "AA"}]},
        {"annotation": "y3", "analytes": [{"reference": "1", "sequence": "AA"}]},
        {"annotation": "y3", "properties": ["unsupported"]},
        {"annotation": "y3", "ignored_field": True},
        {"annotation": "a" * (MAX_TEXT_BYTES + 1)},
    ],
)
def test_invalid_schemas(arguments):
    assert call("calculate_ion", arguments).is_error


def test_application_limits():
    assert call("parse_annotations", {"records": ["y3"] * 101}).is_error
    result = call("parse_annotations", {"records": ["a" * MAX_TEXT_BYTES] * 20})
    assert result.structured_content["error"]["code"] == "request_limit"
    result = call("parse_annotations", {"records": ["é" * MAX_TEXT_BYTES]})
    assert result.structured_content["error"]["code"] == "request_limit"
    result = call("serialize_annotation", {"annotation": {"nested": "a" * MAX_REQUEST_BYTES}})
    assert result.is_error
    # A small input can expand to many structured annotations.
    result = call("parse_annotations", {"records": [",".join(["b2"] * 100)] * 20})
    assert result.structured_content["error"]["code"] == "result_limit"
    result = call("parse_annotations", {"records": ["_{" + "a" * 2000 + "}"] * 100})
    assert result.structured_content["error"]["code"] == "result_limit"
    assert "512 KiB" in result.structured_content["error"]["message"]


def test_smiles_isotope_calculation():
    pytest.importorskip("pysmiles")
    result = call("calculate_ion", {"annotation": "s{[13CH4]}", "properties": ["mass", "composition"]})
    assert not result.is_error
    assert result.structured_content["data"]["composition"]["13C"] == 1


def test_missing_smiles_dependency(monkeypatch):
    monkeypatch.setitem(sys.modules, "pysmiles", None)
    result = call("calculate_ion", {"annotation": "s{C}"})
    assert result.is_error
    assert result.structured_content["data"]["property_errors"]["mass"]["code"] == "missing_dependency"


def test_offline_and_stdout(monkeypatch, capsys):
    def no_network(*args, **kwargs):
        raise AssertionError("Network access attempted")

    monkeypatch.setattr(socket, "create_connection", no_network)
    original = handlers.calculate_ion

    def noisy(request):
        print("dependency diagnostic")
        return original(request)

    monkeypatch.setattr(handlers, "calculate_ion", noisy)
    result = call("calculate_ion", {"annotation": "y3^2", "analyte": "PEPTIDE"})
    assert not result.is_error
    captured = capsys.readouterr()
    assert "dependency diagnostic" not in captured.out
    assert "dependency diagnostic" in captured.err


def test_unexpected_failure_is_sanitized(monkeypatch):
    def fail():
        raise RuntimeError("private diagnostic")

    monkeypatch.setattr(handlers, "capabilities", fail)
    result = call("get_capabilities")
    assert result.is_error
    assert result.structured_content["error"]["code"] == "internal_error"
    assert "private diagnostic" not in result.content[0].text


@pytest.mark.parametrize("launch", ["module", "command"])
@pytest.mark.parametrize("mode", ["auto", "legacy"])
def test_stdio_connection(launch, mode, tmp_path):
    async def run():
        command = sys.executable
        args = ["-I", "-m", "paftacular.mcp"]
        if launch == "command":
            executable = Path(sys.executable).parent / ("paftacular-mcp.exe" if os.name == "nt" else "paftacular-mcp")
            command = str(executable) if executable.exists() else shutil.which("paftacular-mcp")
            assert command is not None
            args = []
        params = StdioServerParameters(command=command, args=args, cwd=tmp_path)
        with anyio.fail_after(25):
            async with Client(params, read_timeout_seconds=10, mode=mode) as client:
                assert len((await client.list_tools()).tools) == 9
                bad = await client.call_tool("calculate_ion", {"request": {"annotation": "y3"}})
                assert bad.is_error
                good = await client.call_tool("calculate_ion", {"request": {"annotation": "y3^2", "analyte": "PEPTIDE"}})
                assert not good.is_error
                assert good.structured_content["data"]["ion"]["sequence"] == "IDE"
                assert len((await client.list_resources()).resources) == 4
                assert len((await client.list_prompts()).prompts) == 2

    asyncio.run(run())
