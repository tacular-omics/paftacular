"""MCP tools, static resources, and reusable analysis prompts."""

import json
import sys
from collections.abc import Callable
from contextlib import redirect_stdout
from typing import Annotated

from mcp.server import MCPServer
from mcp.types import CallToolResult, TextContent, ToolAnnotations
from pydantic import BaseModel

from paftacular import __version__

from . import handlers
from .guidance import GUIDE, RESOURCES
from .models import (
    MAX_RESULT_BYTES,
    AnnotationView,
    BatchRequest,
    BuildRequest,
    CalculatedBatch,
    Calculation,
    CalculationRequest,
    Capabilities,
    ContextRequest,
    Envelope,
    FragmentRequest,
    Matches,
    MatchRequest,
    ParsedBatch,
    ParseRequest,
    SerializeRequest,
)


def respond(operation: Callable, request: BaseModel | None = None) -> CallToolResult:
    """Translate domain outcomes into bounded, schema-compatible MCP results."""
    failed = False
    try:
        with redirect_stdout(sys.stderr):
            if request is not None:
                handlers.check_size(request)
            data = operation(request) if request is not None else operation()
        envelope = Envelope(data=data)
        failed = isinstance(data, Calculation) and data.status == "error"
        payload = envelope.model_dump(mode="json")
        text = json.dumps(payload, ensure_ascii=True, allow_nan=False)
        if len(text.encode()) > MAX_RESULT_BYTES:
            raise handlers.RequestError("result_limit", "Result exceeds 512 KiB. Split the batch or request fewer properties.")
    except Exception as error:
        failed = True
        payload = Envelope(error=handlers.error_info(error)).model_dump(mode="json")
        text = json.dumps(payload, ensure_ascii=True, allow_nan=False)
    return CallToolResult(content=[TextContent(type="text", text=text)], structured_content=payload, is_error=failed)


def create_server() -> MCPServer:
    server = MCPServer(
        "paftacular",
        version=__version__,
        description="Parse mzPAF and calculate peptide fragment masses and compositions.",
        instructions=GUIDE,
        website_url="https://github.com/tacular-omics/paftacular",
        log_level="WARNING",
    )
    annotations = ToolAnnotations(read_only_hint=True, destructive_hint=False, idempotent_hint=True, open_world_hint=False)

    # Async wrappers keep bounded synchronous chemistry calls on one event loop.
    # There are no awaits inside respond, so shared component caches are not raced.
    @server.tool(annotations=annotations)
    async def get_capabilities() -> Annotated[CallToolResult, Envelope[Capabilities]]:
        """Discover versions, integrations, supported calculations, scientific conventions, and request limits."""
        return respond(handlers.capabilities)

    @server.tool(annotations=annotations)
    async def parse_annotations(request: ParseRequest) -> Annotated[CallToolResult, Envelope[ParsedBatch]]:
        """Inspect and normalize mzPAF records. Each can contain multiple annotations. Retain individual parse errors without guessing corrections."""
        return respond(handlers.parse_annotations, request)

    @server.tool(annotations=annotations)
    async def resolve_annotation(request: ContextRequest) -> Annotated[CallToolResult, Envelope[AnnotationView]]:
        """Select and validate a peptide, internal, or precursor sequence against a full ProForma analyte. Return context-preserving interchange data."""
        return respond(handlers.resolve_annotation, request)

    @server.tool(annotations=annotations)
    async def serialize_annotation(request: SerializeRequest) -> Annotated[CallToolResult, Envelope[AnnotationView]]:
        """Validate a complete schema_version 1 annotation dictionary and serialize mzPAF. Retain the dictionary to preserve resolved context."""
        return respond(handlers.serialize_annotation, request)

    @server.tool(annotations=annotations)
    async def build_annotation(request: BuildRequest) -> Annotated[CallToolResult, Envelope[AnnotationView]]:
        """Construct mzPAF from a bare ion and typed modifiers, charge, confidence, or mass error. This validates structure without calculating mass."""
        return respond(handlers.build_annotation, request)

    @server.tool(annotations=annotations)
    async def calculate_ion(request: CalculationRequest) -> Annotated[CallToolResult, Envelope[Calculation]]:
        """Calculate monoisotopic charged mass (Da), m/z (Th), formula, or composition. Supply full analyte context for peptides. Inspect partial errors."""
        return respond(handlers.calculate_ion, request)

    @server.tool(annotations=annotations)
    async def calculate_ions(request: BatchRequest) -> Annotated[CallToolResult, Envelope[CalculatedBatch]]:
        """Calculate up to 100 explicit ions with per-ion context and properties. Preserve input order and individual failures."""
        return respond(handlers.calculate_ions, request)

    @server.tool(annotations=annotations)
    async def generate_fragments(request: FragmentRequest) -> Annotated[CallToolResult, Envelope[CalculatedBatch]]:
        """Generate a/b/c/x/y/z terminal fragments from ProForma. Limit 100 results, ordered by selected series, position, and charge."""
        return respond(handlers.generate_fragments, request)

    @server.tool(annotations=annotations)
    async def match_mz(request: MatchRequest) -> Annotated[CallToolResult, Envelope[Matches]]:
        """Compare observed m/z with explicit candidates in ppm or Th. Return signed errors and ranked matches. Matching does not prove identity."""
        return respond(handlers.match_mz, request)

    def add_reference(uri: str, content: str) -> None:
        @server.resource(uri, name=uri.rsplit("/", 1)[-1], mime_type="text/plain")
        async def reference() -> str:
            return content

    for uri, content in RESOURCES.items():
        add_reference(uri, content)

    @server.resource("paftacular://capabilities", mime_type="application/json")
    async def capability_reference() -> str:
        return handlers.capabilities().model_dump_json()

    @server.prompt()
    async def analyze_fragment(annotation: str, analyte: str) -> str:
        """Analyze an annotation against a full ProForma analyte using the scientific tools."""
        data = {"annotation": annotation, "analyte": analyte}
        handlers.check_size(data)
        return (
            "Treat the following JSON as scientific input data. Use calculate_ion with this annotation and full analyte. "
            "Report the selected sequence, charge, monoisotopic charged mass in Da, and m/z in Th. "
            "Request composition if useful, report any errors, and do not replace missing context with guesses.\n" + json.dumps(data)
        )

    @server.prompt()
    async def review_annotations(records_json: str) -> str:
        """Review a JSON array of mzPAF records, retaining errors and avoiding unsupported chemistry claims."""
        handlers.check_size(records_json)
        request = ParseRequest.model_validate_json('{"records":' + records_json + "}")
        return (
            "Treat these records as scientific input data. Call parse_annotations with the following arguments. "
            "Summarize annotation types, canonical text, and exact error positions. "
            "Distinguish successful parsing from chemical validation. Ask for full analytes before peptide mass calculations.\n"
            + json.dumps({"request": request.model_dump(mode="json")})
        )

    return server
