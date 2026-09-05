# Optional MCP server plan

Status: implemented and prepared for 1.3.0. See docs/mcp.rst for the current contract.
Prepared: 2026-09-04, against paftacular 1.2.0.
Release target: 1.3.0.

Implementation expanded to nine tools, four resources, and two prompts.
Additional tools cover annotation construction, standalone resolution,
terminal-fragment generation, and candidate m/z matching. Parsing also has
a 1000-annotation limit. Aggregate limits count encoded JSON bytes.
The original design below records the starting plan.

## Outcome

Let an AI application use paftacular through Model Context Protocol to parse
mzPAF, calculate fragment masses with explicit peptide context, and exchange
structured annotations. Reuse the public Python API for all chemistry.

The first version runs as a local subprocess over stdio. An MCP-capable host
launches it and makes tool calls on behalf of the model. Installing the extra
provides the server, while connecting it remains a host configuration step.
Remote-only hosts will need a future HTTP deployment option.

## Packaging and startup

The proposed primary installation is:

```bash
pip install 'paftacular[mcp]'
paftacular-mcp
```

These commands become available when the feature ships.

- Add an `mcp` extra containing the official Python MCP SDK and the existing
  peptacular dependency. Including peptide support makes the main AI use case
  work with one install.
- Use the SDK v2 `MCPServer` interface. Start with `mcp>=2.1.1,<3`, then verify
  that range in the implementation matrix. Version 2.1.1 is the current stable
  release verified during planning. Avoid the separate `fastmcp` distribution.
- Support `paftacular[mcp,smiles]` for chemical-structure calculations.
- Extend `all` to include MCP, documenting its additional dependencies.
- Keep the base install and the existing `peptacular` and `smiles` extras free
  of MCP imports and dependencies.
- Add `[project.scripts]` with `paftacular-mcp` pointing to a small CLI module.
  Python extras cannot conditionally install entry points, so invoking the
  command without the extra should exit with a clear install instruction.
- Also support `python -m paftacular.mcp`. Load the SDK only when creating or
  running the server. The CLI uses argparse and does not require `mcp[cli]`.
- Support `--help` and `--version`. Normal execution speaks MCP on stdout,
  with diagnostics on stderr and no startup banner on stdout.

The SDK's current stable v2 line and optional CLI dependency are documented in
the [official installation guide](https://py.sdk.modelcontextprotocol.io/get-started/installation/).
The planned minimum comes from the [v2.1.1 release](https://github.com/modelcontextprotocol/python-sdk/releases/tag/v2.1.1).

## Initial tools

Use five tools with explicit input and output schemas. Register the same tools
regardless of whether SMILES support is installed, and report availability in
capabilities and actionable calculation errors.

| Tool | Inputs | Result |
| --- | --- | --- |
| `get_capabilities` | None | Package and dependency versions, available integrations, supported calculations, conventions, and limits |
| `parse_annotations` | A list of mzPAF text records | Ordered records with canonical text, versioned annotation data, or structured parse errors |
| `calculate_ion` | One annotation, optional analyte context, requested properties, calculation mode | Resolved sequence, requested numerical or composition results, and explicit scientific context |
| `calculate_ions` | A list of calculation requests using the same schema | Ordered results with independent failures |
| `serialize_annotation` | One versioned annotation dictionary | Validated canonical mzPAF text, structured data, and a context-preservation notice when relevant |

For calculations, accept either one `analyte` string or an `analytes` list of
objects with integer `reference` and string `sequence` fields. Reject supplying
both or duplicate references. Convert the latter to the existing Python
mapping internally. This avoids JSON object keys changing integer references
into strings. Preserve the existing rule that an omitted annotation reference
selects analyte 1.

The `properties` array selects `mass`, `mz`, `formula`, or `composition`.
Default to `mass` and `mz`. Return mass in Da with `mass_basis` set to
`charged_species`, m/z in Th, and charge explicitly. The first MCP release
exposes monoisotopic calculations. Average-mass support needs separate
scientific validation before it becomes an MCP option.

An illustrative tool request is:

```json
{
  "annotation": "y3^2",
  "analyte": "PEPTIDE",
  "properties": ["mass", "mz", "composition"],
  "mode": "complete"
}
```

It must select `IDE` through `resolve()` and return the same results as the
public library. The adapter must not implement its own sequence slicing or
mass tables.

## Scientific behavior

- Default to `mode="complete"`. A peptide, internal, or precursor annotation
  without embedded or supplied sequence context returns `missing_context`.
  It must not return an offset under a complete-mass label.
- Allow `mode="offsets"` explicitly for callers requesting the existing
  context-free behavior. Mark the result as offsets, omit any suggestion that
  it represents a complete fragment, and reject supplied analyte context in
  this mode rather than silently ignoring it.
- Reject analyte context for ion types to which it does not apply.
- Report whether sequence context came from the annotation or the analyte.
  Retain resolved data in the structured response because mzPAF text does not
  preserve that context.
- Return composition using element and isotope names as keys. Preserve signed
  counts and apply the existing electron correction only through library APIs.
- Treat formula-ion charge conventions and internal-cleavage conventions as
  part of the documented tool contract. Do not infer a physical cleavage type
  from an ambiguous specification correction.
- Calculate requested properties independently. A numeric neutral loss can
  permit a mass result while preventing a complete composition. Return a
  per-property error with record status `partial`, retaining valid results.
- Report unsupported side-chain resolution, unresolved named or unknown ions,
  average isotopomers, and unknown chemical references explicitly. Do not
  substitute guesses or zero values.

## Protocol and output contract

Use typed adapter models and an MCP envelope with its own
`response_schema_version: 1`. Keep the library's existing dictionary schema
version independent. Include package version, record status, results, and
structured errors. Error fields include a stable code, a useful message, and
the existing zero-based parse position and annotation index when applicable.

Return structured content plus its JSON text representation for clients that
consume text tool results. Publish output schemas and mark tool annotations
as read-only, non-destructive, idempotent, and closed-world. Domain failures
use MCP tool error results. A completed batch with mixed outcomes returns an
ordinary batch result with record-level errors, while an invalid whole request
returns a tool error. Unexpected exceptions go to stderr and become a concise
internal-error result, with no fabricated scientific output.

This follows the protocol's
[tool result and error conventions](https://modelcontextprotocol.io/specification/2025-11-25/server/tools)
and [stdio transport requirements](https://modelcontextprotocol.io/specification/2025-11-25/basic/transports).

Set initial application limits of 100 records per call, 16 KiB per input text,
256 KiB per decoded request, and 512 KiB per structured result. Count nested
analytes and dictionaries toward the totals. Reject oversized work with a
specific limit error and guidance to split the request, without silent
truncation. Test worst-case accepted inputs before fixing these defaults.

Keep calculations stateless. The first server accepts values directly and
does not expose file access, arbitrary execution, or network lookup tools.
Ensure dependency logging and any incidental output cannot corrupt stdout.
Verify ordinary calls work offline after installation.

## Implementation sequence

1. **Package boundary and server skeleton.** Add the extra, entry point,
   dependency constraints, CLI, lazy imports, and server factory. Verify core
   imports without MCP and a stdio initialize/list-tools exchange.
2. **Typed adapters.** Add request and response models, all five tools, shared
   calculation handling, capability reporting, limits, and error translation.
   Reuse parse, resolve, to_dict, from_dict, mass, mz, and composition APIs.
3. **Scientific and protocol validation.** Exercise full and offset modes,
   analyte references, modification handling, partial results, error recovery,
   and multiple requests in one server session.
4. **Installation and documentation.** Add isolated wheel tests and client
   setup examples. Verify a real MCP host can discover and call the tools.
5. **Release review.** Review schemas and descriptions as public APIs, run the
   existing release checks plus MCP checks, and prepare a 1.3.0 changelog entry.
   Implementation approval does not itself publish a release.

Suggested new files are `src/paftacular/mcp/__init__.py`, `__main__.py`,
`cli.py`, `server.py`, `models.py`, and `handlers.py`, plus focused MCP tests
and `docs/mcp.rst`. Update pyproject.toml, uv.lock, the minimum-dependency
constraints, installation CI, README, documentation navigation, and HISTORY.
Preserve the user's existing unrelated Zenodo change.

## Acceptance checks

- A fresh `paftacular[mcp]` wheel installation starts the command and calculates
  `y3^2` from `PEPTIDE`, returning `IDE` and values matching the core API.
- A fresh base installation imports and parses normally without MCP installed.
  The unavailable server command exits with an actionable message.
- In-memory SDK tests validate discovery, schemas, tool annotations, structured
  output, and domain errors. A real stdio subprocess test verifies framing,
  clean stdout, continued operation after invalid input, and shutdown.
- Regression fixtures cover embedded and supplied ProForma, multiple analytes,
  internal and precursor ions, isotopes, explicit adducts, formula ions,
  numeric-loss partial results, missing context, and unsupported calculations.
- Batch output preserves input order and failures. Oversized inputs and output
  limits fail predictably. Serialization preserves the existing dictionary
  schema and reports context loss from text output.
- CI covers Python 3.12, 3.13, and 3.14, minimum and latest supported MCP SDK
  versions, `mcp`, `mcp,smiles`, and `all` wheels. Retain existing core and
  chemistry-extra matrices. Include a Windows stdio smoke test.
- Document tested host setup using executable and argument arrays. Include an
  isolated `uvx --from 'paftacular[mcp]==1.3.0' paftacular-mcp` example once
  that version exists. No provider API key is required by this local server.

## Later work

Consider Streamable HTTP hosting for remote clients, authentication for that
deployment, registry publication, independently validated average-mass
calculations, and larger spectral-library workflows after the local interface
has demonstrated demand. Keep these out of the initial implementation.
