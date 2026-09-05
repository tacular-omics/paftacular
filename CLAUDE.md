# paftacular agent guide

## Project and architecture

paftacular parses and serializes HUPO-PSI mzPAF peak annotations and calculates
ion masses and elemental compositions. Python 3.12 or newer is required.
`tacular` supplies chemistry lookup data. `peptacular` and `pysmiles` are
optional integrations, installed from public package indexes without local
source overrides.

- `annotation.py`: frozen `PafAnnotation`, factories, calculation and serialization
- `parser.py`: single and multiple parsing, lazy batches, `ParseResult`
- `syntax.py`: delimiter-aware annotation boundaries, including nested ProForma
- `errors.py`: `PafParseError`, a `ValueError` subclass with zero-based locations
- `resolution.py`: explicit full-analyte context for peptide, internal, precursor ions
- `serialization.py`: versioned component interchange, separate from legacy `as_dict()`
- `constants.py`: enums, grammar patterns, specification internal-cleavage table
- `comps/ions.py`: ion components and physical backbone corrections
- `comps/modifiers.py`: neutral losses, isotopes, adducts, mass errors
- `comps/base.py`: composition and serialization interfaces, pickle reconstruction
- `util.py` and `comps/util.py`: decimal formatting and formula conversion
- `mcp/`: optional local MCP server, typed contracts, public-API adapters, resources, and prompts

## Commands

```bash
just install       # uv sync --all-extras
just install-prod  # uv sync --no-dev --frozen
just lint          # ruff on source, tests, benchmarks
just format        # format source, tests, benchmarks
just check         # lint, format check, ty, tests
just test-cov      # branch coverage
just docs          # Sphinx HTML with warnings as errors
just docs-test     # executable Sphinx examples
just benchmark     # timing and retained-memory measurements
just build         # source distribution and wheel
just pre-release   # format, checks, documentation, package build
```

## Scientific conventions

- `mass()` returns charged species mass. `comp()` counts nuclei. Compare summed
  element masses after subtracting charge times electron mass. Upstream ion
  offsets are rounded, so independent checks need an absolute tolerance around
  one microdalton, with relative tolerance disabled.
- Without embedded or resolved context, peptide and precursor calculations
  return only offsets and modifiers. Preserve this existing behavior. Use
  `resolve()` to select the complete fragment sequence from an analyte.
- Embedded sequences contribute residue mass and composition using peptacular
  `ion_type="n"`. Default precursor composition adds an incorrect extra water.
- Formula ions already specify every atom in the charged species. Their adducts
  label charge carriers and add no atoms. Other ion types add adduct atoms and
  subtract electron mass. Implicit protonation and explicit `[M+H]` must agree.
- A generic isotope uses the carbon-13 minus carbon-12 shift. Add contributions
  to counters rather than overwriting them when isotope keys coincide. At the
  annotation level, consume ordinary monoisotopic atoms when available. Keep
  genuine deficits when sequence context is absent. Average isotopomers raise
  on mass and composition calculations.
- Counter addition and unary plus discard negative entries. Accumulate with
  `update()` and preserve genuine negative counts from modifications.
- Immonium ions use tacular's `"i"` lookup key, not the internal `"by"` key.
- Explicit internal cleavage fields must be supplied together, in a/b/c and
  x/y/z combinations. Their physical composition comes from their own tacular
  key. Serialize it as signed elemental gains/losses, preserving mass.
- The mzPAF 1.0.1 section 4.4.4 correction table disagrees with physical tacular
  definitions for several combinations. Public specification tables and
  `make_internal()` retain the specification convention. Conversion and explicit
  cleavage fields preserve the source's physical composition. Do not silently
  equate these conventions. Their distinction is documented in `docs/usage.rst`.

## Parser and interchange rules

- Locate commas outside labels and sequences before matching whole annotations.
  Never split blindly on commas or let a sequence consume another annotation.
- Preserve possessive quantifiers in `_ATOM_TOKEN`. Removing them can introduce
  catastrophic backtracking on invalid long element runs.
- Optional leading plus signs in mass errors and spaces/commas in named labels
  are intentional compatibility choices. Keep these accepted.
- Use `format_number()` for grammar-compatible decimal output without truncating
  float precision. Scientific notation is not part of the supported grammar.
- Component parsers consume the entire input. Zero charge, nonintegral positions,
  and reversed ranges must fail. Resolution checks bounds against the analyte.
- `to_dict()` and `from_dict()` use schema version 1 and reject unknown fields.
  Keep `as_dict()` compatible as a display format. Resolved context is stored in
  structured export and is deliberately absent from mzPAF text serialization.
- Optional dependency paths must fail with actionable installation instructions.
  Importing and parsing the core package must work without either extra.

## Caching and tests

Selected frozen components cache constructor instances, bounded by
`MAX_CACHE_SIZE`. Validate arguments before caching to prevent invalid calls
from mutating a previously cached object. Shared `Serializable.__reduce__`
reconstructs dataclasses through constructors for pickle/copy compatibility.
Measure changes with `benchmarks/parse.py` before expanding caching.

Tests cover examples, scientific invariants, component round trips, error
locations, resolution, interchange, and optional installations. CI runs locked
Python 3.12/3.13/3.14 environments plus built-wheel installations with minimum
and newer dependencies for each extra. The full suite requires all extras.
The dedicated installation and public-API suites run with optional skips.

## Optional MCP integration

`paftacular[mcp]` includes the official SDK v2 and peptacular. `all` includes MCP.
Keep SDK and Pydantic imports inside `mcp/`, loaded only when creating the server.
The CLI and importing `paftacular.mcp` must work without SDK imports. The server
runs on stdio, with protocol output on stdout and diagnostics on stderr.

MCP exposes tools, static reference resources, and prompts. Reuse core chemistry
APIs. Complete peptide calculations require sequence context. Offset mode must
be explicit. Monoisotopic results use charged-species mass in Da and m/z in Th.
Preserve partial property errors and structured context across calls. Request
and result contracts are versioned independently from core interchange.

Tool wrappers are async with bounded synchronous handlers and no inner awaits.
This serializes chemistry calls within a server event loop and avoids racing
the core component caches. Test with the official MCP Client in memory and with
real stdio subprocesses from outside the checkout.

## Release preparation

Version is in `src/paftacular/__init__.py`, sourced by hatchling. Maintain
`HISTORY.md` under Unreleased until assigning a release version and date.
Run `just pre-release` and installation checks before publication. The publish
workflow requires checks and a release tag matching the package version.

Remote: https://github.com/tacular-omics/paftacular
Specification: https://www.psidev.info/mzpaf
