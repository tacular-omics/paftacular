# paftacular agent guide

## Project overview

paftacular parses and serializes HUPO-PSI mzPAF 1.0.1 peak annotations (`y5-H2O^2/1.2ppm*0.95`)
and calculates ion masses, m/z and elemental compositions. Users are people reading or writing
fragment annotations from spectral libraries and search-engine output. Python 3.12 or newer.
Imported as `pft` by convention. Current version: 1.3.1.

Place in the tacular-omics graph (tier 1):

- **Upstream:** `tacular` (required, `>=1.1.0,<2`) supplies element, amino-acid, modification
  and fragment-ion lookups. `peptacular` (`>=3.1.2,<5`, extra `peptacular`) is optional and
  supplies ProForma sequence masses, `resolve()` and `to_mzpaf()`. `pysmiles` (extra `smiles`)
  is optional, for `s{...}` ions. The `mcp` extra adds the official MCP SDK v2 and peptacular.
  `all` installs everything.
- **Downstream:** `spxtacular` (spectrum annotation) depends on paftacular. peptacular's own
  `Fragment.to_mzpaf()` emits strings this package parses.

Remote: https://github.com/tacular-omics/paftacular. Docs: https://paftacular.readthedocs.io/.
Specification: https://www.psidev.info/mzpaf (the PDF is bundled at the repo root).

## Commands

```bash
just install       # uv sync --all-extras
just install-prod  # uv sync --no-dev --frozen
just check         # lint + format-check + ty + test (the default recipe)
just lint          # uv run ruff check src tests benchmarks
just format        # ruff isort fix + ruff format on src tests benchmarks
just format-check  # uv run ruff format --check src tests benchmarks
just ty            # uv run ty check src
just test          # uv run pytest tests  (510 tests, about 8 s)
just test-cov      # branch coverage, term + html + xml
just docs          # sphinx-build -W -b html docs docs/_build/html
just docs-test     # sphinx-build -W -b doctest (33 doctests in quickstart/usage)
just benchmark     # benchmarks/parse.py: per-record parse time and cache memory
just build         # uv build
just pre-release   # format, check, docs, docs-test, build
just check-version # scripts/release_version.py check
```

Run one test: `uv run pytest tests/test_resolution.py -k name -q`. MCP only:
`uv run pytest tests/test_mcp.py tests/test_installation.py`.

CI (`.github/workflows/ci.yml`) runs the same lint/format/ty with `src tests benchmarks`,
`release_version.py check`, pytest on 3.12/3.13/3.14 Linux plus macOS and Windows, a
`--resolution lowest-direct` job, the built wheel, coverage, both Sphinx builds, and built-wheel
installs with each extra combination at minimum and newer dependency versions.

## Architecture

```
src/paftacular/
  __init__.py      public names and __version__ (hatch reads the version from here)
  annotation.py    frozen PafAnnotation: make_* factories, mass/mz/comp/formula, serialize, resolve, to/from_dict
  parser.py        mzPAFParser, parse/parse_single/parse_multi, iter_parse/parse_batch, ParseResult
  syntax.py        annotation_spans(): delimiter-aware comma splitting, incl. nested ProForma in {...}
  errors.py        PafParseError(ValueError) with zero-based position and annotation_index
  resolution.py    resolve(): select the fragment sequence from a full ProForma analyte
  serialization.py versioned to_dict/from_dict (schema_version 1), strict field validation
  conversion.py    to_mzpaf(): peptacular Fragment -> PafAnnotation
  constants.py     enums, grammar regexes (_ATOM_TOKEN, FULL_PAF_PATTERN), InternalSeries, INTERNAL_MASS_DIFFS, MAX_CACHE_SIZE
  util.py          validate_number, format_number, validate_integer, parse_formula
  comps/base.py    Serializable (shared __reduce__), MassProvider, CompositionProvider, ScalableComposition
  comps/ions.py    PeptideIon, InternalFragment, ImmoniumIon, ReferenceIon, NamedCompound, ChemicalFormula, SMILESCompound, UnknownIon, PrecursorIon
  comps/modifiers.py MassError, IsotopeSpecification, NeutralLoss, Adduct
  comps/util.py    formula <-> Counter[ElementInfo] helpers
  mcp/             optional MCP server: cli.py, server.py (tools/resources/prompts), handlers.py, models.py (pydantic), guidance.py
```

Data flow: text -> `syntax.annotation_spans` finds annotation boundaries -> each span is
matched whole against `FULL_PAF_PATTERN` -> `mzPAFParser._build_annotation` builds the ion
component plus modifier tuples -> frozen `PafAnnotation`. Calculation sums the ion component
(offset only, or sequence residues when an embedded or resolved sequence exists), neutral
losses, isotopes and adducts, then charge. `serialize()` rebuilds text from the components.

## Public API

Everything is exported from `paftacular/__init__.py`:

- **Parsing:** `parse(s)` returns one `PafAnnotation` or a list (a list for comma input and
  for `""`). `parse_single(s)` requires exactly one. `parse_multi(s)` always returns a list.
  `iter_parse(records)` yields `ParseResult` lazily. `parse_batch(records)` collects them.
  `mzPAFParser` is the class behind these. `PafParseError`, `ParseResult` (`.ok`, `.index`,
  `.text`, `.annotations`, `.error`).
- **Annotation:** `PafAnnotation` (frozen, hashable). Factories `make_peptide`,
  `make_internal`, `make_immonium`, `make_reference`, `make_named_compound`, `make_formula`,
  `make_smiles`, `make_unknown`, `make_precursor`, all taking `CommonAnnotationParams`
  (`neutral_losses`, `isotopes`, `adducts` as strings, `charge`, `mass_error`,
  `mass_error_unit`, `confidence`, `is_auxiliary`, `analyte_reference`). Methods `mass`, `mz`,
  `comp`, `dict_composition`, `formula`, `proforma_formula`, `serialize`, `as_dict`,
  `to_dict`, `from_dict`, `resolve`, `parse`. Properties `sequence`, `peptacular_ion_type`.
- **Ion components** (`ann.ion_type`): `PeptideIon`, `InternalFragment`, `ImmoniumIon`,
  `ReferenceIon`, `NamedCompound`, `ChemicalFormula`, `SMILESCompound`, `UnknownIon`,
  `PrecursorIon`. `IonType` is their union type alias.
- **Modifiers:** `NeutralLoss`, `IsotopeSpecification`, `Adduct`, `MassError`.
- **Enums and tables:** `IonSeries` (a b c d v w x y z da db wa wb), `BackboneCleavageType`,
  `AnnotationName`, `AminoAcids`, `INTERNAL_MASS_DIFFS` (spec section 4.4.4 table).
  `InternalSeries` (ax..cz) lives in `paftacular.constants` and is not exported.
- **Other:** `resolve(annotation, analytes)`, `to_mzpaf(fragment, ...)` (needs peptacular).

## Conventions

- Docstrings are one-line summaries, sometimes with an `Examples:` block. There is no
  enforced Sphinx or Google style. `util.parse_formula` is the only `Args:`/`Returns:` one.
- Full type annotations, checked with `ty check src`. Ships `py.typed`. Ruff line length 160,
  rules E W F I B UP (E741 ignored).
- Components and `PafAnnotation` are frozen dataclasses. Transformations return new objects.
- Errors: parse failures raise `PafParseError` (a `ValueError`). Semantic failures raise
  `ValueError`. Unsupported calculations (unknown ions, named compounds) raise
  `NotImplementedError`. Missing optional dependencies raise `ImportError` naming the extra.
- No logging in the core. `mcp/handlers.py` logs unexpected failures to stderr only.
- Tests live in `tests/test_*.py`, flat. Optional-dependency tests use
  `pytest.importorskip`, so run with all extras (`just install`) to exercise everything.
  MCP tests use the SDK's in-memory Client and real stdio subprocesses outside the checkout.
- Never use em dashes or semicolons in authored text or code comments.

## Gotchas

Scientific conventions:

- `mass()` returns the charged-species mass. `comp()` counts nuclei. To compare summed element
  masses with `mass()`, subtract charge times the electron mass. Upstream ion offsets are
  rounded, so independent checks need an absolute tolerance near one microdalton, with
  relative tolerance disabled.
- Without an embedded or resolved sequence, peptide, internal and precursor calculations
  return only the ion offset and modifiers (`y5` gives 19.0178, `b5` gives 1.0073). Preserve
  this. Use `resolve()` to select the complete fragment sequence from an analyte.
- The `{...}` in `y3{PEP}` is the fragment's own sequence. It is used verbatim and is not
  checked against the position: `y3{PEPTIDE}` computes the mass of all seven residues.
- Embedded sequences contribute residue mass and composition through peptacular
  `ion_type="n"`, because the default precursor composition would add an extra water on top
  of paftacular's own ion offset.
- Formula ions already specify every atom in the charged species. Their adducts label charge
  carriers and add no atoms. Other ion types add adduct atoms and subtract electron mass.
  Implicit protonation and explicit `[M+H]` must agree.
- A generic `+i` uses the carbon-13 minus carbon-12 shift. Add contributions to counters
  rather than overwriting them when isotope keys coincide. At the annotation level, consume
  ordinary monoisotopic atoms when available and keep genuine deficits when sequence context
  is absent. `+iA` (average isotopomer) raises `ValueError` on mass and composition.
- `formula()` raises `ValueError` when the composition mixes signs (for example an unresolved
  `y5+i` or `y5-[Adenine]`). `dict_composition()` returns the signed counts and
  `proforma_formula()` writes negative counts (`C5H11N2O-1`).
- Counter addition and unary plus discard negative entries. Accumulate with `update()` and
  preserve genuine negative counts from modifications.
- Immonium ions use tacular's `"i"` fragment lookup key, not the internal `"by"` key.
- Explicit internal cleavage fields must be supplied together, in a/b/c and x/y/z
  combinations. Their physical composition comes from their own tacular key. Serialize it as
  signed elemental gains or losses, preserving mass.
- The mzPAF 1.0.1 section 4.4.4 correction table disagrees with the physical tacular
  definitions for several combinations. `INTERNAL_MASS_DIFFS` and `make_internal(ion_type=)`
  keep the specification convention (`bx` gives `m2:4+CO`). Conversion and explicit cleavage
  fields keep the source's physical composition (physical `ax` gives `m2:4-H2`). Do not
  silently equate the two. `docs/usage.rst` documents the difference.

Parser and interchange:

- Locate commas outside labels and sequences (`syntax.annotation_spans`) before matching whole
  annotations. Never split blindly on commas or let a sequence consume another annotation.
- Keep the possessive quantifiers (`*+`) in `_ATOM_TOKEN`. Removing them brings back
  catastrophic backtracking on invalid long element runs. A regression test uses SIGALRM and
  is skipped on Windows.
- A leading `+` on a mass error (`y5/+1.2ppm`) and spaces or commas inside named labels
  (`_{A, B}`) are intentional compatibility choices. Keep them accepted.
- Use `format_number()` for grammar-compatible decimal output without truncating float
  precision. Scientific notation is not part of the grammar (`y5/1e-3` is rejected). Mass
  losses print at least five decimals (`y5-17.03` serializes as `y5-17.03000`).
- Component parsers consume the entire input. Zero charge, nonintegral positions and reversed
  ranges must fail. Resolution checks bounds against the analyte.
- `to_dict()` / `from_dict()` use schema version 1 and reject unknown fields. Keep `as_dict()`
  compatible as a display format. Resolved context is stored in structured export and is
  deliberately absent from mzPAF text.
- Optional dependency paths must fail with actionable install instructions. Importing and
  parsing the core must work without any extra.

Caching:

- `ImmoniumIon`, `ReferenceIon`, `NamedCompound`, `UnknownIon`, `PrecursorIon`,
  `IsotopeSpecification`, `NeutralLoss` and `Adduct` cache constructor instances in a
  class-level dict bounded by `MAX_CACHE_SIZE` (10 000, FIFO eviction). `PeptideIon`,
  `InternalFragment` and `PafAnnotation` are not cached. Validate arguments before caching so
  an invalid call cannot mutate a cached object. `Serializable.__reduce__` rebuilds through the
  constructor for pickle and copy. Measure with `just benchmark` before expanding caching.

MCP (`paftacular[mcp]`, console script `paftacular-mcp`, also `python -m paftacular.mcp`):

- Keep SDK and pydantic imports inside `mcp/`, loaded only in `create_server()`. `--help` and
  `--version` and `import paftacular.mcp` must work without the SDK. There is no `--check`.
- Nine tools (`get_capabilities`, `parse_annotations`, `build_annotation`,
  `resolve_annotation`, `serialize_annotation`, `calculate_ion`, `calculate_ions`,
  `generate_fragments`, `match_mz`), four resources (`paftacular://guide`, `conventions`,
  `examples`, `capabilities`) and two prompts (`analyze_fragment`, `review_annotations`).
- stdio only: protocol on stdout, diagnostics on stderr. Tool wrappers are async with bounded
  synchronous handlers and no inner awaits, which serializes chemistry calls and avoids
  racing the component caches. Request/response contracts are versioned separately from core
  interchange. Limits: 100 records, 16 KiB per text, 256 KiB input, 512 KiB output, 1000
  parsed annotations per call. `docs/mcp.rst` is the contract. `MCP_PLAN.md` is historical.

## Releasing

Only the tacular-omics overseer bumps versions or publishes. See `just --list` (`set-version`,
`sync-version`, `check-version`, `pre-release`) and the workspace release checklist. The
version lives in `src/paftacular/__init__.py` (`__version__`, read by hatchling) and is copied
to `CITATION.cff` by `scripts/release_version.py`. Keep `CHANGELOG.md` entries under
`[Unreleased]` until a release is cut. `publish.yml` uploads to PyPI on a published GitHub
release whose tag matches the package version.

## Workspace note

This repo is also developed inside the tacular-omics uv workspace. There `uv run` uses the
shared `.venv` and the local sibling checkouts (tacular, peptacular). See the workspace
CLAUDE.md. Use `just isolated paftacular` from the workspace root to reproduce this repo's CI
with its own `uv.lock` and PyPI siblings.
