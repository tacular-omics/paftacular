# paftacular agent guide

## Project overview

paftacular parses and serializes HUPO-PSI mzPAF 1.0.1 peak annotations (`y5-H2O^2/1.2ppm*0.95`)
and calculates ion masses, m/z and elemental compositions. Users are people reading or writing
fragment annotations from spectral libraries and search-engine output. Python 3.12 or newer.
Imported as `pft` by convention. Current release: 1.4.0. `main` carries the unreleased 2.0.0
breaking major (see `docs/migration.rst`).

Place in the tacular-omics graph (tier 1):

- **Upstream:** `tacular` (required, `>=2.0,<3`) supplies element, amino-acid, modification
  and fragment-ion lookups. `peptacular` (`>=5.0,<6`, extra `peptacular`) is optional and
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
  parser.py        parse/parse_multi, iter_parse, ParseResult, bounded per-substring component caches
  syntax.py        annotation_spans(): delimiter-aware comma splitting, incl. nested ProForma in {...}
  errors.py        PaftacularError(ValueError) base, PafParseError (zero-based position), PafUnknownReferenceError
  resolution.py    resolve(): select the fragment sequence from a full ProForma analyte
  serialization.py versioned to_dict/from_dict (schema_version 1), strict field validation
  conversion.py    to_mzpaf(): peptacular Fragment -> PafAnnotation
  constants.py     enums, grammar regexes (_ATOM_TOKEN, FULL_PAF_PATTERN), InternalSeries, _INTERNAL_MASS_DIFFS, MAX_CACHE_SIZE
  util.py          validate_number, format_number, validate_integer, parse_formula, to_enum
  comps/base.py    Serializable (shared __reduce__), MassProvider, CompositionProvider, ScalableComposition
  comps/ions.py    PeptideIon, InternalFragment, ImmoniumIon, ReferenceIon, NamedCompound, ChemicalFormula, SMILESCompound, UnknownIon, PrecursorIon
  comps/modifiers.py MassError, IsotopeSpecification, NeutralLoss, Adduct
  comps/util.py    formula <-> Counter[ElementInfo] helpers
  mcp/             optional MCP server: cli.py, server.py (tools/resources/prompts), handlers.py, models.py (pydantic), guidance.py
```

Data flow: text -> `syntax.annotation_spans` finds annotation boundaries -> each span is
matched whole against `FULL_PAF_PATTERN` -> `parser._build_annotation` builds the ion
component plus modifier tuples -> frozen `PafAnnotation`. Calculation sums the ion component
(offset only, or sequence residues when an embedded or resolved sequence exists), neutral
losses, isotopes and adducts, then charge. `serialize()` rebuilds text from the components.

## Public API

Everything is exported from `paftacular/__init__.py`:

- **Parsing:** `parse(s)` returns exactly one `PafAnnotation` (raises otherwise).
  `parse_multi(s)` always returns a list. `iter_parse(records)` yields `ParseResult` lazily.
  `PaftacularError`, `PafParseError`, `ParseResult` (`.ok`, `.index`,
  `.text`, `.annotations`, `.error`).
- **Annotation:** `PafAnnotation` (frozen, hashable). Factories `make_peptide`,
  `make_internal`, `make_immonium`, `make_reference`, `make_named_compound`, `make_formula`,
  `make_smiles`, `make_unknown`, `make_precursor`, all taking `CommonAnnotationParams`
  (`neutral_losses`, `isotopes`, `adducts` as strings, `charge`, `mass_error`,
  `mass_error_unit`, `confidence`, `is_auxiliary`, `analyte_reference`). Methods `get_mass`,
  `mz`, `comp`, `formula`, `proforma_formula`, `serialize`, `to_dict`, `from_dict`,
  `resolve`, `parse`. Optional arguments are keyword-only everywhere. Properties `sequence`, `peptacular_ion_type`.
- **Ion components** (`ann.ion_type`): `PeptideIon`, `InternalFragment`, `ImmoniumIon`,
  `ReferenceIon`, `NamedCompound`, `ChemicalFormula`, `SMILESCompound`, `UnknownIon`,
  `PrecursorIon`. `IonType` is their union type alias.
- **Modifiers:** `NeutralLoss`, `IsotopeSpecification`, `Adduct`, `MassError`.
- **Enums and tables:** `IonSeries` (a b c d v w x y z da db wa wb), `BackboneCleavageType`,
  `AnnotationName`. Immonium amino acids are `tacular.AminoAcid` (20 standard codes only). The spec section 4.4.4 table is private
  (`constants._INTERNAL_MASS_DIFFS`).
  `InternalSeries` (ax..cz) lives in `paftacular.constants` and is not exported.
- **Other:** `resolve(annotation, analytes)`, `to_mzpaf(fragment, ...)` (needs peptacular).

## Conventions

- Docstrings are one-line summaries, sometimes with an `Examples:` block. There is no
  enforced Sphinx or Google style. `util.parse_formula` is the only `Args:`/`Returns:` one.
- Full type annotations, checked with `ty check src`. Ships `py.typed`. Ruff line length 160,
  rules E W F I B UP (E741 ignored).
- Components and `PafAnnotation` are frozen dataclasses. Transformations return new objects.
- Errors: every error from user input is a `PaftacularError` (a `ValueError`). Parse failures
  raise its subclass `PafParseError`. Wrap errors from tacular and peptacular. Unsupported calculations (unknown ions, named compounds) raise
  `PafUnsupportedCalculationError`. Missing optional dependencies raise `ImportError` naming the extra.
- No logging in the core. `mcp/handlers.py` logs unexpected failures to stderr only.
- Tests live in `tests/test_*.py`, flat. Optional-dependency tests use
  `pytest.importorskip`, so run with all extras (`just install`) to exercise everything.
  MCP tests use the SDK's in-memory Client and real stdio subprocesses outside the checkout.
- Never use em dashes or semicolons in authored text or code comments.

## Gotchas

Scientific conventions:

- `get_mass()` returns the charged-species mass. `comp()` counts nuclei. To compare summed element
  masses with `get_mass()`, subtract charge times the electron mass. Upstream ion offsets are
  rounded, so independent checks need an absolute tolerance near one microdalton, with
  relative tolerance disabled.
- Named modifications (Unimod, PSI-MOD, ...) add their listed database mass through
  peptacular's default mass path, never `calculate_with_composition=True`, the same rule as
  peptacular. Plain fragment and precursor ions match peptacular to 1e-9 Da
  (`tests/test_peptacular_agreement.py`), and so do ions with neutral losses or isotope
  peaks (peptacular 5 keeps listed masses there too). Unimod reference names in `lookup_reference` keep the listed mass too; mzPAF
  reference-list entries keep exact formula masses. Composition is only the fallback, and
  is used under a global isotope label (`<13C>`) in both packages. Labile mods count only
  for precursor ions (fragments lose them, like peptacular).
- Without an embedded or resolved sequence, peptide, internal and precursor calculations
  return only the ion offset and modifiers (`y5` gives 19.0178, `b5` gives 1.0073). Preserve
  this. `mz()` raises for those instead, because an offset has no meaningful m/z. Use `resolve()` to select the complete fragment sequence from an analyte.
- The `{...}` in `y3{PEP}` is the fragment's own sequence. It is used verbatim. A length
  that differs from the position (`y3{PEPTIDE}`) makes `get_mass()`/`comp()` emit a `UserWarning`
  (mzPAF 4.4.3: MUST NOT be shorter, SHOULD NOT be longer) and still uses all seven residues.
- Embedded sequences contribute residue mass and composition through peptacular
  `ion_type="n"`, because the default precursor composition would add an extra water on top
  of paftacular's own ion offset.
- Formula ions already specify every atom in the charged species. Their adducts label charge
  carriers and add no atoms. Other ion types add adduct atoms and subtract electron mass.
  Implicit protonation and explicit `[M+H]` must agree exactly: an `H` carrier whose sign
  matches the charge is charged as `PROTON_MASS` (monoisotopic) or natural-abundance H less
  an electron (average). `[M+H]^-1` is a hydride, an H atom plus an electron.
- A global isotope label in the embedded sequence relabels the ion offset and formula deltas
  too (peptacular 5 does the same), and so do atoms removed by carriers: the H of a default
  negative charge and negative adducts (`[M-H]`). Mass-only deltas, isotope shifts, added
  adducts and the positive charge proton stay unlabelled. Global fixed modifications on v/w/d follow the explicit-mod rule.
- A generic `+i` uses the carbon-13 minus carbon-12 shift. Add contributions to counters
  rather than overwriting them when isotope keys coincide. At the annotation level, consume
  ordinary monoisotopic atoms when available and keep genuine deficits when sequence context
  is absent. `+iA` (average isotopomer) raises `PaftacularError` on mass and composition.
- `formula()` raises `PaftacularError` when the composition mixes signs (for example an
  unresolved `y5+i` or `y5-[Adenine]`). `comp()` returns the signed counts and
  `proforma_formula()` writes negative counts (`C5H11N2O-1`).
- Counter addition and unary plus discard negative entries. Accumulate with `update()` and
  preserve genuine negative counts from modifications.
- Immonium ions use tacular's `"i"` fragment lookup key, not the internal `"by"` key.
- Explicit internal cleavage fields must be supplied together, in a/b/c and x/y/z
  combinations. Their physical composition comes from their own tacular key. Serialize it as
  signed elemental gains or losses, preserving mass.
- The mzPAF 1.0.1 section 4.4.4 correction table disagrees with the physical tacular
  definitions for several combinations. `_INTERNAL_MASS_DIFFS` and `make_internal(ion_type=)`
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
- Charge is a nonzero integer. Negative charge removes protons and serializes as `^-n`
  (`serialize(signed_charge=False)` writes `^n`). `mz()` divides by the absolute charge.
- Component parsers consume the entire input. Zero charge, nonintegral positions and reversed
  ranges must fail. Resolution checks bounds against the analyte.
- `to_dict()` / `from_dict()` use schema version 1 and reject unknown fields. Resolved context is stored in structured export and is
  deliberately absent from mzPAF text.
- Optional dependency paths must fail with actionable install instructions. Importing and
  parsing the core must work without any extra.

Caching:

- Constructors build fresh objects. There is no `__new__` interning. The parser keeps bounded
  per-substring caches (`parser._ION_CACHE`, `_LOSS_CACHE`, `_ISOTOPE_CACHE`, `_ADDUCT_CACHE`,
  `_MASS_ERROR_CACHE`, each `MAX_CACHE_SIZE` = 10 000 with FIFO eviction), so repeated text
  shares one immutable component. `parser._clear_caches()` empties them (tests use it). Never
  mutate a cached component: `SMILESCompound.composition` is a `cached_property` on a shared
  object. `conversion.py` caches per-ion-type plans, peptide ions and loss units with
  `lru_cache`. `Serializable.__reduce__` rebuilds through the constructor for pickle and copy.
  Measure with `just benchmark` before expanding caching.

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
