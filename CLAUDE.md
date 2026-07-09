# paftacular — Claude Code Guide

## Project Overview

paftacular is a Python library for parsing and serializing **mzPAF**, the
HUPO-PSI Peak Annotation Format — a compact, human-readable notation for
mass spectrometry fragment-ion annotations used in peptide/proteomics
analysis (e.g. `y5-H2O^2/1.2ppm*0.95`, `b2{PEP}`, `IK[Acetyl]`). It handles
both directions: parsing mzPAF strings into structured `PafAnnotation`
objects, and building/serializing annotations back to mzPAF text, plus
computing their mass/elemental composition/chemical formula.

paftacular requires `tacular` (same author, sibling repo) as its only
runtime dependency — `tacular` supplies the element, amino-acid, and
fragment-ion-mass/composition lookup tables (`ELEMENT_LOOKUP`, `AA_LOOKUP`,
`FRAGMENT_ION_LOOKUP`, `REFMOL_LOOKUP`) that all mass/composition
calculations in this package are built on top of. `peptacular` (also same
author) is an **optional** dependency, only needed for peptide-sequence-aware
features: resolving an embedded ProForma `sequence=` on a fragment ion into
its own mass/composition contribution, and `conversion.to_mzpaf()`, which
converts a `peptacular` `Fragment` object directly into a `PafAnnotation`.
Every module that touches `peptacular` guards it with the same pattern —
see "Optional peptacular dependency" below — so importing paftacular without
`peptacular` installed still works for pure mzPAF parsing/serialization.

## Commands

```bash
just install         # uv sync --all-extras (dev + all optional extras)
just install-prod     # uv sync --no-dev --frozen (production install)
just lint              # ruff check src
just format             # ruff isort + format (src, tests)
just ty                  # ty type check src
just check                # lint + ty + test
just test                  # pytest tests
just test-cov                # pytest with branch coverage (term + html + xml)
just codecov-tests             # pytest --cov + junitxml, for CI/Codecov upload
just docs                        # sphinx-build docs -> docs/_build/html
just docs-test                     # sphinx doctest build (runs every testcode/testoutput block)
just docs-clean                      # rm -rf docs/_build
just docs-open                         # build docs and open in browser
just pre-release                         # format + lint + check + test + docs-test
```

## Architecture

```
src/paftacular/
  __init__.py          # public API re-exports: PafAnnotation, ion/modifier types,
                        # parse/parse_multi/parse_single, to_mzpaf
  annotation.py          # PafAnnotation: the central frozen dataclass, wraps an
                          # IonType plus common modifiers (neutral_losses, isotopes,
                          # adducts, charge, mass_error, confidence). Factory methods:
                          # make_precursor/make_peptide/make_internal/make_immonium/
                          # make_reference/make_named_compound/make_formula/make_smiles/
                          # make_unknown. Instance methods: mass()/mz()/comp()/formula()/
                          # proforma_formula()/serialize()/as_dict()/dict_composition()
  parser.py                # mzPAFParser (singleton via __new__, cls._instance) — the
                            # regex-driven parser built on constants.FULL_PAF_PATTERN;
                            # module-level parse()/parse_multi()/parse_single() wrap it
  conversion.py               # to_mzpaf(): converts a peptacular Fragment into a
                              # PafAnnotation (gated behind _require_peptacular())
  constants.py                  # IonSeries/InternalSeries/AminoAcids/BackboneCleavageType
                                 # enums; MAX_CACHE_SIZE; the regex fragments assembled
                                 # into FULL_PAF_PATTERN; INTERNAL_SERIES_TO_DIFF /
                                 # INTERNAL_MASS_DIFFS (the ax/bx/cx/ay/by/... backbone
                                 # cleavage -> neutral-loss-formula table from the spec)
  util.py                        # parse_formula(): low-level "C6H12O6" / "[13C2]H6"
                                  # string -> Counter[str] parser (isotope-aware)
  comps/
    base.py                        # Serializable / MassProvider / CompositionProvider /
                                    # ScalableComposition — the shared protocols every
                                    # ion type and modifier type implements
    ions.py                          # PeptideIon, InternalFragment, ImmoniumIon,
                                      # ReferenceIon, NamedCompound, ChemicalFormula,
                                      # SMILESCompound, UnknownIon, PrecursorIon
                                      # (the mzPAF "ion_type" alternatives; IonType is a
                                      # union of all of them)
    modifiers.py                       # MassError, IsotopeSpecification, NeutralLoss,
                                        # Adduct
    util.py                              # formula_to_composition(),
                                          # composition_to_formula_string(),
                                          # composition_to_proforma_formula_string()
```

## Gotchas (read before touching mass/composition code or adding a new ion type)

- **`FRAGMENT_ION_LOOKUP` key correctness is easy to get wrong, and wrong keys
  fail silently (wrong number, not an exception).** Each ion type must look
  up *its own* backbone-cleavage/ion-series key, not a convenient default.
  This session fixed two real instances: `ImmoniumIon.mass()`/`.composition`
  was reading `FRAGMENT_ION_LOOKUP["by"]` (the internal by-fragment shift,
  `0`) instead of the immonium-specific `"i"` key; and `InternalFragment`
  ignored its own `nterm_ion_type`/`cterm_ion_type` fields entirely, always
  using the default `by` shift instead of deriving the key via
  `_fragment_ion_key` (see `comps/ions.py`). When adding a new ion type or
  wiring up a new lookup, double-check the lookup key matches the ion's
  actual series/cleavage type — there is no test that catches a
  plausible-but-wrong key by construction, only tests that assert specific
  expected values.
- **`InternalFragment.nterm_ion_type`/`cterm_ion_type` must be set together
  or not at all** — `__post_init__` raises `ValueError` if only one is set.
  If you add a new caller that sets one of these fields, always set the
  other (see `_fragment_ion_key` in `comps/ions.py`, which defaults both to
  `IonSeries.B`/`IonSeries.Y` only when *both* are `None`).
- **`Counter` arithmetic (`+`, `+=`, unary `+`) silently drops any key whose
  total is `<=0`** — including net-negative totals that are physically real
  when a modification removes more of an element than the base residue
  supplies (e.g. `IK[Dehydrated]`). Composition accumulation must use
  `.update()` (which never filters) throughout and **must NOT finish with
  unary `+`** — that reintroduces a silent `.composition`/`.mass()`
  disagreement. `ImmoniumIon.composition` (`comps/ions.py`) strips only
  exact-zero entries at the end (`{el: n for el, n in c.items() if n != 0}`)
  and keeps net-negatives; `PafAnnotation.comp()` returns the raw
  `.update()`-accumulated `Counter` (negatives kept). A mixed-sign
  composition can't be rendered by `composition_to_formula_string` (it
  raises, correctly — no plain formula can express `-1` of an element), but
  `composition_to_proforma_formula_string` handles negatives (`C5H11N2O-1`).
- **Instance-caching pattern used throughout `comps/`**: `ImmoniumIon`,
  `ReferenceIon`, `NamedCompound`, `UnknownIon`, `PrecursorIon`, `NeutralLoss`,
  `Adduct`, and `IsotopeSpecification` are all frozen dataclasses that
  override `__new__` to memoize instances in a `ClassVar[dict]` keyed by
  their constructor args, evicting the oldest entry once
  `constants.MAX_CACHE_SIZE` (10,000) is reached. This makes repeated
  parsing/construction of the same annotation components cheap and lets
  equal annotations compare `is`-identical. If you add a new ion/modifier
  type that should participate in this caching, copy the exact
  `_cache: ClassVar[dict[tuple, "T"]] = {}` + `__new__` + eviction pattern
  from an existing class (e.g. `ImmoniumIon`) rather than inventing a new
  scheme.
- **Optional `peptacular` dependency pattern**: `annotation.py`,
  `comps/ions.py`, and `conversion.py` each do
  `try: import peptacular as pt \n except ImportError: pt = None` at module
  level, and gate any code path that actually needs `pt` behind a local
  `_require_peptacular()` helper (raises a clear `ImportError` pointing at
  `pip install paftacular[peptacular]`). Follow this exact pattern for any
  new peptacular-dependent code — don't import `peptacular` unconditionally
  at module scope, and don't assume it's present without calling
  `_require_peptacular()` first.
- **`pyproject.toml` currently pins `peptacular` to a local editable path**
  (`[tool.uv.sources] peptacular = { path = "../peptacular", editable = true }`)
  even though the version floor (`peptacular>=3.1.2`) is set correctly. This
  is intentional/temporary: 3.1.2 fixes composition/mass bugs paftacular
  depends on but hadn't been published to PyPI as of the `1.1.0` release —
  the override should be removed once that version is actually on PyPI.
  Don't "fix" this by deleting the override without checking PyPI first.
- **Isotope mass/composition must agree, and a generic `+i` means a `13C`
  substitution.** Per mzPAF §4.6 a generic isotope (`+i`, `+2i`, no element)
  has the theoretical mass of the `13C`−`12C` difference; `IsotopeSpecification`
  therefore resolves `element is None` to `{13C: +count, 12C: -count}` in
  **both** `.mass()` and `.composition` (they must never disagree — a prior
  bug had `.mass()` return the shift while `.composition` raised). An average
  isotope (`+iA`) has no defined monoisotopic value, so both `.mass()` and
  `.composition` raise — consistently. `PafAnnotation.mass()`/`.comp()` apply
  an isotope loop over `self.isotopes`; any new `IsotopeSpecification` state
  must keep mass and comp consistent (both return a value, or both raise).
- **`ChemicalFormula` (`f{...}`) ions are charged-species, not neutral.** Per
  mzPAF §4.4.9 the braced formula already counts every nucleus in the charged
  ion, so `PafAnnotation.mass()`/`.comp()` special-case `ChemicalFormula`:
  charge is corrected by **one electron mass per charge** (`0.000548579909`),
  **not** a proton, and adducts add **no** mass/atoms (they only label which
  atoms carry the charge). Don't apply the default proton-per-charge path to
  formula ions. Verified against the spec's worked `f{...}` m/z examples.
- **`comp()` is neutral-basis; `mass()` is the charged species.** `comp()`
  adds a full `H` proton per charge, so its summed element masses are one
  electron *heavier* per charge than `mass()`. The invariant is
  `mass() == sum(comp element masses) - charge*electron_mass`. This gap is
  by-design for every ion type (except the electron-only `ChemicalFormula`
  path) — don't "fix" it, and account for it in any mass-vs-comp test.
- **`_ATOM_TOKEN` uses possessive quantifiers (`*+`) on purpose** — its
  `[A-Z][A-Za-z0-9]*+` under the surrounding `_ATOM_TOKEN+` is a classic
  `(a+)+` catastrophic-backtracking (ReDoS) shape, and a long single-letter
  run through the anchored `mzPAFParser.parse` (e.g. `y1+HHHH…!`) would hang.
  Don't relax `*+` back to `*`. (The module-level `parse`/`parse_multi`, which
  use the unanchored `PARTIAL_PAF_PATTERN`, were never vulnerable.)
- **Comma-separated annotations are split via the spec's Appendix A greedy
  match, not `str.split(",")`.** `parse_multi` matches one annotation at the
  current offset with the unanchored `PARTIAL_PAF_PATTERN`, then requires a
  comma or end-of-string (see `parser.py`); `parse()` delegates to it and both
  build via `_build_annotation`. This is required because commas legitimately
  appear inside bracketed labels (references `r[Ref,X]`, named compounds).
  Don't reintroduce a naive comma split.
- **Some parser leniencies are deliberate spec choices — don't "tighten" them
  back.** The mzPAF spec's normative prose and its published reference regex
  disagree in a few places, and paftacular follows the reference regex: mass
  error accepts an optional leading `+` (`/+0.5ppm`), and named compounds
  `_{...}` allow spaces (matching the spec's own `_{Urocanic Acid}` example)
  and commas (benign — the `}` boundary stops the greedy match before any real
  separator). These are asserted by tests; changing them will break those.
- **mzPAF only standardizes the default internal fragment ion**
  (`m<start>:<end>`, the `by` backbone cleavage). The other 8 combinations in
  `constants.INTERNAL_MASS_DIFFS`/`INTERNAL_SERIES_TO_DIFF` come from the
  spec's own table of neutral-loss corrections relative to `by` — if you
  touch these, cross-check against the actual mzPAF grammar/spec rather than
  guessing values.
- This repo's remote is `https://github.com/tacular-omics/paftacular` (moved
  from a prior `pgarrett-scripps/paftacular` remote — badge/link URLs across
  docs were updated to match in the same session that produced this file;
  watch for stale `pgarrett-scripps` URLs creeping back into new docs/README
  edits).

## Testing

- Tests live in `tests/`, one file roughly per module/feature area:
  `test_annot.py` (PafAnnotation), `test_caching.py` (the `__new__`-based
  instance caching described above), `test_comps_parse.py` (ion/modifier
  `.parse()` methods), `test_conversion.py` (`to_mzpaf()`),
  `test_make_methods.py` (the `PafAnnotation.make_*` factories),
  `test_parser_types.py` (`mzPAFParser`/`parse`/`parse_multi`/`parse_single`),
  `test_util.py` (`parse_formula` and `comps/util.py` conversions).
- `just docs-test` doubles as an executable-example test suite — every
  `.. testcode::`/`.. testoutput::` block in `docs/*.rst` is run and its
  output compared exactly (no float rounding tolerance), so a mass/formula
  calculation change that shifts a printed value will show up there too.

## Release process

- Version lives in `src/paftacular/__init__.py` (`__version__`), sourced by
  `[tool.hatch.version]` in `pyproject.toml`.
- Changelog is `HISTORY.md` (`## X.Y.Z (YYYY-MM-DD)` sections, terse bullet
  points describing user-visible fixes/changes).
- `just pre-release` runs format + lint + check + test + docs-test — run
  this (not just `just check`) before cutting a release, since a docs
  example can go stale without any src/tests change flagging it.
