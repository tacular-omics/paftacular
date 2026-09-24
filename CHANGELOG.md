# Changelog

## [Unreleased]

### 2.0.0 (unreleased, breaking)

Every rename and removal, with the replacement, is in the
[migration guide](https://paftacular.readthedocs.io/en/latest/migration.html)
(`docs/migration.rst`).

#### Removed

- `parse_single` (use `parse`), `parse_batch` (use `list(iter_parse(...))`), `mzPAFParser`
  and `paftacular.parser.MZ_PAF_PARSER` (use the module functions).
- `mass()` on `PafAnnotation`, every ion component and every modifier (use `get_mass()`).
- `dict_composition()` (use `comp()`) and `as_dict()` on `PafAnnotation`, `NeutralLoss`,
  `IsotopeSpecification` and `Adduct` (use `to_dict()`).
- `mz(calculate_sequence=...)`. `mz()` now takes only `monoisotopic`.
- The `hill_order` argument of `composition_to_proforma_formula_string`, which had no effect.
- The class-level `_cache` dicts and `__new__` interning on ion components and modifiers.
- `AminoAcids` (top level and `paftacular.constants`). Use `tacular.AminoAcid`.
  `ImmoniumIon.amino_acid` is now a `tacular.AminoAcid`, and immonium ions still accept
  only the 20 standard codes.
- `INTERNAL_MASS_DIFFS` (top level and `paftacular.constants`) and
  `paftacular.constants.INTERNAL_SERIES_TO_DIFF` are now private. `make_internal(ion_type=)`
  applies the table.

#### Changed

- MCP `match_mz`: `tolerance_unit` is tacular's `ToleranceUnit`, `"da"` or `"ppm"`
  (default `"ppm"`). `"Th"` is now `"da"` (an absolute m/z difference) and is rejected.
  The MCP `response_schema_version` is 2.
- `parse(s)` returns exactly one `PafAnnotation` and raises `PafParseError` for comma input
  or empty text. `parse_multi(s)` always returns a list.
- `get_mass(*, monoisotopic=True)` replaces `mass()`, matching tacular 2.0. Optional
  arguments are keyword-only on component constructors, the `make_*` factories, `comp`,
  `formula`, `proforma_formula`, `serialize`, `NeutralLoss.serialize`, `format_number`,
  `validate_integer` and `to_mzpaf`.
- `to_mzpaf(include_annotation=)` is renamed `include_sequence=`.
- `to_mzpaf(mass_error_type=)` is renamed `mass_error_unit=`.
- `PafAnnotation(...)` fields after `ion_type` are keyword-only.
- The mass or composition of a `?` or `_{...}` ion raises the new
  `PafUnsupportedCalculationError` (a `PaftacularError`) instead of `NotImplementedError`.
- Formula parsing accepts ASCII digits only. `ChemicalFormula("H²O").get_mass()` raises
  `PaftacularError`. The annotation grammar does too, so `y٢{DE}` raises `PafParseError`.
- An adduct or neutral-loss formula with a token that is not an element (`IK[M+Methyl]`,
  `y2{DE}[M+Methyl]`, `y2{DE}-Methyl`) raises `PafParseError` when parsed, not later in
  `get_mass()`.
- `IK[M+K]` parses as the K immonium ion with a K+ adduct, so `serialize()` round-trips.
  1.x read `M+K` as the immonium modification. `ImmoniumIon(modification=)` rejects adduct text.
- A negative-charge `to_dict()` keeps `schema_version` 1, and paftacular 1.4.0 `from_dict`
  rejects it.
- `mz()` raises `PaftacularError` for a peptide, internal or precursor ion without a sequence,
  instead of dividing the ion-type offset by the charge.
- New `PaftacularError(ValueError)` base class. `PafParseError`, `PafUnknownReferenceError`
  and every other error caused by user input use it, including errors from tacular and
  peptacular that 1.x let escape.
- Negative charge is allowed (`y2{DE}^-2`). `serialize()` writes `^-n`, and
  `serialize(signed_charge=False)` writes the magnitude only.
- `to_mzpaf` output, for peptacular 5 fragments: negative charge is kept, known neutral
  deltas use canonical names (`-NH3`, `+HCOOH`, `-HCONH2`), equal mass deltas are folded
  and rounded to 6 decimals, charge carriers are sorted, and a terminal modification on an
  immonium residue becomes the immonium modification (`IP[Acetyl]`) so the mass is kept.
  A global fixed modification is kept the same way (`<[Oxidation]@P>` gives `IP[Oxidation]`),
  and a global isotope label becomes isotope shifts (`<13C>` gives `IP+4i13C`). Internal ion
  offsets use canonical names too (`-NH3`, `-HCONH2`, not `-H3N`, `-CH3NO`). Labels match
  peptacular 5 `Fragment.to_mzpaf()`.
- Every ion-type offset (all peptide series, immonium, internal, precursor) is summed from
  exact element masses, not tacular's 6-decimal constants. Offset masses move by up to
  ~4e-7 Da (y by 3.2e-7, immonium by 3.8e-7).
- A named modification (Unimod, PSI-MOD, RESID, XLMOD, GNO) in an embedded or resolved
  sequence counts at its listed database mass (Oxidation 15.994915), the same rule as
  peptacular. Plain fragment and precursor ions agree with peptacular to 1e-9 Da. Ions with
  neutral losses or isotope peaks will agree once the matching peptacular fix lands.
  Composition is used only for modifications with no listed mass (formulas, glycans), and
  under a global isotope label (`<13C>`, `<15N>`) in both packages.
  A Unimod name used as a reference ion or loss (`r[Hex]`, `-[Hex]`) also uses the listed
  6-decimal mass. mzPAF reference-list entries (`r[TMT6plex]`, `r[iTRAQ4plex]`) keep their
  exact formula masses. `comp()` and `formula()` are unchanged, so for a named modification the mass summed
  from `comp()` can differ from `get_mass()` by up to ~1e-6 Da, because the listed mass is
  rounded. A modified sequence mass is about twice as fast after the first call.
- Monoisotopic charge uses tacular's CODATA `PROTON_MASS` for the default charge and for an
  `H` carrier, so `y2{DE}[M+H]` equals `y2{DE}` exactly (1.x charged `[M+H]` as H less an
  electron, 1.4e-8 Da lighter). An `H` carrier of the opposite sign (`[M+H]^-1`) is a hydride:
  an H atom plus an electron. The average charge is natural-abundance H less an electron.
  For average mass 1.x added the monoisotopic proton, 1.16e-4 Da per charge too light.

#### Fixed

- Labile modifications (`{Glycan:Hex}PEPTIDEK`) are lost on fragmentation, as ProForma
  defines them and peptacular computes them. Fragment ions no longer add their mass
  (`y7^2` of `{Glycan:Hex}PEPTIDEK` was 81 Da too heavy). Precursor ions keep it.
- A global isotope label (`<13C>`) in an embedded sequence now replaces its element in the
  ion offset and in formula deltas too, not only in the residues, like peptacular 5.
  `a2{<13C>RY}` has 14 13C, not 15 (1.x labelled the residues and left the offset C
  light, off by one label shift per offset atom for a, c, x, z, v, w, d, internal and
  precursor ions). Mass-only deltas, isotope shifts, added adducts and the positive charge
  proton stay unlabelled.
- A global fixed modification on the residue whose side chain a side-chain ion loses is
  handled like an explicit one. A v ion loses it with the side chain
  (`v3{<[Carbamidomethyl]@C>CFQ}` is 349.151, 1.x gave 406.172), and w and d ions raise
  `PaftacularError`.
- `to_mzpaf` counts an immonium ion's global isotope label on the final ion, after formula
  deltas and removed carrier atoms. `<15N>K` with `-NH3` gives `IK-NH3+i15N` (one 15N left),
  not `+2i15N`, and `<2H>P` at `^-1` gives `IP+6i2H^-1`.
- A negative charge on a labelled ion removes a labelled atom. Under `<2H>` the default
  negative charge and a `[M-H]` carrier remove a deuteron, so `b2{<2H>PE}^-1` is
  `C10[2H13]N2O4` with no negative H, as in peptacular 5.
- Requires `tacular>=2.0,<3`. The `peptacular`, `mcp` and `all` extras require
  `peptacular>=5.0,<6`.

#### Performance

- Parsing is about 2x faster (a leaner parser, with components shared by
  substring in bounded parser caches).
- `to_mzpaf(...).serialize()` is about 11x faster with the sequence embedded and about 4x
  faster without it (per-ion-type conversion plans, cached peptide ions and loss units).

## [1.4.0] (2026-09-23)

### Added

- `PafUnknownReferenceError`, raised when calculating an `r[...]` ion or a `-[...]` loss whose
  name is in neither the mzPAF reference list nor Unimod. It subclasses both `ValueError` and
  `KeyError` and carries the name in `name`.

### Changed

- Requires `tacular>=1.2,<2`, and the `peptacular`, `mcp` and `all` extras require
  `peptacular>=4.2,<5` (were `tacular>=1.1.0`, `peptacular>=3.1.2`).

Behaviour changes since 1.3.2 that callers can notice:

- `PafAnnotation.peptacular_ion_type` returns `pt.IonType.Z_RADICAL` for a `z` ion (was
  `pt.IonType.Z`), and `to_mzpaf` labels a peptacular `z` fragment as `z3-H` (was `z3`).
- `PeptideIon("d"|"v"|"w", n)` without a sequence returns the section 4.4.3 series
  composition (`C2H4N`, `C2H3NO2`, `C3H4O2`), so `mass()`, `composition` and `formula` differ
  from 1.3.2 (`C2H3N`, `C2H2NO`, `C3H3O`).
- `PeptideIon("da"|"db"|"wa"|"wb", n)` without a sequence still raises `ValueError` from
  `mass()`, `composition` and `formula`, now with "needs a sequence" in the message. With a
  sequence they return the residue-specific value.
- Unknown reference names in `r[...]`, `ReferenceIon` and `-[...]` losses raise
  `PafUnknownReferenceError` from `mass()`, `composition` and `formula`. It is a `ValueError`
  and still a `KeyError`, so `except KeyError` code written for 1.3.2 keeps working. The MCP
  server reports these as `unknown_reference` with the plain message.
  Unimod names that used to raise, such as `r[Hex]`, now resolve.
- `d` and `w` ions on a residue where the series is undefined (G, A, P; plain `d`/`w` on
  T or I; a modified residue n) raise `ValueError` when a sequence is given.
- The source distribution ships only `src/`, `tests/`, `docs/usage.rst` (the tests run its
  examples), `README.md`, `LICENSE`, `CHANGELOG.md`, `CITATION.cff` and `pyproject.toml`
  (about 95 KB, was about 800 KB). The specification PDF, logo, other docs, `uv.lock` and
  agent notes are no longer included. The unused `MANIFEST.in` is removed.
- `scripts/release_version.py sync --set X.Y.Z` also sets `date-released` in `CITATION.cff`.

### Fixed

- Bracketed neutral losses and gains accept reference molecule names containing `-` or `_`,
  such as `p-[TMT126-ETD]`, `a1+[TMTpro_zero]` and `y2-[sidechain_A]` (34 of the 71 mzPAF
  Appendix B names). They raised `PafParseError`. mzPAF 1.0.1 section 4.5 allows any reference
  molecule name as a loss and the section 6.2 grammar allows both characters; the section 6.1
  regex omits them. `r[...]` reference ions already accepted these names.
- Bracketed neutral losses and gains accept names with balanced parentheses, such as
  `y2-[HexNAc(2)]` and `p+[Hex(1)HexNAc(2)]`. They raised `PafParseError`. mzPAF 1.0.1
  section 4.4.7 uses `HexNAc(2)` and the section 6.2 grammar allows `(` and `)`. Formula and
  mass losses parse as before.
- `to_mzpaf` converts peptacular `d`, `v`, `da`, `db`, `wa` and `wb` fragments, including the
  residue-specific `d-valine`, `w-valine`, `da-threonine` and similar types. It raised
  "Cannot convert fragment" because peptacular also sets `AA_SPECIFIC_FWD`/`AA_SPECIFIC_BWD`
  on these series.
- `PeptideIon`, `InternalFragment`, `PrecursorIon` and `ReferenceIon` return a copy from
  `composition`. They returned tacular's cached `Counter`, so mutating the result changed
  every later mass and formula in the process.
- `z` ions now follow mzPAF 1.0.1 section 4.4.3 (the z-dot radical, sum + H2O - NH2).
  They were 1.007825 Da light. `peptacular_ion_type` maps `z` to `IonType.Z_RADICAL`, and
  `to_mzpaf` writes peptacular `z`, `z+H` and `c-H` fragments as `z3-H`, `z3+H` and `c3-H`
  (these used to raise or give the wrong mass).
- `d`, `v` and `w` side-chain ions now use the section 4.4.3 formulas (the other n-1
  residues plus the kept part of residue n). They were off by most of a residue.
  `da`/`db`/`wa`/`wb` apply to T and I only. `d` and `w` raise `ValueError` for G, A, P,
  for plain `d`/`w` on T or I, and for a modified residue n. `v` drops a modification on
  residue n with its side chain. Analyte resolution now supports these series.
- Reference names in `r[...]` and `-[...]` fall back to Unimod entry names (`r[Hex]`,
  `r[HexNAc(2)]`, `p-[Hex]`) as sections 4.4.7 and 4.5 allow. Unknown names raise
  `PafUnknownReferenceError`, a `ValueError` that is also a `KeyError`.
- Tests: every worked example in the mzPAF 1.0.1 specification is parsed and round-tripped,
  and 532 m/z values are checked against a frozen reference built with pyteomics
  (`tests/reference/`).
- Tests: Hypothesis property tests (`tests/test_properties.py`) for the serialise and parse
  round trip, mass and m/z across charge states, adducts and isotopes, and clean errors on
  invalid strings. `HYPOTHESIS_PROFILE=thorough` runs 5000 examples per property.

## [1.3.2] (2026-09-23)

### Fixed

- `mass()`, `mz()`, `comp()` and the formula methods now emit a `UserWarning` when an embedded
  sequence's residue count differs from the ion position (`y3{PEPTIDE}`, `m2:5{PEPTIDE}`).
  mzPAF 1.0.1 section 4.4.3 says it MUST NOT be shorter and SHOULD NOT be longer. The
  calculation is unchanged and still uses every residue.
- README: the MCP-with-SMILES install extra is `paftacular[mcp,smiles]`, not `peptacular[...]`.
- `docs/usage.rst`: the `make_internal(2, 5, sequence=...)` example now embeds the four-residue
  `EPTI` instead of `PEPTIDE`.
- CI: removed the redundant `windows-mcp` job (the `test` matrix already runs the full suite
  with all extras on Windows). Dependabot ignores major updates of the sibling packages
  `tacular` and `peptacular`, which stopped its uv job failing with "Expected lockfile to change!".

## [1.3.1] (2026-09-23)

- Capped sibling requirements below their next major version (`tacular>=1.1.0,<2`, `peptacular>=3.1.2,<5`). Verified against peptacular 4.0.0.
- Publish from GitHub Actions with PyPI trusted publishing (`publish.yml`);
  release metadata is checked against the tag.
- Keep `__version__` and `CITATION.cff` in sync with `scripts/release_version.py`
  (`just set-version X.Y.Z`).
- CI tests Python 3.12-3.14 on Linux plus macOS and Windows, the lowest direct
  dependency versions, and the built wheel.

## [1.3.0] (2026-09-04)

- Added an optional local MCP server with nine tools, four resources, and two analysis prompts
- Added peptide-context calculations, annotation construction, fragment-series generation, and candidate m/z matching through MCP
- Added typed MCP input and output schemas, scientific context labels, partial results, and bounded batches
- Added the `mcp` extra and `paftacular-mcp` command, with peptide support included and MCP included in `all`
- Added MCP protocol, isolated installation, and Windows connection checks

## [1.2.0] (2026-09-04)

### Correctness and maintenance

- Corrected embedded-sequence compositions, explicit-adduct electron accounting, isotope substitution, and SMILES isotope labels
- Preserved physical internal-fragment masses through conversion and serialization, documenting the difference from the specification's cleavage table
- Fixed counted charge-adduct conversion and retained precursor sequence context during peptacular conversion
- Fixed annotation boundaries around embedded ProForma, required complete component parsing, and rejected zero charge and invalid positions
- Preserved float precision in grammar-compatible numeric serialization
- Bounded the reference-ion cache and repaired pickle and copy reconstruction of cached components
- Added scientific regression cases, executable documentation examples, installation matrices, and a parser benchmark
- Gated publishing on repository checks and release-version verification
- Updated the bundled specification to the unmodified ratified mzPAF 1.0.1 document

### New APIs

- Added explicit analyte-context resolution for peptide, internal, and precursor annotations
- Added structured parse errors, lazy record parsing, and batch results that retain individual failures
- Added versioned dictionary/JSON interchange with strict structural validation and preservation of resolved context
- Kept the existing as_dict() representation and context-free mass behavior compatible

## [1.1.1] (2026-08-15)

- Added machine-readable citation metadata and publication guidance
- Linked the existing Zenodo concept DOI for version-independent citation
- Added complete PyPI project metadata, community guidelines, and third-party notices
- Expanded CI coverage across supported Python versions and publication artifacts
- Prepared the GitHub release for archival in Zenodo

## [1.1.0] (2026-07-09)

- Fixed `ImmoniumIon.mass()`/`.composition` reading the wrong tacular lookup key (the internal by-fragment shift, `0`) instead of the immonium-specific `-CO` shift, which had been masked by an equivalent bug in tacular<1.1.0
- Fixed `ImmoniumIon.composition` silently dropping atom-removing modifications (e.g. Deamidated, Dehydrated): it now keeps net-negative element totals (stripping only exact zeros) so `.composition`/`.formula` stay consistent with `.mass()`
- Fixed `InternalFragment` ignoring its own `nterm_ion_type`/`cterm_ion_type` fields; non-default backbone cleavage types (e.g. `ax`, `cz`) now use the correct tacular mass/composition shift instead of silently reusing the default `by` (`0` shift) value, and setting only one of the two fields now raises a clear error
- Isotope annotations now contribute to `mass()`/`comp()` instead of being silently ignored: a generic isotope (`+i`, `+2i`, `-i`) applies the 13C−12C shift per mzPAF §4.6, so generic `+i` and explicit `+i13C` agree and `comp()` no longer raises for generic isotopes; average isotopes (`+iA`) raise consistently on both `mass()` and `comp()`
- Rejected element-specified isotopes with no nucleon count (e.g. `+iN`), per mzPAF §4.6
- Chemical-formula ions (`f{...}`) now compute `mass()`/`mz()` treating the formula as the fully charged species per mzPAF §4.4.9 (an electron-mass correction per charge, not an added proton; adducts label the charge but add no mass), matching the spec's worked `f{...}` m/z examples
- Fixed count-prefixed neutral-loss formulas such as `-2H2O` collapsing to a bare mass loss with the formula dropped; they now parse as the intended multiple-formula loss
- Neutral losses and adducts now accept isotope-labeled atoms mixed with plain atoms (e.g. `-H2[18O1]`, `[M+[2H2]]`, `[M+[15N1]H4]`), per mzPAF §4.5/§4.7
- Mass-error values now accept an optional leading `+` sign (matching the mzPAF reference regex)
- Named compounds (`_{...}`) now allow spaces, matching the spec's own `_{Urocanic Acid}` example
- `parse()`/`parse_multi()` now split comma-separated annotations using the mzPAF Appendix A greedy-match strategy instead of a naive comma split, so commas inside bracketed reference/named-compound labels no longer break splitting
- Fixed catastrophic backtracking (ReDoS) in the annotation regex on long single-element runs
- Bumped the `tacular` dependency floor to `>=1.1.0`, which fixes a systematic error in tacular's own internal-fragment-ion offset table
- Bumped the optional `peptacular` dependency floor to `>=3.1.2`, which fixes upstream composition/mass consistency bugs (e.g. atom-removing modifications silently dropped from `comp()`)

## [1.0.0] (2026-03-17)

- Added `to_mzpaf()` function to convert `peptacular` Fragment objects to `PafAnnotation`
- Added optional `peptacular` integration for sequence-aware calculations and fragment conversion
- Improved mass calculation for immonium ions with modifications
- Added caching for serialization and parsing results

## [0.1.0] (2026-01-14)

- First release on PyPI.
