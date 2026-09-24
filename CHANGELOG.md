# Changelog

## [Unreleased]

### Fixed

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
  `ValueError` instead of `KeyError`.
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
