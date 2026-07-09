# History

## 0.1.0 (2026-01-14)

- First release on PyPI.

## 1.0.0 (2026-03-17)

- Added `to_mzpaf()` function to convert `peptacular` Fragment objects to `PafAnnotation`
- `peptacular` is now a required dependency (previously optional)
- Added Zenodo DOI for citation
- Improved mass calculation for immonium ions with modifications
- Added caching for serialization and parsing results

## 1.1.0 (2026-07-09)

- Fixed `ImmoniumIon.mass()`/`.composition` reading the wrong tacular lookup key (the internal by-fragment shift, `0`) instead of the immonium-specific `-CO` shift, which had been masked by an equivalent bug in tacular<1.1.0
- Fixed `ImmoniumIon.composition` silently dropping atom-removing modifications (e.g. Deamidated) due to `Counter` arithmetic discarding non-positive running totals mid-accumulation
- Fixed `InternalFragment` ignoring its own `nterm_ion_type`/`cterm_ion_type` fields; non-default backbone cleavage types (e.g. `ax`, `cz`) now use the correct tacular mass/composition shift instead of silently reusing the default `by` (`0` shift) value, and setting only one of the two fields now raises a clear error
- Bumped the `tacular` dependency floor to `>=1.1.0`, which fixes a systematic error in tacular's own internal-fragment-ion offset table
- Bumped the optional `peptacular` dependency floor to `>=3.1.2`, which fixes upstream composition/mass consistency bugs (e.g. atom-removing modifications silently dropped from `comp()`)
