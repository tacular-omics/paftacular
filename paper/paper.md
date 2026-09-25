---
title: 'paftacular: A Python implementation of the HUPO-PSI mzPAF peak annotation format'
tags:
  - Python
  - Proteomics
  - Mass Spectrometry
  - mzPAF
  - Peak Annotation
  - HUPO-PSI
authors:
  - name: Patrick T. Garrett
    orcid: 0000-0002-8434-9693
    affiliation: 1
  - given-names: John R.
    surname: Yates
    suffix: III
    orcid: 0000-0001-5267-1672
    corresponding: true
    affiliation: 1
affiliations:
  - name: The Scripps Research Institute, United States
    index: 1
date: 24 September 2026
bibliography: paper.bib
---

# Summary

Tandem mass spectrometry identifies peptides by fragmenting them and matching the observed
fragment peaks to predicted ions. Each matched peak carries an explanation: which ion series
and position it belongs to, its charge, any neutral loss or isotope peak, and how far the
observed m/z lies from the theoretical value. The HUPO Proteomics Standards Initiative (PSI)
standardized these explanations as mzPAF, the Peak Annotation Format [@klein-2024;
@mzpaf-spec]. An mzPAF annotation such as `y4-H2O^2/1.2ppm*0.85` is a compact string that
spectral libraries, search engines and viewers can exchange.

**paftacular** is a Python library that parses mzPAF text into typed, immutable objects,
serializes those objects back to text, and computes the mass, m/z and elemental composition
of the annotated ions. It covers the ion types of mzPAF 1.0.1, including peptide series,
internal fragments, immonium, precursor, reference, formula, SMILES, named-compound and
unknown ions, together with neutral losses and gains, isotopes, adducts, charge, analyte
references, mass errors and confidence values. Peptide ions can be resolved against a full
ProForma [@leduc-2022] sequence to give complete masses. An optional Model Context Protocol
(MCP) server exposes the same operations to AI clients.

# Statement of need

mzPAF strings pack a great deal of chemistry into a few characters, and the grammar has
many corner cases: bracketed reference-molecule names that contain `-`, `_` or parentheses,
isotope shifts for named elements, adducts with explicit electron accounting, and embedded
ProForma sequences. Ad hoc string splitting handles common cases and fails silently on the
rest. Software that reads or writes spectral libraries, annotates spectra or validates
annotations needs a parser that either returns a complete structured object or reports
exactly where the text is invalid, and a serializer whose output parses back to the same
object.

Structure alone is not enough for most uses. To check an annotation against a spectrum, a
program must compute the theoretical m/z of the ion it describes, and that value must agree
with the fragment masses produced by the sequence tools that generated the annotation in the
first place. paftacular targets developers of spectral-library, spectrum-annotation and
proteomics-pipeline software who need mzPAF parsing, serialization and mass calculation in
one tested package.

# State of the field

The PSI maintains a reference Python implementation, published on PyPI as `mzpaf`
[@mzpaf-python]. It parses annotations with a regular expression into typed objects,
serializes them, and converts them to a JSON data model. The HUPO-PSI `mzspeclib-py`
library for the companion mzSpecLib spectral-library format [@klein-2024] imports its
annotation classes from this package [@mzspeclib-py]. In the latest PyPI release
(0.2.0b0), the reference implementation computes the mass of a neutral loss but not of an
annotated ion, and it rejects bracketed losses such as `b2-[HexNAc(2)]` and
`p-[TMT126-ETD]`, which the mzPAF 1.0.1 grammar allows. We checked both behaviours by
running that release.

spectrum_utils [@bittremieux-2023] generates theoretical fragments for spectrum plotting and
writes their labels in a notation derived from a draft of the PSI peak-annotation
specification. Its `FragmentAnnotation` class builds and prints labels but does not parse
annotation text, and it requires a strictly positive charge. Pyteomics
[@goloborodko-2013; @levitsky-2019] provides the mass calculations that both of these
packages use, but its current release (5.0.1) contains no mzPAF module. In Rust, the
mzcore project [@Schulte_mzcore] reads and writes mzPAF in its `mzannotate` crate as part of
a spectrum annotator [@schulte-2025]; its Python bindings return the mzPAF label of a
generated fragment.

paftacular differs from these in three ways. It computes monoisotopic and average mass, m/z
and composition for every ion type that carries chemistry, with or without a resolved
peptide sequence. Its calculations are tested against an independent reference and against
the peptacular fragment generator. And it reports invalid text with a zero-based character
position and annotation index, so that a tool can point a user to the problem.

# Software design

**Parsing and serialization.** `parse()` accepts exactly one annotation and `parse_multi()`
a comma-separated list; `iter_parse()` processes many records lazily and keeps each failure
with its index instead of stopping the batch. Parsed components are frozen dataclasses, and
repeated ion and modifier text shares components through bounded, thread-safe caches.
`serialize()` writes canonical mzPAF, and `to_dict()`/`from_dict()` provide a versioned JSON
interchange with strict structural validation. The test suite parses and round-trips 100
annotation strings printed in the mzPAF 1.0.1 specification, and Hypothesis
[@hypothesis] property tests check that generated annotations survive the serialize and
parse round trip.

**Mass and m/z.** Composition and masses use element and modification data from tacular
[@garrett-tacular], which bundles Unimod and the other ontologies, so no data are downloaded
at runtime. Ion-type offsets are summed from exact element masses. With the optional
peptacular [@garrett-peptacular] extra, `resolve()` selects the fragment's residues from a
full ProForma analyte, and `to_mzpaf()` converts peptacular fragments into `PafAnnotation`
objects. Named modifications count at their listed database mass, the same rule peptacular
uses, so the two packages agree: `tests/test_peptacular_agreement.py` compares m/z and mass
for a, b, c, x, y and z ions at charges 1 to 3, precursor ions, average masses, neutral
losses and isotope peaks across nine modified peptides, including labile glycans, fixed
modifications and PSI-MOD accessions, and requires agreement within 1e-9 Da (525 test
cases). Independently, 532 m/z values are checked, within 1e-6 to 5e-6, against a frozen fixture that
`tests/reference/generate_reference.py` builds from the specification's formulas with
Pyteomics atomic data, without importing paftacular, tacular or peptacular.

**Negative charges.** mzPAF 1.0.1 writes a negative-mode charge without the minus sign and
takes the polarity from the spectrum identification. paftacular also accepts a signed charge
such as `y2{DE}^-2`, computes it as a deprotonated ion, and by default writes the sign back
so that a stored annotation keeps its polarity outside its spectrum. `serialize(signed_charge=False)`
writes the unsigned form that the specification prescribes.

**Validation and errors.** Every error caused by user input is a `PaftacularError`, a
subclass of `ValueError`. `PafParseError` carries the input text, a zero-based character
position, the annotation index and a reason. Calculations that the annotation cannot define,
such as the mass of an unknown `?` ion or a named compound, raise
`PafUnsupportedCalculationError`, and unknown reference-molecule names raise
`PafUnknownReferenceError`. Adduct and neutral-loss formulas with a token that is not an
element are rejected when parsed, not later during mass calculation, and errors raised by tacular or peptacular are re-raised
as `PaftacularError`.

**MCP server.** The optional `paftacular[mcp]` extra installs a local stdio server
(`paftacular-mcp`) for the Model Context Protocol [@mcp-spec]. It provides nine read-only
tools for parsing, constructing, resolving, serializing and calculating annotations,
generating fragment series and matching candidates to an observed m/z within a ppm or Da
tolerance, plus four reference resources and two prompts. The tools call the public Python
API, and the base installation does not import the MCP SDK.

The package requires Python 3.12 or later and tacular; peptacular, SMILES support
(pysmiles) and MCP are optional extras. Continuous integration runs the tests on Python
3.12 to 3.14 on Linux, and on Python 3.13 on macOS and Windows. The full suite (1940 passed, 4 expected
failures) was run for this paper against the development versions of paftacular, peptacular
and tacular.

# Example usage

The script below parses an annotation, resolves it against a modified peptide, compares the
result with peptacular's own fragment, writes a negative-mode ion both ways, and shows a
parse error. The output shown was produced by running it.

```python
import paftacular as pft
import peptacular as pt

ann = pft.parse("y4-H2O^2/1.2ppm*0.85")
print(ann.ion_type.series, ann.ion_type.position, ann.charge, ann.confidence)
print(ann.serialize())

ion = ann.resolve("PEM[Oxidation]TIDEK")
print(ion.sequence, round(ion.mz(), 6))
frag = pt.parse("PEM[Oxidation]TIDEK").frag(
    ion_type="y", position=4, charge=2, deltas={"H2O": 1})
print(round(frag.mz, 6))

neg = pft.parse("y2{DE}^-2")
print(neg.serialize(), neg.serialize(signed_charge=False), round(neg.mz(), 6))

try:
    pft.parse_multi("b3,y2x")
except pft.PafParseError as err:
    print(err.annotation_index, err.position, err.reason)
```

```text
y 4 2 0.85
y4-H2O^2/1.2ppm*0.85
IDEK 243.631558
243.631558
y2{DE}^-2 y2{DE}^2 130.032774
1 5 Unexpected or missing annotation content
```

# Research impact statement

Within the authors' software, the spectrum-processing package spxtacular depends on
paftacular to represent matched peaks as mzPAF annotations, and peptacular's fragment
generator emits mzPAF strings that paftacular parses directly.

<!-- TODO(author): add any use of paftacular outside the tacular-omics packages (papers,
pipelines, spectral libraries, in-house use). None is recorded in the repository, so none
is claimed here. -->

# AI usage disclosure

<!-- TODO(author): AI-usage disclosure required by JOSS. Describe which AI tools were used,
for what (code, tests, documentation, this paper) and how their output was reviewed. -->

# Acknowledgements

<!-- TODO(author): acknowledgements and funding. -->

# References
