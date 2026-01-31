# paftacular

<div align="center">
  <img src="paftacular_logo.png" alt="Paftacular Logo" width="400" style="margin: 50px;"/>

  # paftacular

  A Python library for parsing and serializing **mzPAF** (Peak Annotation Format), a standardized format for annotating mass spectrometry fragment ions in peptide/proteomics analysis. mzPAF is a specification from the [Proteomics Standards Initiative (PSI)](https://www.psidev.info/) that provides a compact, human-readable notation for describing fragment ion types, chemical modifications, charge states, mass errors, and confidence scores.

    
  [![Python package](https://github.com/tacular-omics/ppaftacular/actions/workflows/python-package.yml/badge.svg)](https://github.com/tacular-omics/paftacular/actions/workflows/python-package.yml)
  [![codecov](https://codecov.io/github/tacular-omics/paftacular/graph/badge.svg?token=1CTVZVFXF7)](https://codecov.io/github/tacular-omics/paftacular)
  [![PyPI version](https://badge.fury.io/py/paftacular.svg)](https://badge.fury.io/py/paftacular)
  [![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
  [![License: MIT](https://img.shields.io/badge/License-MIT-g.svg)](https://opensource.org/licenses/MIT)
  
</div>

## Features

- Parse and serialize mzPAF annotations with full specification support
- Handle fragment ions, neutral losses, isotopes, adducts, charge states, mass errors, and confidence scores
- Calculate monoisotopic and average masses with elemental composition tracking
- Type-safe with comprehensive type hints and dataclasses

## Installation

```bash
   pip install paftacular
   pip install paftacular[sequence] # for peptide sequence support
   pip install paftacular[smiles] # with smiles support
   pip install paftacular[all] # with all optional dependencies
```

## Quick Start


There are 3 parsing methods available:
* ``parse``: Parses a single or multiple comma-separated mzPAF annotations. Returns a single ``PafAnnotation`` or a list of them.
* ``parse_multi``: Parses multiple comma-separated mzPAF annotations. Always returns a list of ``PafAnnotation``.
* ``parse_single``: Parses a single mzPAF annotation. Returns a single ``PafAnnotation``. Raises ValueError if multiple annotations are provided.


```python
import paftacular as pft

# Parse a simple peptide ion
ann = pft.parse("y5")
print(ann.ion_type.series)  # IonSeries.Y
print(ann.ion_type.position)  # 5

# Calculate masses
print(ann.monoisotopic_mass)   # Calculated mass
print(ann.serialize())         # Round-trip back to string

# Parse multiple ions
anns = pft.parse("y5-H2O^2/1.2ppm*0.95,b3^2")
for ann in anns:
  print(ann.charge)
  print(ann.mass_error.value)
  print(ann.confidence)
```

## Documentation

Full documentation is available at [Read the Docs](https://paftacular.readthedocs.io/).

## mzPAF Format

The mzPAF format uses compact notation:

```
[&][analyte@]ion_type[modifications][^charge][/mass_error][*confidence]
```

Examples: `y5`, `b2{PEP}`, `y5-H2O^2`, `y5/1.2ppm*0.95`

See the [PSI mzPAF specification](https://www.psidev.info/) for full details.

## License

MIT

## Contributing

Contributions welcome! Please submit a Pull Request.

**Author:** Patrick Garrett (pgarrett@scripps.edu)
