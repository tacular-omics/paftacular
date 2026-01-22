# paftacular

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)

A Python library for parsing and serializing **mzPAF** (Peak Annotation Format), a standardized format for annotating mass spectrometry fragment ions in peptide/proteomics analysis.

mzPAF is a specification from the [Proteomics Standards Initiative (PSI)](https://www.psidev.info/) that provides a compact, human-readable notation for describing fragment ion types, chemical modifications, charge states, mass errors, and confidence scores.

## Features

- **Comprehensive mzPAF support**: Handles all major annotation types from the specification
- **Parse complex annotations**: Fragment ions, neutral losses, isotopes, adducts, charge states, mass errors, and confidence scores
- **Mass calculations**: Supports both monoisotopic and average mass calculations
- **Composition analysis**: Full elemental composition tracking with ProForma-style formulas
- **Bidirectional**: Parse mzPAF strings to Python objects and serialize back
- **Type-safe**: Built with dataclasses, enums, and comprehensive type hints
- **Memory-efficient**: Uses slotted dataclasses for optimal performance
- **Export capabilities**: Convert annotations to dictionaries for JSON serialization
- **Multiple ion types**: Peptide fragments, precursors, immonium ions, reference compounds, chemical formulas, SMILES, and more

## Installation

```bash
pip install paftacular
```

Or for development:

```bash
git clone https://github.com/yourusername/paftacular.git
cd paftacular
just install
```

## Quick Start

```python
from paftacular import mzPAFParser

parser = mzPAFParser()

# Parse a simple peptide ion
ann = parser.parse_single("y5")
print(ann.ion_type.series)  # IonSeries.Y
print(ann.ion_type.position)  # 5

# Parse with modifications
ann = parser.parse_single("y5-H2O^2/1.2ppm*0.95")
print(ann.charge)              # 2
print(ann.mass_error.value)    # 1.2
print(ann.confidence)          # 0.95
```

## Usage Examples

### Basic Peptide Ions

```python
from paftacular import mzPAFParser

parser = mzPAFParser()

# Primary fragment ions (a, b, c, x, y, z series)
ann = parser.parse_single("b2")
ann = parser.parse_single("y3")

# With peptide sequence
ann = parser.parse_single("y3{PEP}")
print(ann.ion_type.sequence)  # "PEP"
print(ann.peptide_sequence)   # "PEP"

# Internal fragments
ann = parser.parse_single("m2:5{PEPTIDE}")
print(ann.ion_type.start_pos)  # 2
print(ann.ion_type.end_pos)    # 5
```

### Other Ion Types

```python
# Precursor ion
ann = parser.parse_single("p^2")  # Doubly-charged precursor

# Immonium ions
ann = parser.parse_single("IK")          # Lysine immonium
ann = parser.parse_single("IK[Acetyl]")  # Modified

# Reference compounds
ann = parser.parse_single("r[TMT126]")

# Chemical formula
ann = parser.parse_single("f{C6H12O6}")

# SMILES notation
ann = parser.parse_single("s{CN=C=O}")

# Unannotated/unknown peaks
ann = parser.parse_single("?")
ann = parser.parse_single("?42")  # With reporter number
```

### Modifications

```python
# Neutral losses (formula-based)
ann = parser.parse_single("y5-H2O")
ann = parser.parse_single("b3-NH3")
ann = parser.parse_single("b2-H2O-NH3")  # Multiple losses

# Neutral gains
ann = parser.parse_single("y5+NH3")

# Mass-based losses/gains
ann = parser.parse_single("y5-17.03")
ann = parser.parse_single("b2+22.98")

# Reference-based
ann = parser.parse_single("y5-[Adenine]")
```

### Isotopes

```python
# Generic isotope peak
ann = parser.parse_single("y5+i")

# Specific element
ann = parser.parse_single("y5+i13C")
ann = parser.parse_single("b2+2i13C")  # Multiple isotopes

# Average isotopomer
ann = parser.parse_single("y5+iA")
```

### Adducts

```python
# Single adduct
ann = parser.parse_single("y5[M+Na]")

# Multiple adducts
ann = parser.parse_single("y5[M+2H+Na]")
ann = parser.parse_single("p[M+NH4]")
```

### Charge States

```python
# Singly charged (default)
ann = parser.parse_single("y5")
print(ann.charge)  # 1

# Doubly charged
ann = parser.parse_single("y5^2")
print(ann.charge)  # 2

# Triply charged
ann = parser.parse_single("b3^3")
print(ann.charge)  # 3
```

### Mass Errors and Confidence

```python
# Mass error in ppm
ann = parser.parse_single("y5/1.2ppm")
print(ann.mass_error.value)  # 1.2
print(ann.mass_error.unit)   # "ppm"

# Mass error in daltons
ann = parser.parse_single("y5/-0.003da")

# Confidence score (0.0 to 1.0)
ann = parser.parse_single("y5*0.95")
print(ann.confidence)  # 0.95

# Combined
ann = parser.parse_single("y5/1.2ppm*0.95")
```

### Complex Annotations

```python
# Everything at once
ann = parser.parse_single("&2@y5-H2O+i13C[M+H+Na]^2/-0.55ppm*0.85")
print(ann.is_auxiliary)        # True (& flag)
print(ann.analyte_reference)   # 2 (@ reference)
print(ann.ion_type.series)     # IonSeries.Y
print(ann.ion_type.position)   # 5
print(ann.neutral_mods)        # Tuple with H2O loss
print(ann.isotopes)            # Tuple with 13C isotope
print(ann.adducts)             # Tuple with H and Na
print(ann.charge)              # 2
print(ann.mass_error.value)    # -0.55
print(ann.confidence)          # 0.85
```

### Mass Calculations

```python
# Get calculated masses
ann = parser.parse_single("y5")
print(ann.monoisotopic_mass)   # Monoisotopic mass
print(ann.average_mass)        # Average mass
print(ann.mass)                # Defaults to monoisotopic

# With modifications
ann = parser.parse_single("y5-H2O")
print(ann.monoisotopic_mass)   # Mass adjusted for water loss

# Get composition
ann = parser.parse_single("f{C6H12O6}")
comp = ann.composition  # Counter with ElementInfo keys
print(ann.formula)              # "C6H12O6"
print(ann.proforma_formula)     # ProForma-style formula
```

### Parsing Multiple Annotations

```python
# Parse comma-separated list
annotations = parser.parse("b2, y3-H2O, p^2")
for ann in annotations:
    print(ann.serialize())
# Output:
# b2
# y3-H2O
# p^2
```

### Serialization (Round-trip)

```python
# Parse and serialize back
original = "y5-H2O^2/1.2ppm*0.95"
ann = parser.parse_single(original)
serialized = ann.serialize()
print(serialized)  # "y5-H2O^2/1.2ppm*0.95"
```

### Export to Dictionary

```python
# Export annotation as dictionary (e.g., for JSON)
ann = parser.parse_single("y5-H2O^2/1.2ppm*0.95")
data = ann.as_dict()
print(data)
# {
#   'ion': {'ion_type': 'PeptideIon', 'series': IonSeries.Y, 'position': 5, ...},
#   'neutral_losses': [{'count': -1, 'base_formula': 'H2O', ...}],
#   'charge': 2,
#   'mass_error': {'value': 1.2, 'unit': 'ppm'},
#   'confidence': 0.95,
#   ...
# }
```

## API Reference

### `mzPAFParser`

Main parser class for mzPAF annotations.

#### Methods

- **`parse_single(annotation: str) -> PafAnnotation`**  
  Parse a single mzPAF annotation string.

- **`parse(annotations: str) -> list[PafAnnotation]`**  
  Parse comma-separated mzPAF annotation strings.

### `PafAnnotation`

Dataclass representing a complete fragment annotation.

#### Properties

- `ion_type`: The ion type (peptide, precursor, immonium, etc.)
- `analyte_reference`: Optional integer reference to specific analyte
- `is_auxiliary`: Boolean flag for auxiliary peaks
- `neutral_mods`: Tuple of neutral loss/gain modifications
- `isotopes`: Tuple of isotope specifications
- `adducts`: Tuple of adduct modifications
- `charge`: Charge state (default=1)
- `mass_error`: Mass error with unit (da or ppm)
- `confidence`: Confidence score (0.0-1.0)

#### Methods

- **`mass`**: Calculate ion mass (defaults to monoisotopic)
- **`monoisotopic_mass`**, **`average_mass`**: Get specific mass types
- **`composition`**: Get elemental composition as `Counter[ElementInfo]`
- **`formula`**, **`proforma_formula`**: Get chemical formula strings
- **`peptide_sequence`**: Get peptide sequence if applicable
- **`serialize() -> str`**: Convert back to mzPAF string
- **`as_dict() -> dict`**: Convert to dictionary representation for JSON serialization

### Ion Types

- `PrimaryFragmentIon`: a, b, c, x, y, z series
- `InternalFragmentIon`: Internal fragments (m notation)
- `ImmoniumIon`: Immonium ions (I notation)
- `PrecursorIon`: Precursor/parent ions (p notation)
- `ReferenceIon`: Reference compounds (r notation)
- `FormulaIon`: Formula-based ions (f notation)
- `NamedCompoundIon`: Named compounds (_ notation)
- `SmilesIon`: SMILES notation (s notation)
- `UnannotatedIon`: Unannotated peaks (? notation)

## Requirements

- Python ≥3.12
- [tacular](https://github.com/yourusername/tacular) - Element and reference molecule lookups

### Optional Dependencies

- `pysmiles>=2.0.1` - For SMILES notation support

## Development

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/paftacular.git
cd paftacular

# Install with dev dependencies
just install
```

### Commands

```bash
# Run linter
just lint

# Format code
just format

# Type check
just check

# Run tests
just test
```

## mzPAF Format Reference

The mzPAF format uses a compact notation for annotations:

```
[&][analyte@]ion_type[modifications][^charge][/mass_error][*confidence]
```

**Examples:**
- `y5` - Y-series ion at position 5
- `b2{PEP}` - B-series ion with sequence PEP
- `y5-H2O` - Y5 with water loss
- `y5^2` - Doubly-charged Y5
- `y5/1.2ppm` - Y5 with 1.2 ppm mass error
- `y5*0.95` - Y5 with 95% confidence
- `&2@y5-H2O^2/1.2ppm*0.95` - Full annotation

For more details, see the [PSI mzPAF specification](https://www.psidev.info/).

## License

See project metadata for license information.

## Author

Patrick Garrett (pgarrett@scripps.edu)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
