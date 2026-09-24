Migrating from 1.x to 2.0
=========================

paftacular 2.0 is a breaking release. It requires ``tacular>=2.0,<3`` and, for the
``peptacular``, ``mcp`` and ``all`` extras, ``peptacular>=5.0,<6``. This page lists every
public name that was renamed or removed, and every behaviour change a caller can notice.

Renamed and removed names
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - 1.x
     - 2.0
   * - ``parse(s)`` returning one annotation or a list
     - ``parse(s)`` always returns exactly one ``PafAnnotation`` and raises ``PafParseError``
       otherwise. Use ``parse_multi(s)`` for comma-separated text and for ``""``.
   * - ``parse_single(s)``
     - ``parse(s)``
   * - ``parse_batch(records)``
     - ``list(iter_parse(records))``
   * - ``mzPAFParser``, ``mzPAFParser().parse(s)``, ``.parse_multi(s)``
     - removed. Use the module functions ``parse`` and ``parse_multi``.
   * - ``paftacular.parser.MZ_PAF_PARSER``
     - removed
   * - ``ann.mass(monoisotopic, calculate_sequence)``
     - ``ann.get_mass(*, monoisotopic=True, calculate_sequence=True)``
   * - ``component.mass(monoisotopic)`` on every ion component and modifier
     - ``component.get_mass(*, monoisotopic=True)``
   * - ``ann.mz(monoisotopic, calculate_sequence)``
     - ``ann.mz(*, monoisotopic=True)``. It has no ``calculate_sequence`` and raises
       ``PaftacularError`` for a peptide, internal or precursor ion without a sequence.
       Call ``resolve()`` or embed the sequence first. For the old offset value use
       ``ann.get_mass(calculate_sequence=False) / abs(ann.charge)``.
   * - ``ann.dict_composition()``, ``component.dict_composition()``
     - ``{str(element): n for element, n in ann.comp().items()}``
   * - ``ann.as_dict()``, ``NeutralLoss.as_dict()``, ``IsotopeSpecification.as_dict()``,
       ``Adduct.as_dict()``
     - ``ann.to_dict()`` (versioned and reversible with ``PafAnnotation.from_dict``)
   * - ``to_mzpaf(frag, include_annotation=...)``
     - ``to_mzpaf(frag, include_sequence=...)``
   * - ``paftacular.constants.INTERNAL_SERIES_TO_DIFF``
     - private. Use ``INTERNAL_MASS_DIFFS[(nterm, cterm)]``.
   * - ``INTERNAL_MASS_DIFFS`` as a ``dict``
     - a read-only mapping. Copy it with ``dict(INTERNAL_MASS_DIFFS)`` to change it.
   * - ``composition_to_proforma_formula_string(comp, hill_order)``
     - ``composition_to_proforma_formula_string(comp)``. The ``hill_order`` argument had no
       effect and is gone.
   * - ``ReferenceIon._cache`` and the other class-level ``_cache`` dicts, instance
       interning in ``__new__``
     - removed. Constructors build fresh objects. The parser shares components for repeated
       substrings through bounded caches.

Keyword-only arguments
----------------------

Optional arguments must now be passed by name. Positional calls raise ``TypeError``.

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - 1.x
     - 2.0
   * - ``PeptideIon("y", 3, "PEP")``
     - ``PeptideIon("y", 3, sequence="PEP")``
   * - ``InternalFragment(2, 4, "EPT", "a", "x")``
     - ``InternalFragment(2, 4, sequence="EPT", nterm_ion_type="a", cterm_ion_type="x")``
   * - ``ImmoniumIon("M", "Oxidation")``
     - ``ImmoniumIon("M", modification="Oxidation")``
   * - ``UnknownIon(42)``
     - ``UnknownIon(label=42)``
   * - ``NeutralLoss(-1, "H2O")``
     - ``NeutralLoss(-1, base_formula="H2O")`` (also ``base_mass``, ``base_reference``)
   * - ``IsotopeSpecification(1, "13C")``
     - ``IsotopeSpecification(1, element="13C")`` (also ``is_average``)
   * - ``MassError(1.2, "ppm")``
     - ``MassError(1.2, unit="ppm")``
   * - ``make_peptide("y", 3, "PEP")``
     - ``make_peptide("y", 3, sequence="PEP")``
   * - ``make_internal(2, 4, "bx", "EPT")``
     - ``make_internal(2, 4, ion_type="bx", sequence="EPT")``
   * - ``make_immonium("M", "Oxidation")``
     - ``make_immonium("M", modification="Oxidation")``
   * - ``make_unknown(42)``
     - ``make_unknown(label=42)``
   * - ``ann.comp(False)``, ``ann.formula(False)``, ``ann.proforma_formula(False)``
     - ``calculate_sequence=False``
   * - ``ann.serialize(False)``
     - ``ann.serialize(include_sequence=False)``
   * - ``loss.serialize("mass", True)``
     - ``loss.serialize(loss_type="mass", monoisotopic=True)``
   * - ``format_number(x, 5)``, ``validate_integer(x, 1)``
     - ``format_number(x, minimum_places=5)``, ``validate_integer(x, minimum=1)``
   * - ``to_mzpaf(frag, 0.9, 1.2)``
     - ``to_mzpaf(frag, confidence=0.9, mass_error=1.2)``

Errors
------

Every error caused by user input is now a ``PaftacularError``, a new base class that
subclasses ``ValueError``. ``PafParseError`` and ``PafUnknownReferenceError`` subclass it.
Code that catches ``ValueError`` keeps working. Errors that 1.x let escape from tacular or
peptacular (an unknown element in a formula, an invalid analyte in ``resolve()``) are now
wrapped in ``PaftacularError`` too. ``NotImplementedError`` is still raised for the mass of
``?`` and ``_{...}`` ions.

Behaviour changes
-----------------

- **Negative charge.** ``y2{DE}^-2`` parses, a negative charge removes protons, and
  ``mz()`` divides by the absolute charge. ``serialize()`` writes ``^-2``.
  ``serialize(signed_charge=False)`` writes ``^2`` for readers that only accept the mzPAF
  1.0.1 positive form. Zero charge is still rejected.
- **``to_mzpaf`` output** (with peptacular 5):

  - negative charges are kept (``b3{PEP}^-1``), where 1.x wrote ``^1`` or raised.
  - known neutral deltas use their canonical mzPAF names (``-NH3``, ``+HCOOH``,
    ``-HCONH2``). Other formula deltas are written in Hill order.
  - mass deltas of the same value are folded into one term and rounded to 6 decimals, like
    peptacular's own label. A delta that rounds to zero is left out.
  - several charge carriers are sorted alphabetically (``[M+H+Na]``), as mzPAF 4.7 asks.
  - a terminal modification on an immonium ion's residue becomes the immonium modification
    (``[Acetyl]-P`` gives ``IP[Acetyl]``), so the mass is kept. More than one modification
    raises ``PaftacularError``.

- **Caching.** Parsing the same substring twice shares one immutable component
  (``parse("y5-H2O").neutral_losses is parse("b3-H2O").neutral_losses``). Constructors no
  longer return interned objects, so compare components with ``==``, never ``is``.
