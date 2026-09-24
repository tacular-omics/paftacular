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
   * - ``paftacular.AminoAcids``, ``paftacular.constants.AminoAcids``
     - ``tacular.AminoAcid``. ``ImmoniumIon.amino_acid`` is a ``tacular.AminoAcid``. Immonium
       ions still accept only the 20 standard codes (not B, J, O, U, X or Z).
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
   * - ``ann.dict_composition()`` (a method), ``component.dict_composition`` (a property)
     - ``{str(element): n for element, n in ann.comp().items()}``, and the same with
       ``component.composition``
   * - ``ann.as_dict()``, ``NeutralLoss.as_dict()``, ``IsotopeSpecification.as_dict()``,
       ``Adduct.as_dict()``
     - ``ann.to_dict()`` (versioned and reversible with ``PafAnnotation.from_dict``)
   * - ``to_mzpaf(frag, include_annotation=...)``
     - ``to_mzpaf(frag, include_sequence=...)``
   * - ``to_mzpaf(frag, mass_error_type=...)``
     - ``to_mzpaf(frag, mass_error_unit=...)``, the name ``PafAnnotation.from_components``
       and ``MassError.unit`` already use
   * - ``paftacular.constants.INTERNAL_SERIES_TO_DIFF``
     - private.
   * - ``INTERNAL_MASS_DIFFS`` (exported table of section 4.4.4 corrections)
     - private. ``PafAnnotation.make_internal(start, end, ion_type="bx")`` applies the
       correction (``m2:4+CO``). ``str(make_internal(2, 4, ion_type=key))`` shows it.
   * - ``composition_to_proforma_formula_string(comp, hill_order)``
     - ``composition_to_proforma_formula_string(comp)``. The ``hill_order`` argument had no
       effect and is gone.
   * - ``ReferenceIon._cache`` and the other class-level ``_cache`` dicts, instance
       interning in ``__new__``
     - removed. Constructors build fresh objects. The parser shares components for repeated
       substrings through bounded caches.
   * - MCP ``match_mz`` request ``"tolerance_unit": "Th"``
     - ``"tolerance_unit": "da"`` (tacular's ``ToleranceUnit``, ``"da"`` or ``"ppm"``). The
       value is still an absolute m/z difference. ``"Th"`` is rejected. The matched
       candidate field is still ``delta_th`` (observed minus theoretical m/z, in Th). MCP
       responses carry ``response_schema_version`` 2.
   * - MCP ``build_annotation`` ``charge`` / ``generate_fragments`` ``charges``: positive only
     - any nonzero integer, negative for negative mode (``-1`` writes ``^-1``). Zero is
       rejected. ``mz_th`` is mass / \|charge\|.
   * - MCP unknown top-level tool arguments silently ignored
     - rejected with an error, as unknown fields inside ``request`` already were.
   * - MCP ``get_capabilities`` fields ``ion_series``, ``resolvable_series``
     - ``ion_types``, ``resolvable_ion_types`` (plural lists, as in peptacular's MCP). Same
       values. No alias.

Keyword-only arguments
----------------------

Optional arguments must now be passed by name. Positional calls raise ``TypeError``.

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - 1.x
     - 2.0
   * - ``PafAnnotation(ion, None, False, (loss,))``
     - ``PafAnnotation(ion, neutral_losses=(loss,))``. Every field after ``ion_type`` is
       keyword-only.
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
   * - ``format_number(x, 5)``, ``validate_integer(x, "charge", 1)``
     - ``format_number(x, minimum_places=5)``, ``validate_integer(x, "charge", minimum=1)``
   * - ``to_mzpaf(frag, 0.9, 1.2)``
     - ``to_mzpaf(frag, confidence=0.9, mass_error=1.2)``

Errors
------

Every error caused by user input is now a ``PaftacularError``, a new base class that
subclasses ``ValueError``. ``PafParseError`` and ``PafUnknownReferenceError`` subclass it.
Code that catches ``ValueError`` keeps working. Errors that 1.x let escape from tacular or
peptacular (an unknown element in a formula, an invalid analyte in ``resolve()``) are now
wrapped in ``PaftacularError`` too. The mass or composition of a ``?`` or ``_{...}`` ion raises
``PafUnsupportedCalculationError``, a ``PaftacularError``, where 1.x raised
``NotImplementedError``.

Behaviour changes
-----------------

- **Negative charge.** ``y2{DE}^-2`` parses, a negative charge removes protons, and
  ``mz()`` divides by the absolute charge. ``serialize()`` writes ``^-2``.
  ``serialize(signed_charge=False)`` writes ``^2`` for readers that only accept the mzPAF
  1.0.1 positive form. Zero charge is still rejected. ``to_dict()`` output with a negative
  charge keeps ``schema_version`` 1, but paftacular 1.4.0 ``from_dict`` rejects it, so
  structured data is compatible forward only (1.x output loads in 2.0, not always the reverse).
- **Immonium adducts.** ``IK[M+K]`` is the K immonium ion with a K+ adduct. 1.x read
  ``M+K`` as the immonium modification, so ``serialize()`` output did not parse back.
  ``ImmoniumIon("K", modification="M+K")`` now raises ``PaftacularError``.
- **``to_mzpaf`` output** (with peptacular 5):

  - negative charges are kept (``b3{PEP}^-1``), where 1.x wrote ``^1`` or raised.
  - known neutral deltas and internal ion offsets use their canonical mzPAF names (``-NH3``,
    ``+HCOOH``, ``-HCONH2``). Other formula deltas are written in Hill order.
  - an immonium ion keeps a global fixed modification (``<[Oxidation]@P>`` gives
    ``IP[Oxidation]``) and a global isotope label (``<13C>`` gives ``IP+4i13C``). 1.x dropped
    both, so the label was off by the modification or label mass.
  - mass deltas of the same value are folded into one term and rounded to 6 decimals, like
    peptacular's own label. A delta that rounds to zero is left out.
  - several charge carriers are sorted alphabetically (``[M+H+Na]``), as mzPAF 4.7 asks.
  - a terminal modification on an immonium ion's residue becomes the immonium modification
    (``[Acetyl]-P`` gives ``IP[Acetyl]``), so the mass is kept. More than one modification
    raises ``PaftacularError``.

- **Exact offsets.** Every ion-type offset (all peptide series, immonium, internal,
  precursor) is summed from exact element masses, not tacular's 6-decimal constants. Offset
  masses move by up to ~4e-7 Da (y by 3.2e-7, immonium by 3.8e-7).
- **Listed modification masses.** A named modification (Unimod, PSI-MOD, RESID, XLMOD, GNO)
  in an embedded or resolved sequence counts at its listed database mass (Oxidation
  15.994915), the same rule as peptacular. Fragment and precursor ions, including neutral losses
  and isotope peaks, agree with peptacular 5 to 1e-9 Da. Composition is used only when there is no listed mass
  (formula modifications, glycans), and under a global isotope label (``<13C>``, ``<15N>``)
  in both packages. Unimod reference names (``r[Hex]``, ``-[Hex]``) use the listed 6-decimal
  mass too. mzPAF reference-list entries (``r[TMT6plex]``) keep their exact formula masses. ``comp()`` is
  unchanged, so the mass summed from ``comp()`` can differ from ``get_mass()`` by up to
  ~1e-6 Da for a named modification, because the listed mass is rounded.
- **Labile modifications.** A labile modification (``{Glycan:Hex}PEPTIDEK``) is lost on
  fragmentation, as ProForma defines it and peptacular computes it. Fragment ions no longer
  add its mass. Precursor ions keep it.
- **Charge carrier mass.** Monoisotopic charge is tacular's CODATA ``PROTON_MASS``, for the
  default charge and for an ``H`` carrier alike, so ``y2{DE}[M+H]`` equals ``y2{DE}``.
  An ``H`` carrier of the opposite sign (``[M+H]^-1``) is a hydride, an H atom plus an
  electron. Average charge is natural-abundance H less an electron, 1.16e-4 Da per charge
  heavier than 1.x.
- **Global isotope labels** (``<13C>``) in an embedded sequence replace their element in
  the ion offset and formula deltas too, like peptacular 5. ``a2{<13C>RY}`` has 14 13C,
  where 1.x counted 15. Atoms removed by a negative charge are labelled too: under ``<2H>``
  a ``^-1`` charge or ``[M-H]`` removes a deuteron. Mass-only deltas, isotope shifts, added
  adducts and the positive charge proton stay unlabelled. ``to_mzpaf`` counts an immonium
  label on the final ion (``IK-NH3+i15N``, ``IP+6i2H^-1``).
- **Global fixed modifications on side-chain ions.** A v ion loses one on its side-chain
  residue (``v3{<[Carbamidomethyl]@C>CFQ}`` is 349.151, not 406.172), and w and d ions
  raise ``PaftacularError``, as for explicit modifications.
- **ASCII digits only.** Formulas and the annotation grammar accept ASCII digits only.
  ``ChemicalFormula("H²O").get_mass()`` raises ``PaftacularError`` and ``y٢{DE}`` raises
  ``PafParseError``.
- **Formula tokens are checked at parse time.** ``IK[M+Methyl]``, ``y2{DE}[M+Methyl]`` and
  ``y2{DE}-Methyl`` raise ``PafParseError`` from ``parse()``. 1.x parsed them and failed in
  ``get_mass()``.
- **Caching.** Parsing the same substring twice shares one immutable component
  (``parse("y5-H2O").neutral_losses is parse("b3-H2O").neutral_losses``). Constructors no
  longer return interned objects, so compare components with ``==``, never ``is``.
