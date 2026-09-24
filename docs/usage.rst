Usage
=====

.. contents:: Table of Contents
   :local:
   :depth: 2

.. testsetup:: *

   import paftacular as pft

This guide covers the full mzPAF grammar. For a five-minute introduction, see
:doc:`quickstart`.

Basic Peptide Ions
------------------

.. testcode::

   import paftacular as pft

   # Primary fragment ions (a, b, c, x, y, z series)
   ann = pft.parse("b2")
   print(ann.ion_type.series)
   
   ann = pft.parse("y3")
   print(ann.ion_type.series)

.. testoutput::

   b
   y

.. testcode::

   # With peptide sequence
   ann = pft.parse("y3{PEP}")
   print(ann.ion_type.sequence)

.. testoutput::

   PEP

.. testcode::

   # Internal fragments
   ann = pft.parse("m2:5{PEPTIDE}")
   print(ann.ion_type.start_position)
   print(ann.ion_type.end_position)

.. testoutput::

   2
   5

Other Ion Types
---------------

.. code-block:: python

   import paftacular as pft

   # Precursor ion
   ann = pft.parse("p^2")  # Doubly-charged precursor

   # Immonium ions
   ann = pft.parse("IK")          # Lysine immonium
   ann = pft.parse("IK[Acetyl]")  # Modified

   # Reference compounds
   ann = pft.parse("r[TMT126]")

   # Chemical formula
   ann = pft.parse("f{C6H12O6}")

   # SMILES notation
   ann = pft.parse("s{CN=C=O}")

   # Unannotated/unknown peaks
   ann = pft.parse("?")
   ann = pft.parse("?42")  # With reporter number

Modifications
-------------

.. code-block:: python

   import paftacular as pft

   # Neutral losses (formula-based)
   ann = pft.parse("y5-H2O")
   ann = pft.parse("b3-NH3")
   ann = pft.parse("b2-H2O-NH3")  # Multiple losses

   # Neutral gains
   ann = pft.parse("y5+NH3")

   # Mass-based losses/gains
   ann = pft.parse("y5-17.03")
   ann = pft.parse("b2+22.98")

   # Reference-based
   ann = pft.parse("y5-[Adenine]")

Isotopes
--------

.. code-block:: python

   import paftacular as pft

   # Generic isotope peak
   ann = pft.parse("y5+i")

   # Specific element
   ann = pft.parse("y5+i13C")
   ann = pft.parse("b2+2i13C")  # Multiple isotopes

   # Average isotopomer
   ann = pft.parse("y5+iA")

Adducts
-------

.. code-block:: python

   import paftacular as pft

   # Single adduct
   ann = pft.parse("y5[M+Na]")

   # Multiple adducts
   ann = pft.parse("y5[M+2H+Na]^3")
   ann = pft.parse("p[M+NH4]")

Charge States
-------------

.. code-block:: python

   import paftacular as pft

   # Singly charged (default)
   ann = pft.parse("y5")
   print(ann.charge)  # 1

   # Doubly charged
   ann = pft.parse("y5^2")
   print(ann.charge)  # 2

   # Triply charged
   ann = pft.parse("b3^3")
   print(ann.charge)  # 3

   # Negative mode: written as ^-2 so the round trip keeps the sign
   ann = pft.parse("y5{PEPTI}^-2")
   print(ann.charge)       # -2
   print(ann.serialize())  # y5{PEPTI}^-2
   print(ann.serialize(signed_charge=False))  # y5{PEPTI}^2 (mzPAF 1.0.1 section 4.8)

Mass Errors and Confidence
---------------------------

.. code-block:: python

   import paftacular as pft

   # Mass error in ppm
   ann = pft.parse("y5/1.2ppm")
   print(ann.mass_error.value)  # 1.2
   print(ann.mass_error.unit)   # "ppm"

   # Mass error in daltons
   ann = pft.parse("y5/-0.003")

   # Confidence score (0.0 to 1.0)
   ann = pft.parse("y5*0.95")
   print(ann.confidence)  # 0.95

   # Combined
   ann = pft.parse("y5/1.2ppm*0.95")

Complex Annotations
-------------------

.. code-block:: python

   import paftacular as pft

   # Everything at once
   ann = pft.parse("&2@y5-H2O+i13C[M+H+Na]^2/-0.55ppm*0.85")
   print(ann.is_auxiliary)        # True (& flag)
   print(ann.analyte_reference)   # 2 (@ reference)
   print(ann.ion_type.series)     # y
   print(ann.ion_type.position)   # 5
   print(ann.neutral_losses)      # Tuple with H2O loss
   print(ann.isotopes)            # Tuple with 13C isotope
   print(ann.adducts)             # Tuple with H and Na
   print(ann.charge)              # 2
   print(ann.mass_error.value)    # -0.55
   print(ann.confidence)          # 0.85

Mass Calculations
-----------------

.. code-block:: python

   import paftacular as pft

   # Without sequence context get_mass() returns the ion offset, not a complete fragment mass
   ann = pft.parse("y5")
   print(ann.get_mass(monoisotopic=True))    # Monoisotopic mass
   print(ann.get_mass(monoisotopic=False))   # Average mass

   # mz() needs a sequence for peptide, internal and precursor ions
   print(pft.parse("y2{DE}").mz())

   # With modifications
   ann = pft.parse("y5-H2O")
   print(ann.get_mass(monoisotopic=True))    # Mass adjusted for water loss

   # Get composition
   ann = pft.parse("f{C6H12O6}")
   comp = ann.comp(calculate_sequence=True)  # Counter with ElementInfo keys
   print(ann.formula())              # "C6H12O6" (formula ions already specify all atoms)
   print(ann.proforma_formula())     # ProForma-style formula

Parsing Multiple Annotations
-----------------------------

.. testcode::

   import paftacular as pft

   # parse() takes exactly one annotation. Use parse_multi() for a comma-separated list
   annotations = pft.parse_multi("b2, y3-H2O, p^2")
   for ann in annotations:
       print(ann.serialize())

.. testoutput::

   b2
   y3-H2O
   p^2

Creating Annotations Programmatically
--------------------------------------

Instead of parsing strings, you can create annotations directly using factory methods.

Basic Factory Methods
~~~~~~~~~~~~~~~~~~~~~

.. testcode::

   from paftacular import PafAnnotation

   # Create a precursor ion
   ann = PafAnnotation.make_precursor()
   print(ann.serialize())

.. testoutput::

   p

.. testcode::

   # Create a peptide fragment ion
   ann = PafAnnotation.make_peptide("y", 5)
   print(ann.serialize())

.. testoutput::

   y5

.. testcode::

   # Create an internal fragment
   ann = PafAnnotation.make_internal(start_position=2, end_position=5)
   print(ann.serialize())

.. testoutput::

   m2:5

.. testcode::

   # Create an immonium ion
   ann = PafAnnotation.make_immonium("K")
   print(ann.serialize())

.. testoutput::

   IK

.. testcode::

   # Create a reference ion
   ann = PafAnnotation.make_reference("TMT126")
   print(ann.serialize())

.. testoutput::

   r[TMT126]

Adding Modifications
~~~~~~~~~~~~~~~~~~~~

All factory methods accept common parameters for modifications.

.. testcode::

   from paftacular import PafAnnotation

   # Add neutral losses
   ann = PafAnnotation.make_peptide("y", 5, neutral_losses=["-H2O", "-NH3"])
   print(ann.serialize())

.. testoutput::

   y5-H2O-NH3

.. testcode::

   # Add isotopes
   ann = PafAnnotation.make_peptide("b", 3, isotopes=["+i13C", "+i15N"])
   print(ann.serialize())

.. testoutput::

   b3+i13C+i15N

.. testcode::

   # Add adducts
   ann = PafAnnotation.make_precursor(adducts=["+H", "+Na"], charge=2)
   print(ann.serialize())

.. testoutput::

   p[M+H+Na]^2

.. testcode::

   # Add charge state
   ann = PafAnnotation.make_peptide("y", 5, charge=2)
   print(ann.serialize())

.. testoutput::

   y5^2

.. testcode::

   # Add mass error (in ppm)
   ann = PafAnnotation.make_peptide("y", 5, mass_error=1.2, mass_error_unit="ppm")
   print(ann.serialize())

.. testoutput::

   y5/1.2ppm

.. testcode::

   # Add confidence score
   ann = PafAnnotation.make_peptide("y", 5, confidence=0.95)
   print(ann.serialize())

.. testoutput::

   y5*0.95

Complex Annotations
~~~~~~~~~~~~~~~~~~~

Combine multiple parameters to create complex annotations.

.. testcode::

   from paftacular import PafAnnotation

   # Everything at once
   ann = PafAnnotation.make_peptide(
       "y", 5,
       neutral_losses=["-H2O"],
       isotopes=["+i13C"],
       adducts=["+H", "+Na"],
       charge=2,
       mass_error=-0.55,
       mass_error_unit="ppm",
       confidence=0.85,
       is_auxiliary=True,
       analyte_reference=2
   )
   print(ann.serialize())

.. testoutput::

   &2@y5-H2O+i13C[M+H+Na]^2/-0.55ppm*0.85

With Sequences
~~~~~~~~~~~~~~

.. testcode::

   from paftacular import PafAnnotation

   # Peptide ion with sequence
   ann = PafAnnotation.make_peptide("y", 3, sequence="PEP")
   print(ann.serialize())

.. testoutput::

   y3{PEP}

.. testcode::

   # Internal fragment with sequence
   ann = PafAnnotation.make_internal(2, 5, sequence="EPTI")
   print(ann.serialize())

.. testoutput::

   m2:5{EPTI}

.. testcode::

   # Modified immonium ion
   ann = PafAnnotation.make_immonium("M", modification="Oxidation")
   print(ann.serialize())

.. testoutput::

   IM[Oxidation]

Other Ion Types
~~~~~~~~~~~~~~~

.. testcode::

   from paftacular import PafAnnotation

   # Chemical formula
   ann = PafAnnotation.make_formula("C6H12O6")
   print(ann.serialize())

.. testoutput::

   f{C6H12O6}

.. testcode::

   # SMILES notation
   ann = PafAnnotation.make_smiles("CN=C=O")
   print(ann.serialize())

.. testoutput::

   s{CN=C=O}

.. testcode::

   # Named compound
   ann = PafAnnotation.make_named_compound("Urocanic Acid")
   print(ann.serialize())

.. testoutput::

   _{Urocanic Acid}

.. testcode::

   # Unknown ion
   ann = PafAnnotation.make_unknown(label=42)
   print(ann.serialize())

.. testoutput::

   ?42

Serialization (Round-trip)
--------------------------

.. testcode::

   import paftacular as pft

   # Parse and serialize back
   original = "y5-H2O^2/1.2ppm*0.95"
   ann = pft.parse(original)
   serialized = ann.serialize()
   print(serialized)

.. testoutput::

   y5-H2O^2/1.2ppm*0.95

Export to Dictionary
--------------------

.. code-block:: python

   import paftacular as pft

   # Export annotation as a versioned dictionary (for JSON)
   ann = pft.parse("y5-H2O^2/1.2ppm*0.95")
   data = ann.to_dict()


Peptacular Integration
----------------------

`paftacular` integrates with `peptacular <https://github.com/tacular-omics/peptacular>`_
to convert its ``Fragment`` objects directly into ``PafAnnotation`` objects via ``to_mzpaf()``.

.. code-block:: python

   import peptacular as pt
   import paftacular as pft

   # Generate fragment ions from a peptide sequence
   fragments = pt.fragment("PEPTIDE", charges=[1, 2], ion_types=["b", "y"])

   # Convert each fragment to a PafAnnotation
   annotations = [pft.to_mzpaf(f) for f in fragments]
   for ann in annotations[:3]:
       print(ann.serialize())
   # b1{P}
   # b2{PE}
   # b3{PEP}

Optional parameters let you attach a mass error and confidence score to each annotation:

.. code-block:: python

   ann = pft.to_mzpaf(
       fragments[0],
       confidence=0.95,
       mass_error=1.2,
       mass_error_unit="ppm",
   )
   print(ann.serialize())   # b1{P}/1.2ppm*0.95

Pass ``include_sequence=False`` to omit the embedded sequence from the annotation
(useful when sequence context is not needed):

.. code-block:: python

   ann = pft.to_mzpaf(fragments[0], include_sequence=False)
   print(ann.serialize())   # b1

For more details, see the `PSI mzPAF specification <https://www.psidev.info/>`_.

Resolving Analyte Context
-------------------------

Without sequence context, ``get_mass()`` of a peptide, internal or precursor
annotation returns only the ion offset and modifiers, and ``mz()`` raises
``PaftacularError``. Use ``resolve()`` for a complete fragment calculation. It returns a
new annotation and leaves the original unchanged.

.. testcode::

   original = pft.parse("y2")
   resolved = original.resolve("PEPTIDE")
   print(resolved.sequence)
   print(round(resolved.mz(), 4))
   print(original.sequence)

.. testoutput::

   DE
   263.0874
   None

A mapping selects analytes by their mzPAF references. An omitted reference
selects key 1. A supplied string is the full analyte for the current annotation.
Resolution supports peptide fragments (including the d, v and w side-chain
series), internal fragments, and precursors. It rejects missing references,
out-of-range positions, internal fragments that include a terminus, and
conflicting embedded sequences.
The analyte's charge is excluded from the selected sequence. The annotation's
own charge and adducts determine the final charged species.

.. testcode::

   resolved = pft.parse("2@b2").resolve({1: "AAAA", 2: "PEPTIDE"})
   print(resolved.sequence)
   print(resolved.serialize())

.. testoutput::

   PE
   2@b2

Resolved context is deliberately separate from mzPAF text. Use ``to_dict()``
when persisting it. Parsing ``resolved.serialize()`` requires resolving the
analyte again before computing a complete mass.

Structured Errors and Batch Parsing
-----------------------------------

Every error caused by bad input is a ``PaftacularError``, which is a
``ValueError``. ``PafParseError`` is a ``PaftacularError`` with ``text``, ``reason``,
``position``, and ``annotation_index`` attributes. Positions and indices start
at zero. Syntax errors point to unexpected content or the end of incomplete
input. Semantic errors identify the start of the affected annotation.

Parsing does not look up reference molecule names. Calculating an ``r[...]``
ion or a ``-[...]`` loss whose name is in neither the mzPAF reference list nor
Unimod raises ``PafUnknownReferenceError``, with the name in ``name``. It is a
``PaftacularError`` and also a ``KeyError``, which 1.3.2 raised for ``r[...]``.

``iter_parse()`` processes an iterable lazily. Wrap it in ``list()`` to
collect the results. Each result preserves the original record text and index.
If any annotation in a record is malformed, that record has an error and no
partial annotations. Later records still run. Empty records succeed with an
empty annotation tuple. Programming errors such as passing non-string records
are not converted into parse errors.

.. testcode::

   results = list(pft.iter_parse(["y2,b3", "y2,b3!", "p"]))
   print([result.ok for result in results])
   print(results[1].index, results[1].error.annotation_index, results[1].error.position)

.. testoutput::

   [True, False, True]
   1 1 5

Versioned Interchange
---------------------

``to_dict()`` exports schema version 1, including an ion type discriminator,
component fields, optional mass error, and resolved sequence context.
``PafAnnotation.from_dict()`` reconstructs that representation, rejecting
unknown versions, missing or extra fields, invalid types, nonfinite numbers,
and malformed mzPAF fields. It does not require optional dependencies to
reconstruct data. Chemical and ProForma interpretation happens when properties
are calculated.

.. testcode::

   import json
   resolved = pft.parse("b2").resolve("PEPTIDE")
   data = resolved.to_dict()
   restored = pft.PafAnnotation.from_dict(json.loads(json.dumps(data)))
   print(data["schema_version"], data["ion"]["type"])
   print(restored == resolved)
   print(restored.sequence)

.. testoutput::

   1 PeptideIon
   True
   PE

Calculation and Serialization Conventions
-----------------------------------------

``get_mass()`` is the charged species mass. ``comp()`` counts nuclei, so summing
its elemental masses requires subtracting ``charge`` electron masses (a negative
charge adds them) to compare with ``get_mass()``. ``mz()`` divides by the absolute charge. Upstream tabulated ion offsets are rounded,
so numerical comparisons should allow approximately one microdalton.

A named modification (Unimod, PSI-MOD, RESID, XLMOD, GNO) adds its listed database mass,
the same rule as peptacular. Fragment and precursor ions, including neutral losses and isotope
peaks, agree with peptacular 5 to 1e-9 Da. A global isotope label (``<13C>``, ``<15N>``) uses composition in both packages.
Unimod reference names (``r[Hex]``) use the listed 6-decimal mass, while mzPAF
reference-list entries (``r[TMT6plex]``) keep their exact formula masses. Labile
modifications (``{Glycan:Hex}``) count only for precursor ions: fragments lose them. The listed
mass is rounded (Oxidation is 15.994915, its composition gives 15.99491462), so the mass of
``comp()`` can differ from ``get_mass()`` by up to about 1e-6 Da for a named modification.
Composition is used only for modifications with no listed mass, such as formulas.
Formula ions already describe the charged species' atoms. Their adducts label
charge carriers without adding atoms. Other supported ions add the specified
adduct atoms, or implicit hydrogens when no adduct is given.

A generic ``+i`` uses the carbon-13 minus carbon-12 mass shift. In compositions,
isotope substitutions consume ordinary or explicit monoisotopic atoms when
available. Genuine deficits remain visible when context is incomplete.
Average-isotopomer masses remain undefined and raise ``PaftacularError``.
SMILES isotope labels are retained during calculation.
Charged SMILES molecules are rejected. Charge carriers belong in mzPAF adducts.

Precursor conversion with ``to_mzpaf()`` retains sequence context when
``include_sequence=True``. This context is preserved by structured export,
while the mzPAF string remains ``p`` with its ordinary modifiers and charge.

There is a known convention difference between the internal-cleavage table in
mzPAF 1.0.1 section 4.4.4 and the physical ion definitions in tacular and
peptacular. ``make_internal(ion_type=...)`` and the public correction tables
retain the specification convention for compatibility. Explicit
``InternalFragment`` cleavage fields and ``to_mzpaf()`` preserve the source
ion's physical composition using explicit neutral gains and losses. For
example, a physical ``ax`` fragment is represented as ``m2:4-H2``. The mass
of an existing mzPAF string is never reinterpreted based on an inferred series.

Numeric output uses decimal notation compatible with the grammar and retains
float precision. Mass-based losses retain at least five decimal places for
compatibility, adding more when necessary. Round trips preserve meaning, not
necessarily the original spelling or redundant zero-isotope annotations.
