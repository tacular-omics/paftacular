Usage
=====

.. testsetup:: *

   import paftacular as pft

Quick Start
-----------

There are 3 parsing methods available:
* ``parse``: Parses a single or multiple comma-separated mzPAF annotations. Returns a single ``PafAnnotation`` or a list of them.
* ``parse_multi``: Parses multiple comma-separated mzPAF annotations. Always returns a list of ``PafAnnotation``.
* ``parse_single``: Parses a single mzPAF annotation. Returns a single ``PafAnnotation``. Raises ValueError if multiple annotations are provided.

.. testcode::

   import paftacular as pft

   # Parse a simple peptide ion
   ann = pft.parse("y5")
   print(ann.ion_type.series)
   print(ann.ion_type.position)

.. testoutput::

   y
   5

.. testcode::

   # Parse with modifications
   ann = pft.parse_single("y5-H2O^2/1.2ppm*0.95")
   print(ann.charge)
   print(ann.mass_error.value)
   print(ann.confidence)

.. testoutput::

   2
   1.2
   0.95

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
   ann = pft.parse("y5[M+2H+Na]")
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

Mass Errors and Confidence
---------------------------

.. code-block:: python

   import paftacular as pft

   # Mass error in ppm
   ann = pft.parse("y5/1.2ppm")
   print(ann.mass_error.value)  # 1.2
   print(ann.mass_error.unit)   # "ppm"

   # Mass error in daltons
   ann = pft.parse("y5/-0.003da")

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
   print(ann.ion_type.series)     # IonSeries.Y
   print(ann.ion_type.position)   # 5
   print(ann.neutral_mods)        # Tuple with H2O loss
   print(ann.isotopes)            # Tuple with 13C isotope
   print(ann.adducts)             # Tuple with H and Na
   print(ann.charge)              # 2
   print(ann.mass_error.value)    # -0.55
   print(ann.confidence)          # 0.85

Mass Calculations
-----------------

.. code-block:: python

   import paftacular as pft

   # Get calculated masses
   ann = pft.parse("y5")
   print(ann.monoisotopic_mass)   # Monoisotopic mass
   print(ann.average_mass)        # Average mass
   print(ann.mass)                # Defaults to monoisotopic

   # With modifications
   ann = pft.parse("y5-H2O")
   print(ann.monoisotopic_mass)   # Mass adjusted for water loss

   # Get composition
   ann = pft.parse("f{C6H12O6}")
   comp = ann.composition  # Counter with ElementInfo keys
   print(ann.formula)              # "C6H12O6"
   print(ann.proforma_formula)     # ProForma-style formula

Parsing Multiple Annotations
-----------------------------

.. testcode::

   import paftacular as pft

   # Parse comma-separated list (can also use parse_multi)
   annotations = pft.parse("b2, y3-H2O, p^2") 
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
   ann = PafAnnotation.make_precursor(adducts=["+H", "+Na"])
   print(ann.serialize())

.. testoutput::

   p[M+H+Na]

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
   ann = PafAnnotation.make_internal(2, 5, sequence="PEPTIDE")
   print(ann.serialize())

.. testoutput::

   m2:5{PEPTIDE}

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

   # Export annotation as dictionary (e.g., for JSON)
   ann = pft.parse("y5-H2O^2/1.2ppm*0.95")
   data = ann.as_dict()


For more details, see the `PSI mzPAF specification <https://www.psidev.info/>`_.
