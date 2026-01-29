Usage
=====

.. testsetup:: *

   import paftacular as pft

Quick Start
-----------

.. testcode::

   import paftacular as pft

   # Parse a simple peptide ion
   ann = pft.parse_single("y5")
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
   ann = pft.parse_single("b2")
   print(ann.ion_type.series)
   
   ann = pft.parse_single("y3")
   print(ann.ion_type.series)

.. testoutput::

   b
   y

.. testcode::

   # With peptide sequence
   ann = pft.parse_single("y3{PEP}")
   print(ann.ion_type.sequence)

.. testoutput::

   PEP

.. testcode::

   # Internal fragments
   ann = pft.parse_single("m2:5{PEPTIDE}")
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
   ann = pft.parse_single("p^2")  # Doubly-charged precursor

   # Immonium ions
   ann = pft.parse_single("IK")          # Lysine immonium
   ann = pft.parse_single("IK[Acetyl]")  # Modified

   # Reference compounds
   ann = pft.parse_single("r[TMT126]")

   # Chemical formula
   ann = pft.parse_single("f{C6H12O6}")

   # SMILES notation
   ann = pft.parse_single("s{CN=C=O}")

   # Unannotated/unknown peaks
   ann = pft.parse_single("?")
   ann = pft.parse_single("?42")  # With reporter number

Modifications
-------------

.. code-block:: python

   import paftacular as pft

   # Neutral losses (formula-based)
   ann = pft.parse_single("y5-H2O")
   ann = pft.parse_single("b3-NH3")
   ann = pft.parse_single("b2-H2O-NH3")  # Multiple losses

   # Neutral gains
   ann = pft.parse_single("y5+NH3")

   # Mass-based losses/gains
   ann = pft.parse_single("y5-17.03")
   ann = pft.parse_single("b2+22.98")

   # Reference-based
   ann = pft.parse_single("y5-[Adenine]")

Isotopes
--------

.. code-block:: python

   import paftacular as pft

   # Generic isotope peak
   ann = pft.parse_single("y5+i")

   # Specific element
   ann = pft.parse_single("y5+i13C")
   ann = pft.parse_single("b2+2i13C")  # Multiple isotopes

   # Average isotopomer
   ann = pft.parse_single("y5+iA")

Adducts
-------

.. code-block:: python

   import paftacular as pft

   # Single adduct
   ann = pft.parse_single("y5[M+Na]")

   # Multiple adducts
   ann = pft.parse_single("y5[M+2H+Na]")
   ann = pft.parse_single("p[M+NH4]")

Charge States
-------------

.. code-block:: python

   import paftacular as pft

   # Singly charged (default)
   ann = pft.parse_single("y5")
   print(ann.charge)  # 1

   # Doubly charged
   ann = pft.parse_single("y5^2")
   print(ann.charge)  # 2

   # Triply charged
   ann = pft.parse_single("b3^3")
   print(ann.charge)  # 3

Mass Errors and Confidence
---------------------------

.. code-block:: python

   import paftacular as pft

   # Mass error in ppm
   ann = pft.parse_single("y5/1.2ppm")
   print(ann.mass_error.value)  # 1.2
   print(ann.mass_error.unit)   # "ppm"

   # Mass error in daltons
   ann = pft.parse_single("y5/-0.003da")

   # Confidence score (0.0 to 1.0)
   ann = pft.parse_single("y5*0.95")
   print(ann.confidence)  # 0.95

   # Combined
   ann = pft.parse_single("y5/1.2ppm*0.95")

Complex Annotations
-------------------

.. code-block:: python

   import paftacular as pft

   # Everything at once
   ann = pft.parse_single("&2@y5-H2O+i13C[M+H+Na]^2/-0.55ppm*0.85")
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
   ann = pft.parse_single("y5")
   print(ann.monoisotopic_mass)   # Monoisotopic mass
   print(ann.average_mass)        # Average mass
   print(ann.mass)                # Defaults to monoisotopic

   # With modifications
   ann = pft.parse_single("y5-H2O")
   print(ann.monoisotopic_mass)   # Mass adjusted for water loss

   # Get composition
   ann = pft.parse_single("f{C6H12O6}")
   comp = ann.composition  # Counter with ElementInfo keys
   print(ann.formula)              # "C6H12O6"
   print(ann.proforma_formula)     # ProForma-style formula

Parsing Multiple Annotations
-----------------------------

.. testcode::

   import paftacular as pft

   # Parse comma-separated list
   annotations = pft.parse("b2, y3-H2O, p^2")
   for ann in annotations:
       print(ann.serialize())

.. testoutput::

   b2
   y3-H2O
   p^2

Serialization (Round-trip)
--------------------------

.. testcode::

   import paftacular as pft

   # Parse and serialize back
   original = "y5-H2O^2/1.2ppm*0.95"
   ann = pft.parse_single(original)
   serialized = ann.serialize()
   print(serialized)

.. testoutput::

   y5-H2O^2/1.2ppm*0.95

Export to Dictionary
--------------------

.. code-block:: python

   import paftacular as pft

   # Export annotation as dictionary (e.g., for JSON)
   ann = pft.parse_single("y5-H2O^2/1.2ppm*0.95")
   data = ann.as_dict()

mzPAF Format Reference
----------------------

The mzPAF format uses a compact notation for annotations::

   [&][analyte@]ion_type[modifications][^charge][/mass_error][*confidence]

**Examples:**

* ``y5`` - Y-series ion at position 5
* ``b2{PEP}`` - B-series ion with sequence PEP
* ``y5-H2O`` - Y5 with water loss
* ``y5^2`` - Doubly-charged Y5
* ``y5/1.2ppm`` - Y5 with 1.2 ppm mass error
* ``y5*0.95`` - Y5 with 95% confidence
* ``&2@y5-H2O^2/1.2ppm*0.95`` - Full annotation

For more details, see the `PSI mzPAF specification <https://www.psidev.info/>`_.
