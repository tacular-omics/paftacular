Quickstart
==========

.. testsetup:: *

   import paftacular as pft

Install paftacular with ``pip install paftacular`` (see :doc:`installation` for
extras).

There are three parsing functions:

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

Compute the mass of an annotated ion:

.. testcode::

   ann = pft.parse("y5")
   print(ann.mass())

.. testoutput::

   19.017841466812

Serialize back to mzPAF:

.. testcode::

   print(pft.parse("y5-H2O^2").serialize())

.. testoutput::

   y5-H2O^2

Next, read :doc:`usage` for every ion type, modifications, isotopes, adducts and
the peptacular integration.
