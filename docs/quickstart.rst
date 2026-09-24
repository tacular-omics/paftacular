Quickstart
==========

.. testsetup:: *

   import paftacular as pft

Install paftacular with ``pip install paftacular`` (see :doc:`installation` for
extras).

There are two parsing functions:

* ``parse``: Parses exactly one mzPAF annotation and returns a ``PafAnnotation``. Raises ``PafParseError`` for zero or several annotations.
* ``parse_multi``: Parses comma-separated mzPAF annotations. Always returns a list of ``PafAnnotation``.

Every error caused by bad input is a ``PaftacularError`` (a ``ValueError``).

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
   ann = pft.parse("y5-H2O^2/1.2ppm*0.95")
   print(ann.charge)
   print(ann.mass_error.value)
   print(ann.confidence)

.. testoutput::

   2
   1.2
   0.95

Compute the mass of an annotated ion. Without a sequence, ``get_mass()`` is the
ion offset only, and ``mz()`` needs a sequence (embedded or from ``resolve()``):

.. testcode::

   ann = pft.parse("y5")
   print(ann.get_mass())

.. testoutput::

   19.017841150651

Serialize back to mzPAF:

.. testcode::

   print(pft.parse("y5-H2O^2").serialize())

.. testoutput::

   y5-H2O^2

Next, read :doc:`usage` for every ion type, modifications, isotopes, adducts and
the peptacular integration.
