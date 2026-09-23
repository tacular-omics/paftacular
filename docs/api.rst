API Reference
=============

Annotation
----------

.. automodule:: paftacular.annotation
   :members:
   :undoc-members:
   :show-inheritance:

Parser
------

.. automodule:: paftacular.parser
   :members:
   :undoc-members:
   :show-inheritance:

Resolution and Errors
---------------------

.. automodule:: paftacular.resolution
   :members:

.. automodule:: paftacular.errors
   :members:

Structured Interchange
----------------------

.. automodule:: paftacular.serialization
   :members: to_dict, from_dict

Annotation Types
----------------

.. automodule:: paftacular.comps
   :members:
   :undoc-members:
   :inherited-members:

.. py:data:: paftacular.IonType

   Type alias for any ion-type component of an annotation: the union of
   :class:`~paftacular.comps.PeptideIon`, :class:`~paftacular.comps.InternalFragment`,
   :class:`~paftacular.comps.ImmoniumIon`, :class:`~paftacular.comps.ReferenceIon`,
   :class:`~paftacular.comps.NamedCompound`, :class:`~paftacular.comps.ChemicalFormula`,
   :class:`~paftacular.comps.SMILESCompound`, :class:`~paftacular.comps.UnknownIon`
   and :class:`~paftacular.comps.PrecursorIon`. Also importable from
   ``paftacular.comps``.
   
Constants
---------

.. automodule:: paftacular.constants
   :members:
   :undoc-members:
   :show-inheritance:

.. py:data:: paftacular.INTERNAL_MASS_DIFFS
   :type: dict[tuple[str, str], str | None]

   Formula offset of an internal fragment for each pair of N-terminal and
   C-terminal cleavage types, for example ``("b", "y")``. ``None`` means no
   offset. Defined in ``paftacular.constants``.

Peptacular Conversion
---------------------

Requires the ``peptacular`` extra.

.. automodule:: paftacular.conversion
   :members: to_mzpaf

Utilities
---------

.. automodule:: paftacular.util
   :members:
   :undoc-members:
   :show-inheritance:
