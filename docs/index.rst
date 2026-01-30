.. paftacular documentation master file

paftacular
==========

A Python library for parsing and serializing **mzPAF** (Peak Annotation Format), 
a standardized format for annotating mass spectrometry fragment ions in peptide/proteomics analysis.

mzPAF is a specification from the `Proteomics Standards Initiative (PSI) <https://www.psidev.info/>`_ 
that provides a compact, human-readable notation for describing fragment ion types, chemical modifications, 
charge states, mass errors, and confidence scores.

Features
--------

* **mzPAF parsing**: Handles all major annotation types from the specification
* **Mass calculations**: Supports both monoisotopic and average mass calculations
* **Composition analysis**: Full elemental composition tracking
* **Type-safe**: Comprehensive type hints
* **Caching**: Object caching for improved performance

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
