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

* **Comprehensive mzPAF support**: Handles all major annotation types from the specification
* **Parse complex annotations**: Fragment ions, neutral losses, isotopes, adducts, charge states, mass errors, and confidence scores
* **Mass calculations**: Supports both monoisotopic and average mass calculations
* **Composition analysis**: Full elemental composition tracking with ProForma-style formulas
* **Type-safe**: Comprehensive type hints

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
