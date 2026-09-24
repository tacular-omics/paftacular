.. image:: _static/paftacular_logo.png
   :alt: Paftacular Logo
   :align: center
   :width: 300px

|

.. raw:: html

   <div style="text-align: center; margin-bottom: 5px;">
      <a href="https://github.com/tacular-omics/paftacular/actions/workflows/ci.yml"><img src="https://github.com/tacular-omics/paftacular/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
      <a href="https://codecov.io/github/tacular-omics/paftacular" > 
         <img src="https://codecov.io/github/tacular-omics/paftacular/graph/badge.svg?token=lZDTvRrnuq"/> 
      </a>
      <a href="https://paftacular.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/paftacular/badge/?version=latest" alt="Documentation Status"></a>
      <a href="https://badge.fury.io/py/paftacular"><img src="https://badge.fury.io/py/paftacular.svg" alt="PyPI version"></a>
      <a href="https://doi.org/10.5281/zenodo.19076277"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.19076277.svg" alt="Zenodo DOI"></a>
      <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.12+-blue.svg" alt="Python 3.12+"></a>
      <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License"></a>
   </div>

|

.. raw:: html

   <div style="text-align: center; font-size: 1.0em; margin-bottom: 20px;">
      Welcome to Paftacular's documentation! Paftacular is a Python library for parsing and serializing <strong>mzPAF</strong> (Peak Annotation Format). mzPAF is a specification from the <a href="https://www.psidev.info/">Proteomics Standards Initiative (PSI)</a>
      that provides a compact, human-readable notation for describing fragment ion types, chemical modifications, 
      charge states, mass errors, and confidence scores.

   </div>


Features
--------

* **mzPAF parsing**: Handles parsing / serializing of mzPAF strings
* **Properties**: Supports calculating mass and composition of annotated ions
* **Type-Annotations**: typed.py file for static type checking
* **Caching**: serialization and parsing results are cached for performance (when applicable)
* **Integrated**: Integrated with peptacular, such that peptacular can output mzPAF annotations for fragment ions

Related packages
----------------

* `tacular <https://tacular.readthedocs.io/>`_ provides the element, amino acid and
  reference-molecule lookups paftacular uses for masses.
* `peptacular <https://peptacular.readthedocs.io/>`_ parses **ProForma** peptide
  sequences. Install ``paftacular[peptacular]`` to compute masses of annotations
  that carry a sequence, and to have peptacular emit mzPAF for its fragment ions.

Quick Example
-------------

.. testcode::

   import paftacular as pft

   ann = pft.parse("y5")
   print(ann.ion_type.series)
   print(ann.get_mass())

.. testoutput::

   y
   19.017841466621

New to paftacular? Start with :doc:`installation` and :doc:`quickstart`.
See :doc:`usage` for the full guide, including creating annotations programmatically,
computing masses/compositions, and round-tripping to mzPAF strings.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   usage
   migration
   mcp
   api
   changelog
   citation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
