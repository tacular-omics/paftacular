.. image:: _static/paftacular_logo.png
   :alt: Paftacular Logo
   :align: center
   :width: 300px

|

.. raw:: html

   <div style="text-align: center; margin-bottom: 5px;">
      <a href="https://github.com/pgarrett-scripps/paftacular/actions/workflows/python-package.yml"><img src="https://github.com/pgarrett-scripps/paftacular/actions/workflows/python-package.yml/badge.svg" alt="Python package"></a>
      <a href="https://codecov.io/github/paftacular-omics/paftacular"><img src="https://codecov.io/github/paftacular-omics/paftacular/graph/badge.svg?token=1CTVZVFXF7" alt="codecov"></a>
      <a href="https://paftacular.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/paftacular/badge/?version=latest" alt="Documentation Status"></a>
      <a href="https://badge.fury.io/py/paftacular"><img src="https://badge.fury.io/py/paftacular.svg" alt="PyPI version"></a>
      <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.12+-blue.svg" alt="Python 3.12+"></a>
      <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License"></a>
   </div>

|

.. raw:: html

   <div style="text-align: center; font-size: 1.0em; margin-bottom: 20px;">
      Welcome to Paftacular's documentation! Paftacular is a Python library for parsing and serializing <strong>mzPAF</strong> (Peak Annotation Format)

   </div>

mzPAF is a specification from the <a href="https://www.psidev.info/">Proteomics Standards Initiative (PSI)</a>
that provides a compact, human-readable notation for describing fragment ion types, chemical modifications, 
charge states, mass errors, and confidence scores.

Features
--------

* **mzPAF parsing**: Handles parsing / serializing of mzPAF strings
* **Properties**: Supports calculating mass and composition of annotated ions
* **Type-Annotations**: typed.py file for static type checking
* **Caching**: serialization and parsing results are cached for performance (when applicable)
* **Integrated**: Integrated with peptacular, such that peptacular can output mzPAF annotations for fragment ions

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
