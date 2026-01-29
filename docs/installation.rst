Installation
============

Requirements
------------

* Python ≥3.12
* `tacular <https://github.com/tacular-omics/tacular>`_ - Element and reference molecule lookups


Installing from PyPI
--------------------

.. code-block:: bash

   pip install paftacular
   pip install paftacular[sequence] # for peptide sequence support
   pip install paftacular[smiles] # with smiles support
   pip install paftacular[all] # with all optional dependencies

Optional Dependencies
---------------------

* ``pysmiles`` - For SMILES notation support
* ``peptacular`` - For Mass/Composition calculations via included sequence

Installing from Source
----------------------

.. code-block:: bash

   git clone https://github.com/tacular-omics/paftacular.git
   cd paftacular
   pip install .

Development Installation
------------------------

For development with all dev dependencies:

.. code-block:: bash

   git clone https://github.com/tacular-omics/paftacular.git
   cd paftacular
   just install
