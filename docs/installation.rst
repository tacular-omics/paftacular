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
   pip install paftacular[peptacular] # with peptacular integration
   pip install paftacular[smiles]     # with SMILES support
   pip install 'paftacular[mcp]'      # with the local MCP server and peptide support
   pip install paftacular[all]        # with all optional dependencies

MCP support requires paftacular 1.3.0 or newer. See :doc:`mcp` for client
configuration. The ``all`` extra includes MCP and its dependencies.

Optional Dependencies
---------------------

* ``pysmiles`` - For SMILES notation support
* ``peptacular`` - For Mass/Composition calculations via included sequence
* ``mcp`` extra - Official MCP SDK and peptacular for local AI client connections

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
