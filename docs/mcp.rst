AI clients through MCP
======================

The optional MCP server exposes paftacular to AI applications through local
stdio. It provides nine tools, four reference resources, and two reusable
prompts. Calculations reuse the public Python API and run locally.

Installation and connection
---------------------------

MCP support is included starting with paftacular 1.3.0. Install the optional extra:

.. code-block:: bash

   pip install 'paftacular[mcp]'
   paftacular-mcp --version

The ``mcp`` extra includes the official Python MCP SDK v2 and peptacular.
Use ``paftacular[mcp,smiles]`` to add SMILES support. The ``all`` extra includes MCP.
Ordinary imports and the base install do not load the SDK or Pydantic.
For a development checkout, use ``pip install '.[mcp]'`` instead.

Configure the host to launch ``paftacular-mcp`` from that environment. Using
an absolute Python path avoids PATH differences in desktop applications.
For hosts that use an ``mcpServers`` JSON configuration:

.. code-block:: json

   {
     "mcpServers": {
       "paftacular": {
         "command": "/absolute/path/to/venv/bin/python",
         "args": ["-m", "paftacular.mcp"]
       }
     }
   }

Replace the path with the environment where you installed the extra.
On Windows, use the environment's ``Scripts/python.exe`` path. Hosts with a
different configuration format need the same executable and argument array.
The host starts the process and connects automatically. Running the command
in a terminal alone waits for MCP messages and does not open an interactive
chat. ``--help`` and ``--version`` work without the SDK installed.

Alternatively, hosts can launch an isolated, version-pinned installation with uv:

.. code-block:: json

   {
     "mcpServers": {
       "paftacular": {
         "command": "uvx",
         "args": ["--from", "paftacular[mcp]==1.3.0", "paftacular-mcp"]
       }
     }
   }

The ``uvx`` executable must be available to the host. Its first launch installs
the package and dependencies, which requires access to the package index.

The server itself requires no provider API key. It opens no HTTP listener.
Remote-only hosts need a separate future HTTP deployment option.

Tools
-----

All tools advertise input and output schemas and read-only, idempotent,
closed-world behavior. Every tool except ``get_capabilities`` takes one
``request`` object.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Tool
     - Purpose
   * - ``get_capabilities``
     - Discover dependency versions, installed integrations, conventions, and limits.
   * - ``parse_annotations``
     - Inspect and normalize a list of text records. Each record can contain multiple annotations.
   * - ``build_annotation``
     - Add typed modifiers, charge, analyte reference, confidence, and mass error to a bare ion.
   * - ``resolve_annotation``
     - Select a peptide, internal, or precursor sequence from full ProForma analyte context.
   * - ``serialize_annotation``
     - Validate a complete core interchange dictionary and return canonical mzPAF text.
   * - ``calculate_ion``
     - Calculate requested mass, m/z, formula, or composition for one annotation.
   * - ``calculate_ions``
     - Calculate a list of requests with independent outcomes and stable ordering.
   * - ``generate_fragments``
     - Generate proper a/b/c/x/y/z terminal fragments for selected positions and charges.
   * - ``match_mz``
     - Compare supplied candidates with an observed m/z using ppm or Th tolerance.

Calculate a fragment
--------------------

Call ``calculate_ion`` with:

.. code-block:: json

   {
     "request": {
       "annotation": "y3^2",
       "analyte": "PEPTIDE",
       "properties": ["mass", "mz", "composition"]
     }
   }

The selected sequence is ``IDE``. The default calculation mode is ``complete``.
Peptide, internal, and precursor annotations need an embedded sequence or
supplied analyte. Missing context produces an actionable error. For mixtures,
use ``analytes`` with integer references instead of ``analyte``:

.. code-block:: json

   {
     "request": {
       "annotation": "2@y3^2",
       "analytes": [{"reference": 2, "sequence": "PEPTIDE"}]
     }
   }

Supply either ``analyte`` or ``analytes``. Duplicate references are invalid.
An omitted annotation reference selects analyte 1. The annotation controls
charge, and an embedded sequence must agree with the selected analyte fragment.

Use ``mode: "offsets"`` explicitly to calculate only ion offsets and modifiers.
This mode applies to peptide, internal, and precursor ions and rejects supplied
analytes. Any embedded sequence is ignored and labeled accordingly in the result.

Scientific outputs
------------------

Results label monoisotopic charged-species mass as ``mass_da`` and m/z as
``mz_th``. Charge and context source are explicit. ``mass_basis`` distinguishes
a complete charged species from offsets and modifiers. Average-mass
calculations are not exposed in this MCP version.

``properties`` defaults to ``["mass", "mz"]``. Formula and composition are
optional. Element and isotope names are composition keys, with signed counts
preserved. A numeric mass loss can allow mass and m/z while preventing
composition. Such a result has ``status: "partial"`` and per-property errors,
while retaining the valid calculated values. Inspect these fields before
presenting an answer.

The response includes the core annotation dictionary. Its ``schema_version``
is independent from the outer ``response_schema_version``. Resolved sequence
context is preserved in that dictionary but omitted from canonical mzPAF
text. Calculation, resolution, and matching requests accept either mzPAF text
or a complete dictionary in ``annotation``. Pass the returned dictionary to
calculate again with its resolved context. A text-only result needs analyte
context supplied again.

Parsing validates notation and structural constraints. It does not establish
that every named reference is known or that every composition is physically
possible. Formula ions already describe the charged species atoms. Internal
cleavage corrections follow the conventions explained in :doc:`usage`.

Generation and matching
-----------------------

``generate_fragments`` defaults to b/y series with charge 1 and all positions
from 1 through analyte length minus 1. It excludes the full-length peptide.
Select ``series``, ``charges``, and ``positions`` to reduce the results.
Duplicate selections are removed while preserving order. Output is ordered
by series, position, then charge. Modified analytes use the same context
resolution as individual calculations. Internal fragments and custom losses
can be supplied through ``calculate_ions``.

``match_mz`` takes ``observed_mz``, a list of candidate annotation/context
objects, ``tolerance``, and ``tolerance_unit`` (``ppm`` or ``Th``). The default
tolerance is 10 ppm. It returns every candidate in input order, including
errors, plus matching indices sorted by absolute ppm error. Signed errors are
observed minus theoretical, with theoretical m/z as the ppm denominator.
Matching is inclusive at the tolerance boundary. It does not identify a
peptide, search a database, or assign a confidence probability.

Resources and prompts
----------------------

The server lists ``paftacular://guide``, ``paftacular://conventions``,
``paftacular://examples``, and ``paftacular://capabilities``. These provide
bundled workflow guidance and current capabilities without remote lookups.

The ``analyze_fragment`` prompt accepts ``annotation`` and ``analyte`` strings
and guides tool-based analysis. The ``review_annotations`` prompt accepts
``records_json``, a JSON array encoded as a string, and guides validation.
Prompts provide instructions to the host and do not themselves call a model.

Errors and limits
-----------------

Domain failures return an MCP error result with a structured error code and
message. Parse errors include zero-based character and annotation positions.
Malformed tool argument schemas are rejected by the SDK. Batches retain
per-record failures without failing the whole call. Partial single-ion
results remain successful MCP responses with explicit property errors.

Calls accept at most 100 records or generated fragments, 16 KiB per text,
256 KiB of JSON-encoded application inputs, and 512 KiB of structured JSON
output. Unicode escapes count toward aggregate JSON limits. These are
application limits, not a transport-level request firewall. Oversized work
returns an error with guidance to split the batch. No results are silently
truncated. The MCP result includes a text JSON copy for compatible clients.
Parsing is additionally limited to 1000 annotations across all records in a call.
Empty text records parse successfully to an empty annotation list.

Diagnostics go to stderr. Protocol messages alone use stdout. The server
accepts data values and does not expose file access or arbitrary execution.
Its chemistry handlers run sequentially on one event loop with bounded
batches. Large spectral-library processing remains better suited to the
Python API. No hard wall-clock deadline is imposed on an individual chemistry
operation.

Development checks
------------------

.. code-block:: bash

   uv sync --all-extras
   uv run pytest tests/test_mcp.py tests/test_installation.py
   just pre-release

Tests use the official SDK Client in memory and launch the actual command and
module as stdio subprocesses from outside the checkout. Installation CI covers
the base package, chemistry extras, MCP, MCP with SMILES, and all extras using
minimum and newer compatible dependencies. A Windows job checks the connection
and entry points as well.
