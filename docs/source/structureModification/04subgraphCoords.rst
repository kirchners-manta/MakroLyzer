Subgraph Coordinates
====================

This modifier writes each connected subgraph (molecule) as a separate XYZ file.

Command line
------------
.. line-block::
  ``-sub``
  ``--subgraph-coords``
      Output filename base.
      *Default: subgraphCoordinates.xyz*

Example
^^^^^^^
.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -sub mySubgraphCoordinates.xyz

Output
------
Each subgraph is written to ``<base>_frag<ID>.xyz`` with an incremented fragment ID
starting at 1. Each file contains element symbols and coordinates.
