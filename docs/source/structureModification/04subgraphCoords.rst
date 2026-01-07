Subgraph Coordinates [1]_
=========================

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

.. [1]  Drysch, K.; Dawer, Y.; Zaby, P.; Buchmüller, K.; Dick, L.; Mutzel, P.; Hollóczki, O.; Kirchner, B. MakroLyzer: A Graph-Based Software to Comb through Molecular Hairballs Using the Example of Nanoplastics. J. Phys. Chem. B 2025
       DOI: 10.1021/acs.jpcb.5c06175