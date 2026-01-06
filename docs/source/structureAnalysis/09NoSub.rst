Connected Substructures
=======================

This analyzer counts the number of connected subgraphs (molecules) and classifies each as a chain or ring. 
A subgraph is considered a chain if its backbone has free ends after removing 1-order nodes, otherwise it is classified as a ring.

Command line
------------
.. line-block::
  ``-NoSub``
  ``--NoSubgraphs``
      Calculate molecule, chain, and ring counts. Optionally provide an output filename.
      *Default: noSubGraphs.csv*

Example
^^^^^^^
.. code-block:: bash

    MakroLyzer -xyz polymer.xyz -NoSub myNoSubOutput.csv

Output
------
The output file contains one row per frame:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Column
     - Description
   * - Frame
     - Frame index in the trajectory 
   * - Molecule count
     - Number of connected subgraphs
   * - Chain count
     - Number of chain subgraphs
   * - Ring count
     - Number of ring subgraphs
