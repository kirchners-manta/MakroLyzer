Chemical Formulas
=================

This analyzer reports the chemical formulas of all connected subgraphs and their counts per frame. 
Element symbols are sorted alphabetically in each formula.

Command line
------------
.. line-block::
  ``-f``
  ``--formula``
      Calculate chemical formulas. Optionally provide an output filename.
      *Default: chemicalFormulas.csv*

Output
------
The output file contains one row per formula and frame:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Column
     - Description
   * - Frame
     - Frame index in the trajectory
   * - Chemical Formula
     - Formula string (e.g. C6H12O6)
   * - Count
     - Number of subgraphs with this formula
