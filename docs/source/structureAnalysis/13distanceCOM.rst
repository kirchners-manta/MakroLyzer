Distance from Center of Mass
============================

This analysis measures the distance from the mass-weighted center of a
selected particle to the center of mass of each selected solvent molecule.
If the particle selection matches several molecules, their atoms contribute
to a single reference center of mass.

The two selections use MakroLyzer's subgraph-selection syntax. For example,
``CHON`` selects molecules containing exactly those four elements in any
proportions, and ``H2O`` selects molecules with two hydrogen atoms and one
oxygen atom. Supply the cubic box length in Å with ``-bs`` to calculate
minimum-image distances under periodic boundary conditions. Without a box
length, the analyzer uses ordinary Euclidean distances.

Command line
------------
.. line-block::
  ``-DistCOM PARTICLE_SEL SOLVENT_SEL``
      Calculate the distance from the particle center of mass to each selected solvent center of mass.
      *Default output: DistanceFromCOM.csv*

  ``-DistCOM-file FILE``
      Set a custom output filename.

Example
^^^^^^^
.. code-block:: bash

   MakroLyzer -xyz solvated_polymer.xyz -DistCOM CHON H2O -bs 50
   MakroLyzer -xyz solvated_polymer.xyz -DistCOM CHON H2O -DistCOM-file water_distances.csv -o results

Output
------
The first CSV line records the reference selection. The second line lists
the solvent molecules selected in the first analyzed frame. MakroLyzer
writes the zero-based frame index followed by the distances in Å as soon
as it analyzes each frame.
Distances follow the order of the header. Molecules are tracked by their
atom IDs. A blank field indicates that the molecule did not match the
solvent selection in that frame. The columns are fixed for the run, so
only molecules selected in the first analyzed frame have columns.

.. code-block:: text

   Reference,CHON
   H2O #1 (Å),H2O #2 (Å)
   0,3.250,5.100
   1,,4.875
