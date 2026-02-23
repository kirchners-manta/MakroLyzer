Surface Atoms [1]_
=================

This analyzer identifies surface atoms of a polymer particle and classifies them into:

- nonpolar surface atoms
- polar surface atoms
- ring surface atoms

Surface atoms are detected with the GITIM method [2]_ as implemented in **pytim** [3]_ on an **MDAnalysis** universe [4]_.
Ring surface atoms are carbon atoms that belong to cycles matching the selected ring-size criterion.

Command line
------------
.. line-block::
  ``-sA``
  ``--surfaceAtoms``
      Identify and count nonpolar, polar, and ring atoms on the particle surface.
      Optionally provide an output filename.
      *Default: surfaceAtoms.csv*

  ``--sA-ring-size``
      Ring-size filter used for ring atom assignment.
      Accepts a single ring size (e.g. ``6``) or an inclusive range (e.g. ``[4,7]``).
      Minimum allowed ring size is ``3``.
      *Default: 6*

Example
^^^^^^^
.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -sA surfaceAtoms.csv --sA-ring-size [4,7]

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
   * - Nonpolar Surface Atoms
     - Number of surface atoms classified as nonpolar
   * - Polar Surface Atoms
     - Number of surface atoms classified as polar
   * - Ring Surface Atoms
     - Number of surface carbon atoms that belong to rings matching ``--sA-ring-size``

.. [1]  Drysch, K.; Dawer, Y.; Zaby, P.; Buchmüller, K.; Dick, L.; Mutzel, P.; Hollóczki, O.; Kirchner, B. MakroLyzer: A Graph-Based Software to Comb through Molecular Hairballs Using the Example of Nanoplastics. J. Phys. Chem. B 2025
       DOI: 10.1021/acs.jpcb.5c06175
.. [2]  Sega, M.; Kantorovich, S.; Jedlovszky, P.; Jorge, M. The generalized identification of truly interfacial molecules (ITIM) algorithm for nonplanar interfaces. J. Chem. Phys. 2013, 138, 044110.
       DOI: 10.1063/1.4774361
.. [3]  Sega, M.; Hantal, G.; Fábián, B.; Jedlovszky, P. Pytim: A Python Package for the Interfacial Analysis of Molecular Simulations. J. Comput. Chem. 2018, 39, 2118-2125.
       DOI: 10.1002/jcc.25384
.. [4]  Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. J. Comput. Chem. 2011, 32, 2319-2327.
       DOI: 10.1002/jcc.21787
