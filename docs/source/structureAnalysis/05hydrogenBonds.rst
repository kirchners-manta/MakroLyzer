Hydrogen Bonds [1]_
===================

.. figure:: ../figures/structureAnalysis/03hb.png
   :width: 300
   :align: center
   :alt: Dihedral plot

This analyzer counts hydrogen bonds using geometric cutoffs.

Command line
------------
.. line-block::
  ``-hb A:b:c:alpha``
  ``--hydrogenBonds``
      Count hydrogen bonds for one or more acceptor types.
      ``A`` is the acceptor element symbol
      ``b`` is the maximum hydrogen atom (H) - acceptor atom (A) distance (float)
      ``c`` maximum donor atom (D) - acceptor atom distance (float)
      ``alpha`` is the maximum angle in degrees between the D-H and H-A vectors (float)

  ``-hb-file``
  ``--hbonds-file``
      Output filename.
      *Default: hydrogenBonds.csv*

Example
^^^^^^^

.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -hb O:2.4:3.4:30 -hb-file myHbondOutput.csv

Output
------
The output file contains one row per frame and cutoff:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Column
     - Description
   * - Frame
     - Frame index in the trajectory
   * - Element Type
     - Acceptor element symbol
   * - H-Acceptor dist
     - Cutoff distance ``b``
   * - Donor-Acceptor dist
     - Cutoff distance ``c``
   * - Angle cutoff
     - Cutoff angle ``alpha`` (degrees)
   * - Number of Hydrogen Bonds
     - Count for this cutoff

.. [1]  Drysch, K.; Dawer, Y.; Zaby, P.; Buchmüller, K.; Dick, L.; Mutzel, P.; Hollóczki, O.; Kirchner, B. MakroLyzer: A Graph-Based Software to Comb through Molecular Hairballs Using the Example of Nanoplastics. J. Phys. Chem. B 2025
       DOI: 10.1021/acs.jpcb.5c06175
