What the user needs to provide
=================================================
The user needs to provide a structure or a trajecory file of a macromolecule in the
``.lmp`` or ``.xyz`` format.
In case periodic boundary conditions are used, the user also needs to provide the box size.

Additionally, the flags fot the parameters that should be calculated need to be set.
One exaple input line is the following:

.. code-block:: bash

   MakroLyzer -xyz PolyEthylene.xyz -bs 80.0 -r -as -af

This command will calculate the radius of gyration and the asphericity parameter of 
the structure in the ``PolyEthylene.xyz`` file with a box size of 80.0 Angstroms.

.. note::
    So far only cubic box sizes are supported.

Parameters that can be calculated
=================================================
The following parameters can be calculated with MakroLyzer:

.. list-table:: Command Line Parameters
   :widths: 35 50 35
   :header-rows: 1

   * - Flag
     - Description
     - Default / Output File
   * - ``-xyz``, ``--xyzFile``
     - Path to the XYZ trajectory file
     - None
   * - ``-lmp``, ``--lmpFile``
     - Path to the LAMMPS trajectory file
     - None
   * - ``-nth``, ``--nthStep``
     - Read every nth step from the trajectory
     - 1
   * - ``-bs``, ``--BoxSize``
     - Box size for periodic boundary conditions
     - required if PBC used
   * - ``-p``, ``--patternFile``
     - Path to TXT file for finding repeating units
     - false
   * - ``--repeatingUnits-file``
     - Output file for repeating units
     - ``repeatingUnits.csv``
   * - ``-s``, ``--saturation``
     - Saturate ends of polymers
     - false
   * - ``--saturation-file``
     - Output file for saturated polymers
     - ``saturatedPolymers.xyz``
   * - ``-f``, ``--formula``
     - Get chemical formulas of the polymer
     - false
   * - ``--formula-file``
     - Output file for chemical formulas
     - ``chemicalFormulas.csv``
   * - ``-noSub``, ``--noSubgraphs``
     - Calculate number of subgraphs in the polymer
     - false
   * - ``--noSub-file``
     - Output file for number of subgraphs
     - ``noSubGraphs.csv``
   * - ``-e2e``, ``--endToEndDistance``
     - Calculate end-to-end distance
     - false
   * - ``--e2e-file``
     - Output file for end-to-end distances
     - ``endToEndDistances.csv``
   * - ``-d``, ``--dihedral``
     - Calculate dihedral angles
     - false
   * - ``--dihedral-range``
     - Dihedral angle range (``abs`` = 0-180, ``nonabs`` = -180-180)
     - ``abs``
   * - ``--dihedral-file``
     - Output file for dihedrals
     - ``dihedrals.csv``
   * - ``-ct``, ``--cisTrans``
     - Calculate cis and trans counts
     - false
   * - ``--CisTrans-file``
     - Output file for cis/trans counts
     - ``CisTrans.csv``
   * - ``-r``, ``--radiusOfGyration``
     - Calculate radius of gyration
     - false
   * - ``--Rg-file``
     - Output file for radius of gyration
     - ``radiusOfGyration.csv``
   * - ``-hb``, ``--hydrogenBonds``
     - List of tuples for H-bond detection (e.g., ``O:2.4:3.4:30 N:2.8:3.9:25``)
     - false
   * - ``--hbonds-file``
     - Output file for hydrogen bonds
     - ``hydrogenBonds.csv``
   * - ``-sub``, ``--subgraph-coords``
     - Get subgraph coordinates
     - false
   * - ``--subgraph-coord-file``
     - Output file for subgraph coordinates
     - ``subgraphCoordinates.csv``
   * - ``-af``, ``--anisotropyFactor``
     - Calculate anisotropy factor
     - false
   * - ``--anisotropy-file``
     - Output file for anisotropy factor
     - ``anisotropyFactor.csv``
   * - ``-as``, ``--asphericityParameter``
     - Calculate asphericity parameter
     - false
   * - ``--asphericity-file``
     - Output file for asphericity parameter
     - ``asphericityParameter.csv``
   * - ``-op``, ``--orderParameter``
     - Calculate order parameter S. Format: ``BoxSize:n:unitSize``
     - false
   * - ``--order-file``
     - Output file for order parameter
     - ``orderParameter.csv``
