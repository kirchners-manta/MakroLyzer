What the user needs to provide
=================================================
The user needs to provide a structure or a trajecory file of a macromolecule in the
``.lmp`` or ``.xyz`` format.
In case periodic boundary conditions are used, the user also needs to provide the box size.

Additionally, the flags fot the parameters that should be calculated need to be set.
One example input line is the following:

.. code-block:: bash

   MakroLyzer -xyz PolyEthylene.xyz -bs 80.0 -r -as -af

This command will calculate the radius of gyration and the asphericity parameter of 
the structure in the ``PolyEthylene.xyz`` file with a box size of 80.0 Angstroms.

.. note::
    So far only cubic box sizes are supported.

Parameters that can be calculated
=================================================
The following parameters need to be provided or can be calculated with MakroLyzer:

XYZ or LMP file
------------------
- ``-xyz``: Path to the XYZ file containing the macromolecule structure.
- ``-lmp``: Path to the LAMMPS trajectoryfile containing the macromolecule structure.
One of these two parameters is **required**.

Box Size
------------------
- ``-bs``: Box size of the macromolecule in Angstroms. This is required if periodic boundary conditions are used.

Timesteps to calculate
------------------
- ``-nth``: In case the user provides a trajectory file, this parameter can be used to specify the timesteps to calculate
            the parameters for.
            The default value is 1, which means that every timestep will be calculated.

Pattern file 
------------------
- ``-p``: Path to the pattern file containing the patterns to search for in the macromolecule structure.
- ``--repeatingUnits-file--``: Output file name for the repeating units found in the macromolecule structure. Default is ``repeatingUnits.csv``.
