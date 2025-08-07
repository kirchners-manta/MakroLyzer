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

The following parameters need to be provided:
==================================================

XYZ or LMP file (required)
------------------
- ``-xyz``: Path to the XYZ file containing the macromolecule structure.
- ``-lmp``: Path to the LAMMPS trajectoryfile containing the macromolecule structure.

Box Size 
------------------
- ``-bs``: Box size of the macromolecule in Angstroms. This is required if periodic boundary conditions are used.

Timesteps to calculate
------------------
- ``-nth``: In case the user provides a trajectory file, this parameter can be used to specify the timesteps to calculate
            the parameters for. The default value is 1, which means that every timestep will be calculated.

Pattern file 
------------------
- ``-p``: Path to the pattern file containing the patterns to search for in the macromolecule structure.
- ``--repeatingUnits-file--``: Output file name for the repeating units found in the macromolecule structure. Default is *repeatingUnits.csv*.


Parameters that can be calculated
=================================

Radius of gyration
------------------
- ``-r``: Calculate the radius of gyration of the macromolecule structure.
- ``--Rg-file``: Output file name for the radius of gyration. Default is *radiusOfGyration.csv*.

The radius of gyration is defined as:

.. math::

   R_g^2 = \sqrt{\frac{1}{N} \sum_{j=1}^{N} (\vec{r_j} - \vec{r_{cm}})^2}

where :math:`\vec{r_j}` is the position vector of atom :math:`j`, :math:`\vec{r_{cm}}` is the center of mass position vector, and :math:`N` is the number of atoms.

Asphericity parameter
------------------
- ``-as``: Calculate the asphericity parameter of the macromolecule structure.
- ``--asphericity-file``: Output file name for the asphericity parameter. Default is *asphericity.csv*.

Anisotropy parameter
------------------
- ``-an``: Calculate the anisotropy parameter of the macromolecule structure.
- ``--anisotropy-file``: Output file name for the anisotropy parameter. Default is *anisotropyFactor.csv*.

Order parameter
------------------
- ``-o``: Calculate the order parameter of the macromolecule structure.
- ``--order-file``: Output file name for the order parameter. Default is *orderParameter.csv*.