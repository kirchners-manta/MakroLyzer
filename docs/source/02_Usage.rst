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

Parameters that need to be provided:
==================================================

XYZ or LMP file (required)
------------------
.. line-block::
  ``-xyz``
  ``-lmp``
      Path to the XYZ or LAMMPS trajectory file containing the macromolecule structure.
      *Required. Either -xyz or -lmp must be provided.*

Box Size 
------------------
.. line-block::
  ``-bs``
      Box size of the macromolecule in Angstroms. This is required if periodic boundary conditions are used.
      *Optional. Default: infinity*

Timesteps to calculate
------------------
.. line-block::
  ``-nth``
      Specifies the for which how many timesteps from the trajectory the given parameters should be calculated.
      *Optional. Default: 1*
      
Pattern file 
------------------
.. line-block::
  ``-p``
  ``--pattern-file``
      Path to the pattern file containing the patterns to search for in the macromolecule structure.
      *Optional. Default: None*


Parameters that can be calculated:
=================================

Radius of gyration
------------------
.. line-block::
  ``-r``
  ``--Rg-file``
      Calculate the radius of gyration of the provided structure. Provide a file name for the output if desired. 
      *Optional. Default: radiusOfGyration.csv*

The radius of gyration is defined as:

.. math::

   R_g^2 = \frac{1}{N} \sum_{j=1}^{N} (\vec{r_j} - \vec{r_{cm}})^2

where :math:`\vec{r_j}` is the position vector of atom :math:`j`, :math:`\vec{r_{cm}}` is the center of mass position vector, and :math:`N` is the number of atoms.

Asphericity parameter
------------------
.. line-block::
  ``-as``
  ``--asphericity-file``
      Calculate the asphericity parameter of the macromolecule structure. Provide a file name for the output if desired. 
      *Optional. Default: asphericity.csv*

Anisotropy parameter
------------------
.. line-block::
  ``-an``
  ``--anisotropy-file``
      Calculate the anisotropy parameter of the macromolecule structure. Provide a file name for the output if desired. 
      *Optional. Default: anisotropyFactor.csv*

Order parameter
------------------
.. line-block::
  ``-op``
  ``--order-file``
      Calculate the order parameter of the macromolecule structure. Provide a file name for the output if desired. 
      *Optional. Default: orderParameter.csv*

End to End Distance
------------------
.. line-block::
  ``-e2e``
  ``e2e-file``
      Calculate the end to end distance of the macromolecule structure. Provide a file name for the output if desired. 
      *Optional. Default: endToEndDistances.csv*


