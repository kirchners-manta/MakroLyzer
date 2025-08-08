The user needs to provide a structure or a trajectory file of a macromolecule in the
``.lmp`` or ``.xyz`` format.
In case periodic boundary conditions are used, the user also needs to provide the box size.

Additionally, the flags for the parameters that should be calculated need to be set.
One example input line is the following:

.. code-block:: bash

   MakroLyzer -xyz PolyEthylene.xyz -bs 80.0 -r -as -af

This command will calculate the radius of gyration and the asphericity parameter of 
the structure in the ``PolyEthylene.xyz`` file with a box size of 80.0 Angstroms.

.. note::
    So far only cubic box sizes are supported.

Parameters that need to be provided
=======================================

XYZ or LMP file
^^^^^^^^^^^^^^^
.. line-block::
  ``-xyz``
  ``-lmp``
      Path to the XYZ or LAMMPS trajectory file containing the macromolecular structure.
      *Required. Either -xyz or -lmp must be provided.*


Box Size
^^^^^^^^^^^^
.. line-block::
  ``-bs``
      Box size of the macromolecule in Angstroms. This is required if periodic boundary conditions are used.
      *Optional. Default: infinity*


Timesteps to calculate
^^^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-nth``
      Specifies the for which how many timesteps from the trajectory the given parameters should be calculated.
      *Optional. Default: 1*
      

Pattern file
^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-p``
  ``--pattern-file``
      Path to the pattern file containing the patterns to search for in the macromolecular structure.
      *Optional. Default: None*


Parameters that can be calculated
=====================================

Radius of gyration
^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-r``
  ``--Rg-file``
      Calculate the radius of gyration of the provided structure. Provide a file name for the output if desired. 
      *Optional. Default: radiusOfGyration.csv*

The radius of gyration is defined as:

.. math::

   R_\mathrm{g}^2 = \frac{1}{N} \sum_{j=1}^{N} \left(\vec{r}_j - \vec{r}_{\mathrm{cm}}\right)^2

where :math:`\vec{r_j}` is the position vector of atom :math:`j`, :math:`\vec{r_{cm}}` is the center of mass position vector, and :math:`N` is the number of atoms.
The radius of gyration is calculated for each connected substructure in the macromolecule and for the whole macromolecule.

The gyration tensor is defined as:

.. math::

   S_{m,n} = \frac{1}{N} \sum_{i=j}^N
             \left(r_{j}^{(m)} - r_{\mathrm{cm}}^{(m)}\right)
             \left(r_{j}^{(n)} - r_{\mathrm{cm}}^{(n)}\right)

where :math:`r_{j}^{(m)}` is the :math:`m`-th component of the position vector of atom :math:`j`, and :math:`r_{\mathrm{cm}}^{(m)}` is the :math:`m`-th component of the center of mass position vector.
:math:`\lambda_1 \ge \lambda_2 \ge \lambda_3` are the eigenvalues of the tensor.
They are used to calculate the asphericity parameter and the anisotropy factor.

Asphericity parameter
^^^^^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-as``
  ``--asphericity-file``
      Calculate the asphericity parameter of the macromolecular structure. Provide a file name for the output if desired. 
      *Optional. Default: asphericity.csv*

The asphericity parameter is defined as:

.. math::

   b = \lambda_1 - \frac{1}{2}\left(\lambda_2 + \lambda_3\right)


Anisotropy factor
^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-af``
  ``--anisotropy-file``
      Calculate the anisotropy parameter of the macromolecular structure. Provide a file name for the output if desired. 
      *Optional. Default: anisotropyFactor.csv*

The anisotropy factor is defined as:

.. math::

   \kappa^2 = 1 - 3 \frac{\lambda_1 \lambda_2 + \lambda_2 \lambda_3 + \lambda_3 \lambda_1}
                   {(\lambda_1 + \lambda_2 + \lambda_3)^2}


Order parameter
^^^^^^^^^^^^^^^^^
.. line-block::
  ``-op b:n:v``
  ``--order-file``
      Calculate the order parameter of the macromolecular structure. This is calculated for each cube defined by **b** and **n**.
      **b** is the box size in Angstroms,
      **n** is the number of cubes into which the box is divided in x, y, and z direction (thus overall n^3 boxes),
      **v** is the number of atoms in a structure backbone, that is used to form one molecular vector.
      Provide a file name for the output if desired. 
      *Optional. Default: orderParameter.csv*

The order parameter is defined as:

.. math::

   S = \frac{1}{2} \langle 3 \cos^2{\theta_m} - 1 \rangle

where :math:`\theta_m` is the angle between the molecular axis and the reference axis (director), 
and :math:`\langle \cdots \rangle` denotes the average over all molecules.

Dihedral angles 
^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-d``
  ``dihedral-range {abs,nonabs}``
  ``--dihedral-file``
      Calculate the dihedral angles of the macromolecular structure along the backbones of the substructures.
      **dihedral-range** specifies whether to calculate the absolute dihedral angles or the non-absolute dihedral angles.
      Provide a file name for the output if desired. 
      *Optional. Default: dihedrals.csv*


Cis trans count
^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-ct``
  ``--CisTrans-file``
      Calculate the cis and trans counts of the macromolecular structure along the backbones of the substructures.
      Provide a file name for the output if desired. 
      *Optional. Default: cisTransCounts.csv*


End-to-End Distance
^^^^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-e2e``
  ``--e2e-file``
      Calculate the end-to-end distance of the macromolecular structure. Provide a file name for the output if desired. 
      *Optional. Default: endToEndDistances.csv*

The end to end distance is defined as:

.. math::

   \vec{R}_n = \sum_{i=1}^n \vec{r}_i

where :math:`\vec{r}_i` is the bond vector of two atoms in the backbone of the macromolecule, and :math:`n` is the number of bond vectors in the backbone.


Number of hydrogen bonds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-hb A:AH:AD:B``
  ``--hbonds-file``
      Calculate the number of hydrogen bonds in the macromolecular structure. 
      **A** is the atom type of the acceptor, **AH** is the maximum hydrogen atom acceptor atom distance, **AD** is the maximum hydrogen atom donor atom distance, and **B** is the maximum acceptor-hydrogen-donor angle in degrees.
      Provide a file name for the output if desired. 
      *Optional. Default: hydrogenBonds.csv*


Number of connected substructures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-noSub``
  ``--noSub-file``
      Calculates the number of connected substructures in the macromolecular structure. Provide a file name for the output if desired.
      *Optional. Default: noSubGraphs.csv*


Formulas of connected substructures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-f``
  ``--formula-file``
      Calculates the formulas of the connected substructures in the macromolecular structure. Provide a file name for the output if desired.
      *Optional. Default: chemicalFormulas.csv*


Ring and strand count
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. line-block::
  ``-RScount``
  ``--RingStrandCount-file``
      Calculates the number of rings and strands in the macromolecular structure. Provide a file name for the output if desired.
      *Optional. Default: ringAndStandCount.csv*