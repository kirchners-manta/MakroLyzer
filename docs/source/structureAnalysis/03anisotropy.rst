Anisotropy Factor [1]_
======================

The anisotropy factor :math:`\kappa^2` is an invariant of the radius of gyration tensor G.
It is a measure of symmetry and dimensionality of a polymer.

Definition
----------
It is defined as

.. math::

   \kappa^2 = 1 - 3 \frac{\lambda_1\lambda_2 + \lambda_2\lambda_3 + \lambda_3\lambda_1}
                      {(\lambda_1 + \lambda_2 + \lambda_3)^2}

where :math:`\lambda_1 \ge \lambda_2 \ge \lambda_3` are the eigenvalues of the
radius of gyration tensor.

Command Line Input
------------
.. line-block::
  ``-af``
  ``--anisotropyFactor``
      Calculate :math:`\kappa^2`. Optionally provide an output filename.
      *Default: anisotropyFactor.csv*

Example
^^^^^^^
.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -af myAnisotropyOutput.csv

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
   * - Anisotropy Factor (:math:`\kappa^2`)
     - Anisotropy factor for the full graph

.. [1]  Drysch, K.; Dawer, Y.; Zaby, P.; Buchmüller, K.; Dick, L.; Mutzel, P.; Hollóczki, O.; Kirchner, B. MakroLyzer: A Graph-Based Software to Comb through Molecular Hairballs Using the Example of Nanoplastics. J. Phys. Chem. B 2025
       DOI: 10.1021/acs.jpcb.5c06175