Asphericity Parameter [1]_
==========================

The asphericity parameter :math:`b` measures how far the structure deviates from a sphere. 
Values close to zero indicate a near-spherical shape.

Definition
----------
It is defined as

.. math::

   b = \lambda_1 - \frac{1}{2}(\lambda_2 + \lambda_3)

where :math:`\lambda_1 \ge \lambda_2 \ge \lambda_3` are the eigenvalues of the
radius of gyration tensor.

Command line
------------
.. line-block::
  ``-as``
  ``--asphericityParameter``
      Calculate :math:`b`. Optionally provide an output filename.
      *Default: asphericityParameter.csv*

Example
^^^^^^^
.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -as myAsphericityOutput.csv

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
   * - Asphericity Parameter (b)
     - Asphericity value for the full graph


.. [1]  Drysch, K.; Dawer, Y.; Zaby, P.; Buchmüller, K.; Dick, L.; Mutzel, P.; Hollóczki, O.; Kirchner, B. MakroLyzer: A Graph-Based Software to Comb through Molecular Hairballs Using the Example of Nanoplastics. J. Phys. Chem. B 2025
       DOI: 10.1021/acs.jpcb.5c06175