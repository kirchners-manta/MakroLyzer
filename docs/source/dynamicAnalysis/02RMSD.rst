Root Mean Square Deviation (RMSD)
================================

The root mean square deviation (RMSD) measures how far the selected atoms
move from their positions in the first analyzed frame. MakroLyzer reports
the raw RMSD and RMSD center of mass (COM) corrected one, where the COM
is calculated for all selected atoms.

Definition
----------
The raw RMSD is defined as

.. math::

   \mathrm{RMSD}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^N
   \left\|\mathbf{r}_i(t) - \mathbf{r}_i(t_0)\right\|^2}

Here, :math:`\mathbf{r}_i(t)` is the position of selected atom :math:`i` at time :math:`t`,
:math:`N` is the number of selected atoms, and :math:`t_0` denotes the first
analyzed frame. For the COM-corrected result, the mass-weighted center of
mass :math:`\mathbf{R}_{\mathrm{COM}}` is subtracted from each frame:

.. math::

   \mathrm{RMSD}_{\mathrm{COM}}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^N
   \left\|[\mathbf{r}_i(t) - \mathbf{R}_{\mathrm{COM}}(t)]
   - [\mathbf{r}_i(t_0) - \mathbf{R}_{\mathrm{COM}}(t_0)]\right\|^2}

Both RMSDs weight atoms equally. The COM correction removes collective translation
of the entire selection, it does not center each selected molecule separately.

Command line
------------
.. line-block::
  ``-RMSD [FILE]`` or ``--RMSD [FILE]``
      Calculate raw and COM-corrected RMSD for each analyzed frame.
      Optionally provide an output filename. *Default: RMSD.csv*


Example
^^^^^^^
For an XYZ trajectory in a cubic box with side length 50 Å, calculate RMSD
for the atoms in the selected water molecules:

.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -bs 50 -RMSD -sel H2O

Output
------
The CSV file contains one row per analyzed frame, including the reference
frame, where both RMSDs are zero. RMSD values are written to three decimal
places.

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Column
     - Description
   * - Frame
     - Zero-based trajectory frame index, retaining the original indices when using ``-nth``.
   * - RMSD / Å
     - Raw atom-averaged RMSD relative to the reference frame.
   * - COM-corrected RMSD / Å
     - Atom-averaged RMSD after subtracting each frame's mass-weighted center of mass.
