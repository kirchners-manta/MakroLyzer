Mean Squared Displacement (MSD)
===============================

The mean squared displacement (MSD) measures the average squared distance that
atoms move over time. MakroLyzer calculates atomic MSD for the selected atoms,
not the MSD of molecular centers of mass. This also applies when selecting
molecules with ``-sel H2O``.

Definition
----------
The MSD is defined as

.. math::

   \mathrm{MSD}(\tau) = \left\langle \frac{1}{N} \sum_{i=1}^N
   \left\|\mathbf{r}_i(t + \tau) - \mathbf{r}_i(t)\right\|^2 \right\rangle_t

Here, :math:`\mathbf{r}_i(t)` is the position of atom :math:`i` at time
:math:`t`, :math:`\tau` is the lag time, and :math:`N` is the number of selected
atoms. The sum divided by :math:`N` averages over atoms; the angle brackets
average over all available time origins for each lag. The selected atom
identities and their ordering must remain unchanged throughout the trajectory.

Command line
------------
.. line-block::
  ``-MSD TIME`` or ``--MSD TIME``
      Calculate atomic MSD up to the specified maximum correlation time.
      Provide a positive value in the same units as ``--timestep``.
  ``--MSD-file FILE``
      Specify an output filename for the MSD results.
      *Default: MSD.csv*
  ``--timestep TIME``
      Required when using MSD. Provide the positive time interval between
      saved trajectory frames, before applying ``-nth``.

With ``-nth n``, the interval between analyzed frames is
:math:`n \times \mathrm{timestep}`. The maximum correlation time is unchanged.

Coordinates are unwrapped before applying the analysis stride. XYZ and wrapped
LAMMPS coordinates require a fixed box size supplied with ``-bs``. Native
LAMMPS ``xu/yu/zu`` coordinates are already unwrapped; GROMACS unwrapping uses
the box information in the trajectory.

Example
^^^^^^^
For an XYZ trajectory in a cubic box with side length 50 Å, saved every 1 ps,
calculate atomic MSD for the selected water molecules up to a lag of 100 ps:

.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -bs 50 -MSD 100 --MSD-file myMSDOutput.csv --timestep 1 -sel H2O

The time unit is supplied by the user: the program does not infer it from the
numerical values in this command.

Output
------
The output file contains one row per positive sampled lag within the requested
maximum correlation time. Lags without available frame pairs are omitted,
including when the trajectory is too short. Lag zero is not written.

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Column
     - Description
   * - Correlation Time
     - Sampled lag time, in the units supplied by the user.
   * - MSD / Å²
     - Squared displacement averaged over selected atoms and available time origins.
