Polyethylene Surface Functionalization 
===========================================

.. figure:: ../figures/structureModification/02bfuncPE3.png
   :width: 800
   :align: center

This modifier functionalizes polyethylene (PE) surfaces by replacing selected CH2 units at
the particle surface with CO, ester, or amide groups. Surface atoms are identified using
the GITIM algorithm implemented in **pytim** [1]_ [2]_ on an **MDAnalysis** universe [3]_. A random subset of
eligible backbone carbons is chosen. End-group carbons are excluded. 

Functionalization types
-----------------------
- ``CO`` replaces CH2 with a carbonyl group. 
- ``COH``replaces CH2 with a hydroxy group.
- ``ester`` replaces CH2-CH2 with COO.
- ``amide`` replace CH2-CH2 with CONH.

Modes
-----
- ``Random`` mode selects a target number of CH2 sites based on the requested percentage and then chooses positions at random while respecting ``neighbor_exclusion`` to keep functional groups separated. A neighbor_exclusion value of three means, that between two functional groups are at least three untouched CH2 units.

Command line
------------
.. line-block::
  ``-funcPEsurf random:percentage:func_type:neighbor_exclusion``
  ``--functionalizePEsurface``
      Random functionalization of the PE surface by replacing CH2 groups with ``func_type``.
      ``func_type`` is ``CO``, ``COH``, ``ester``, or ``amide``.
      ``percentage`` is a float from 0 to 100.
      ``neighbor_exclusion`` is the minimum C-C hop distance between functional groups.

  ``-funcPE-file``
  ``--functionalizePE-file``
      Output filename.
      *Default: functionalizedPE.xyz*

Example
^^^^^^^
.. code-block:: bash

   MakroLyzer -xyz PE.xyz -funcPEsurf random:15:CO:2 -funcPE-file PE_CO.xyz

.. note::

    - Chain-end carbons (with only one carbon neighbor) are excluded from functionalization.

Output
------
An XYZ file containing the modified structure.


.. [1]  Sega, M.; Kantorovich, S.; Jedlovszky, P.; Jorge, M. The generalized identification of truly interfacial molecules (ITIM) algorithm for nonplanar interfaces. J. Chem. Phys. 2013, 138, 044110.
       DOI: 10.1063/1.4774361
.. [2] Sega, M.; Hantal, G.; Fábián, B.; Jedlovszky, P. Pytim: A Python Package for the Interfacial Analysis of Molecular Simulations. J. Comput. Chem. 2018, 39, 2118-2125.
       DOI:  10.1002/jcc.25384
.. [3]  Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. J. Comput. Chem. 2011, 32, 2319-2327.
       DOI: 10.1002/jcc.21787
