Saturation
==========

.. figure:: ../figures/structureModification/03sat3.png
   :width: 800
   :align: center

This modifier saturates polymer chain ends by adding terminal atoms. 
Terminal carboxylic carbon atoms receive OH groups and terminal amide nitrogens receive H.

Command line
------------
.. line-block::
  ``-s``
  ``--saturation``
      Saturate chain ends.
      *Default: saturatedPolymers.xyz*

Example
-------
.. code-block:: bash

   MakroLyzer -xyz polymer.xyz -s -saturatedPolymer.xyz


Output
------
An XYZ file containing the saturated structure. 
For trajectories, multiple output files are created with suffixes if files already exist.

