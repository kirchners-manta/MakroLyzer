Installation
=================

The only requirement to install MakroLyzer is `Git <https://git-scm.com/docs/user-manual>`_.
To install MakroLyzer, first clone the repository to your local machine:

.. code-block:: bash

   # SSH 
   git clone git@github.com:kirchners-manta/MakroLyzer.git

   # HTTPS 
   git clone https://github.com/kirchners-manta/MakroLyzer.git


Then change into the directory into the cloned repository.
For using the code, the dependencies listed in the ``pyproject.toml`` file need to be installed, which is done automatically when
following the steps below.
`Python <https://www.python.org/>`_ 3.10 or higher are supported by MakroLyzer.

To set up the MakroLyzer, using a virtual environment is recommended.
One can either use ``conda`` or ``venv`` to create a virtual environment.

Why use a virtual environment?
------------------------------
- Prevents dependency conflicts between projects
- Lets you manage project-specific dependencies
- Keeps your global Python environment clean


Set up with conda
---------------------
Conda is a package and environment manager.

.. code-block:: bash

   conda create -n makrolyzer python=3.12
   conda activate makrolyzer

   # If you prefer conda-forge:
   conda create -n makrolyzer -c conda-forge python=3.12

*Windows:* activation is the same: ``conda activate makrolyzer``.

A user guide for conda can be found `here <https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html>`_.

   
Set up with venv
--------------------
Venv is built into Python and creates lightweight virtual environments.

.. code-block:: bash

   python3 -m venv makrolyzer-env

   # macOS/Linux
   source makrolyzer-env/bin/activate
   
   # Windows (PowerShell)
   .\makrolyzer-env\Scripts\Activate.ps1

A user guide for venv can be found `here <https://docs.python.org/3/library/venv.html>`_.


Install MakroLyzer
=================
After activating the virtual environment, install the dependencies using pip.
Make sure you are in the ``polymer_analysis`` directory of the cloned repository.

.. code-block:: bash

   pip install .


Call help page
=================
To call the MakroLyzer help, run the following command in the terminal:

.. code-block:: bash

   MakroLyzer -h
