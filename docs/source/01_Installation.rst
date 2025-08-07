Installation
=================
To install Makrolyzer, clone the repository to your local machine

.. code-block:: bash

   git clone git@github.com:kirchners-manta/MakroLyzer.git

The python versions supported by MakroLyzer are 3.10 to 3.12.
For using the code, the dependencies listed in the ``pyproject.toml`` file need to be installed.

To build the MakroLyzer, using a virtual environment is recommended.
One can either use ``conda`` or ``venv`` to create a virtual environment.

Setting up a virtual environment with ``conda``:
------------------
.. code-block:: bash

   conda create -n makrolyzer python=3.10
   conda activate makrolyzer
   
Setting up a virtual environment with ``venv``:
------------------
.. code-block:: bash

   python3 -m venv makrolyzer-env
   source makrolyzer-env/bin/activate

After activating the virtual environment, install the dependencies using pip:

.. code-block:: bash

   pip install .


Call help page
=================
To call the MakroLyzer help, run the following command in the terminal:

.. code-block:: bash

   MakroLyzer -h
