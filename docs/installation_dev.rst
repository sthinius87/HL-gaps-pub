.. _installation:

Development Installation Instructions
=====================================

The sources for `HL-gaps-pub` can be downloaded from the `Gitlab repo`_.
You can clone the public repository as follows:

.. code-block:: console

    $ git clone https://github.com/sthinius87/HL-gaps-pub.git

.. _Reproduce the Environment:

Use the `environment.yml` file to create the environment:

.. code-block:: console

    $ mamba env create -f environment.yml

Install the package in development mode:

.. code-block:: console

    $ cd HL-gaps-pub
    $ mamba activate py310hl_gaps_pub
    $ pip install -e .


.. _mamba: https://github.com/conda-forge/miniforge
.. _pip: https://pip.pypa.io
.. _Installing packages using pip and virtual environments: https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


Optional Installation Steps
---------------------------

Build Documentation
...................

This step requires a `sphinx` installation, shipped with the installation described above. Make sure to activate the environment.
Since we have a Linux package, build it under Linux. Building the documentation under Windos will fail due to the usage of links to tho notebooks.

.. code-block:: console

    $ make docs


Install Jupyter Kernel for you Notebooks
........................................

This step requires an `ipython` installation. If not available on the system, install the development environment
which contains the necessary packages (`pip install -r requirements/dev.txt`).

.. code-block:: console

    $ ipython kernel install --user --name=py310hl_gaps_pub


Project Organization
--------------------

::

    ├── LICENSE
    ├── Makefile           
    ├── README.rst         <- The top-level README for developers using this project.
    ├── data
    │   ├── coconut        <- Data from the COCONUT database.
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   ├── raw            <- The original, immutable data dump.
    │   └── test           <- Data for the code testing.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── notebooks          <- Jupyter notebooks.
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    ├── scripts            <- scripts that can be run on compute cluster.    
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    |
    ├── hl_gaps_pub
    │   └── __init__.py    <- The python package
    ├── tests              <- Testing of the package.
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


.. _Gitlab repo: https://gitlab.cc-asp.fraunhofer.de/ifam418/Cic
.. _Docker with Proxy: https://docs.docker.com/network/proxy
.. _Docker behind a Firewall: https://stackoverflow.com/questions/29630480/apt-get-in-docker-behind-corporate-proxy
.. _Docker Tricks with Firewalls: https://mandie.net/2017/12/10/docker-for-windows-behind-a-corporate-web-proxy-tips-and-tricks

