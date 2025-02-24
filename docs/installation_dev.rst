.. highlight:: shell

Installation (Development Mode)
===============================

.. _Use a Virtual Environment:

Use a Virtual Environment
-------------------------

Use the Standard Library `venv` package to create a Python virtual environment.

.. code-block:: console

    $ python -m venv env
    $ .\env\Scripts\activate  # (Windows)
    $ source env/bin/activate # (Linux)

The Python version of the virtual environment will match the Python version of your base
interpreter, see `Installing packages using pip and virtual environments`_.

.. important::
   The command to activate the environment depends on the virtual environment manager
   of our choice, your operating system, and your shell. On Linux, you must `source`.

Upgrade to the latest version of pip in the virtual environment (Recommended).

.. code-block:: console

    $ python -m pip install --upgrade pip
    $ pip -V # >= 21.3?

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.
As an alternative to using the `venv` package, you may create and activate a Conda-environment.
This step requires a Python-distribution with Conda package management https://www.anaconda.com/.

.. code-block:: console

    $ conda env create --name py310hl_gaps_pub python=3.10
    $ conda activate py310hl_gaps_pub

Or use the `environment.yml` file:

.. code-block:: console

    $ mamba env create -f environment.yml


.. _pip: https://pip.pypa.io
.. _Installing packages using pip and virtual environments: https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


.. _Development Installation Instructions:

Installation from Sources
-------------------------

The sources for `HL-gaps-pub` can be downloaded from the `Gitlab repo`_.
You can clone the public repository as follows:

.. code-block:: console

    $ git clone https://HL_gaps_pub/ifam418/HL-gaps-pub.git

Having a SSL problem?

.. code-block:: console

    $ env GIT_SSL_NO_VERIFY=true git clone https://HL_gaps_pub/ifam418/HL-gaps-pub.git


.. _Reproduce the Environment:

Reproduce the Environment and GitHooks
......................................

.. code-block:: console

    $ make install-dev
    $ make install-githooks


Install the Package in Editable Mode
....................................

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ pip install -e .


Optional Installation Steps
---------------------------

Build Documentation
...................

This step requires a `sphinx` installation. If not available on the system, install the development environment
which contains the necessary packages (`pip install -r requirements/dev.txt`). make sure to activate the environment.

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
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    |
    ├── hl_gaps_pub
    │   └── __init__.py    <- The python package
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


.. _Gitlab repo: https://gitlab.cc-asp.fraunhofer.de/ifam418/Cic
.. _Docker with Proxy: https://docs.docker.com/network/proxy
.. _Docker behind a Firewall: https://stackoverflow.com/questions/29630480/apt-get-in-docker-behind-corporate-proxy
.. _Docker Tricks with Firewalls: https://mandie.net/2017/12/10/docker-for-windows-behind-a-corporate-web-proxy-tips-and-tricks

