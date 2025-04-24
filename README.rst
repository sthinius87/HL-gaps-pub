===========
HL-gaps-pub
===========

High throughput tight binding calculation of electronic HOMO-LUMO gaps and its prediction for natural compounds

.. image:: https://github.com/sthinius87/HL-gaps-pub/actions/workflows/main.yml/badge.svg
    :target: https://github.com/sthinius87/HL-gaps-pub/actions
    :alt: Build Status

.. image:: https://codecov.io/gh/sthinius87/HL-gaps-pub/graph/badge.svg?token=WFJUQSK6B9
    :target: https://codecov.io/gh/sthinius87/HL-gaps-pub
    :alt: Build Status

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: LICENSE
    :alt: License

.. image:: https://img.shields.io/badge/tox-py3.10 | py3.xx-blue.svg
    :target: https://github.com/sthinius87/HL-gaps-pub/blob/main/tox.ini
    :alt: tox environments

.. image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
    :target: https://pycqa.github.io/isort/
    :alt: Import Style: isort

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: Code Style: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit: enabled

.. image:: https://zenodo.org/badge/938057305.svg
    :target: https://doi.org/10.5281/zenodo.15113789
    :alt: Zenodo

.. image:: docs/figures/graphical_abstract_edit.png
    :alt: Graphical Abstract
    :width: 600px

#.. image:: figures/graphical_abstract_edit.png
#    :alt: Graphical Abstract
#    :width: 600px

**Predicting HOMO-LUMO Gaps of Natural Compounds with Machine Learning**

This package provides tools and workflows for predicting the HOMO-LUMO (HL) gap of natural compounds using machine learning models. It leverages the COCONUT database, RDKit for descriptor calculation, xTB for electronic structure calculations, and Common Workflow Language (CWL) for reproducible computational workflows. The models used are a Gradient Boosting Regressor (GBR) and a Multi-layer Perceptron Regressor (MLPR).

**Abstract**

This research develops a high-throughput, machine-learning approach to predict HOMO-LUMO (HL) gaps of natural compounds, a key property in cheminformatics and materials science. Using the COCONUT database (407,000 molecules) and RDKit descriptors, we compare Gradient Boosting Regression (GBR) and Multi-layer Perceptron Regression (MLPR) models. The computational workflow, managed by Toil and CWL on a Slurm cluster, uses xTB for electronic structure calculations with Boltzmann weighting. Molecular polarizability (especially SMR_VSA descriptors), aromatic rings, and functional groups like ketones were found to be crucial for HL-gap prediction. While the MLPR model showed good overall performance, accuracy varied across molecular subsets, particularly for molecules with aliphatic carboxylic acids, alcohols, and amines in complex electronic systems. This work demonstrates the potential of machine learning for HL-gap prediction, highlighting the importance of specific structural features and identifying areas for future model improvement.

**Key Features**

* **Data:** Utilizes a subset of the COCONUT database of natural products.
* **Descriptors:** Calculates molecular descriptors using RDKit.
* **Electronic Structure:** Employs xTB for efficient electronic structure calculations.
* **Workflows:** Uses CWL and Toil for reproducible and scalable computational workflows.
* **Machine Learning:** Implements Gradient Boosting Regression (GBR) and Multi-layer Perceptron Regression (MLPR) models.
* **Analysis:** Provides tools for analyzing feature importance and model performance.

**Installation**

Detailed installation instructions can be found in the :doc:`installation guide <installation_dev>` section of the documentation.

**Usage**

Examples and detailed usage instructions are available in the :doc:`usage Guide <usage>` section of the documentation. This includes information on using the command-line interface (CLI) and running the CWL workflows.

**Tutorial**
A comprehensive tutorial is available in the :doc:`tutorial <tutorial>` section of the documentation. This includes step-by-step instructions on how to set up and run the package, as well as examples of how to use the various features. Further material can be found in ``notebooks/`` directory.

**Contributing**

Contributions are welcome! Please see the :doc:`contributing <contributing>` guidelines for details on how to contribute to the project.

**Authors**

The authors and contributors are listed in the :doc:`authors <authors>` section of the documentation.

**License**

This project is licensed under the MIT License - see the :doc:`license <license>` file for details.

**Citation**

If you use this package in your research, please cite it as follows:

.. code-block:: bibtex

    @software{sthinius87_2025_15113790,
    author       = {sthinius87},
    title        = {sthinius87/HL-gaps-pub: HL-gaps v0.2.1},
    month        = mar,
    year         = 2025,
    publisher    = {Zenodo},
    version      = {v0.2.1},
    doi          = {10.5281/zenodo.15113790},
    url          = {https://doi.org/10.5281/zenodo.15113790},
    }
