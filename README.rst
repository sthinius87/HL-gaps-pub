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

.. image:: https://img.shields.io/badge/DOI-PLACEHOLDER_DOI-blue.svg
    :target: https://doi.org/PLACEHOLDER_DOI
    :alt: DOI

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

.. image:: placeholder_for_graphical_abstract.png
    :alt: Graphical Abstract
    :width: 600px


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

Detailed installation instructions can be found in the :ref:`installation` section of the documentation.

**Usage**

Examples and detailed usage instructions are available in the :ref:`usage` section of the documentation. This includes information on using the command-line interface (CLI) and running the CWL workflows.

**Contributing**

Contributions are welcome! Please see the :ref:`contributing` guidelines for details on how to contribute to the project.

**Authors**

The authors and contributors are listed in the :ref:`authors` section of the documentation.

**License**

This project is licensed under the MIT License - see the `LICENSE` file for details.

**Citation**

If you use this package in your research, please cite it as follows:

.. code-block:: bibtex

    @misc{hl_gaps_pub_2025,
        author = {YOUR NAME(S)},
        title = {{HL-gaps-pub: Predicting HOMO-LUMO Gaps of Natural Compounds with Machine Learning}},
        year = {2025},
        publisher = {Zenodo},
        version = {v0.1.0},
        doi = {PLACEHOLDER_DOI},
        url = {https://zenodo.org/record/PLACEHOLDER_ZENODO_RECORD_ID}
    }

**Replace the following placeholders:**

* **`placeholder_for_graphical_abstract.png`:** The filename of your graphical abstract image. Place the image file in the same directory as your README.rst, or provide a relative path.
* **`main.yml` in the Build Status badge URL:** If your main GitHub Actions workflow file has a different name, update the URL accordingly.
* **`PLACEHOLDER_DOI`:** with the actual DOI you get when deposit your repository.
* **`https://zenodo.org/record/PLACEHOLDER_ZENODO_RECORD_ID`:** The link to your Zenodo deposit.
* **`YOUR NAME(S)`:** Replace this by the author names.
* **`year`:** The correct Year (Updated to 2025 based on current date).
* **`version`:** v0.2.0
  
This improved README provides a good starting point. It clearly explains the project, highlights its key features, and directs users to the detailed documentation for installation and usage instructions. It also includes placeholders for important information like the graphical abstract and citation information. The use of reStructuredText directives (`.. image::`, `.. _HL-gaps-pub:`, `:ref:`) ensures proper formatting and linking within the Sphinx documentation. The inclusion of badges gives a quick overview of project health.




