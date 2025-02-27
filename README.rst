===========
HL-gaps-pub
===========

High throughput tight binding calculation of electronic HOMO-LUMO gaps and its prediction for natural compounds

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/badges/master/pipeline.svg
   :alt:

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/-/badges/release.svg
   :alt:

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/-/jobs/artifacts/master/raw/badges/wheel.svg?job=publish_badges
   :width: 16%
   :target: https://HL_gaps_pub/ifam418/HL-gaps-pub/-/packages
   :alt:

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/-/jobs/artifacts/master/raw/badges/dockerimage.svg?job=publish_badges
   :width: 16%
   :target: https://HL_gaps_pub/ifam418/HL-gaps-pub/container_registry
   :alt:

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/-/jobs/artifacts/master/raw/badges/tox.svg?job=integration_tests
   :alt:

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/-/jobs/artifacts/master/raw/badges/docs.svg?job=pages
   :alt:

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/-/jobs/artifacts/master/raw/badges/pylint.svg?job=pylint
   :alt:

.. image:: https://HL_gaps_pub/ifam418/HL-gaps-pub/badges/master/coverage.svg
   :alt:

.. image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
   :target: https://pycqa.github.io/isort/

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
  :target: https://github.com/ambv/black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt:






.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT


.. image:: placeholder_for_graphical_abstract.png
   :alt: Graphical Abstract
   :width: 600px

.. _HL-gaps-pub:

HL-gaps-pub
===========

|Build Status| |Coverage Status| |License| |DOI|

.. |Build Status| image:: https://img.shields.io/gitlab/pipeline/USER/HL-gaps-pub/main?branch=main
   :target: https://gitlab.com/USER/HL-gaps-pub/-/pipelines
   :alt: Build Status
.. |Coverage Status| image:: https://img.shields.io/codecov/c/gitlab/USER/HL-gaps-pub.svg
   :target: https://codecov.io/gl/USER/HL-gaps-pub
   :alt: Coverage Status
.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: LICENSE
   :alt: License
.. |DOI| image:: https://img.shields.io/badge/DOI-PLACEHOLDER_DOI-blue.svg
   :target: https://doi.org/PLACEHOLDER_DOI
   :alt: DOI

**Predicting HOMO-LUMO Gaps of Natural Compounds with Machine Learning**

This package provides tools and workflows for predicting the HOMO-LUMO (HL) gap of natural compounds using machine learning models.  It leverages the COCONUT database, RDKit for descriptor calculation, xTB for electronic structure calculations, and Common Workflow Language (CWL) for reproducible computational workflows.  The models used are a Gradient Boosting Regressor (GBR) and a Multi-layer Perceptron Regressor (MLPR).

**Abstract**

This research investigates predicting the HOMO-LUMO (HL) gap of natural compounds, a crucial property for understanding molecular electronic behavior relevant to cheminformatics and material science. Addressing the computational expense of traditional methods, this study develops a high-throughput, machine learning-based approach. Using 407,000 molecules from the COCONUT database, RDKit was employed to calculate and select molecular descriptors. The computational workflow, managed by Toil and CWL on a high-performance computing Slurm cluster, utilized xTB for electronic structure calculations with Boltzmann weighting across multiple conformational states. Gradient boosting regression (GBR) and a Multi-layer Perceptron regressor (MLPR) were compared based on their ability to accurately predict HL-gaps in this chemical space. Key findings reveal molecular polarizability, particularly SMR_VSA descriptors, as crucial for HL-gap determination in both models. Aromatic rings and functional groups, such as ketones, also significantly influence the HL-gap prediction. While the MLPR model demonstrated good overall predictive performance, accuracy varied across molecular subsets. Challenges were observed in predicting HL-gaps for molecules containing aliphatic carboxylic acids, alcohols, and amines in molecular systems with complex electronic structure. This work emphasizes the importance of polarizability and structural features in HL-gap predictive modeling, showcasing the potential of machine learning while also highlighting limitations in handling specific structural motifs. These limitations point towards promising perspectives for further model improvements.

**Key Features**

*   **Data:** Utilizes a subset of the COCONUT database of natural products.
*   **Descriptors:** Calculates molecular descriptors using RDKit.
*   **Electronic Structure:**  Employs xTB for efficient electronic structure calculations.
*   **Workflows:**  Uses CWL and Toil for reproducible and scalable computational workflows.
*   **Machine Learning:**  Implements Gradient Boosting Regression (GBR) and Multi-layer Perceptron Regression (MLPR) models.
*   **Analysis:** Provides tools for analyzing feature importance and model performance.

**Installation**

Detailed installation instructions can be found in the :ref:`installation` section of the documentation.

**Usage**

Examples and detailed usage instructions are available in the :ref:`usage` section of the documentation. This includes information on using the command-line interface (CLI) and running the CWL workflows.

**Contributing**

Contributions are welcome!  Please see the :ref:`contributing` guidelines for details on how to contribute to the project.

**Authors**

The authors and contributors are listed in the :ref:`authors` section of the documentation.

**License**

This project is licensed under the MIT License - see the `LICENSE` file for details.

**Citation**

If you use this package in your research, please cite it as follows:

.. code-block:: bibtex

    @misc{hl_gaps_pub,
      author = {YOUR NAME(S)},
      title = {{HL-gaps-pub: Predicting HOMO-LUMO Gaps of Natural Compounds with Machine Learning}},
      year = {2024}, %Update the Year
      publisher = {Zenodo},
      version = {v0.1.0}, % Update as you go
      doi = {PLACEHOLDER_DOI},
      url = {https://zenodo.org/record/PLACEHOLDER_ZENODO_RECORD_ID}
    }


**Replace the following placeholders:**

*   **`placeholder_for_graphical_abstract.png`:**  The filename of your graphical abstract image.  Place the image file in the same directory as your README.rst, or provide a relative path.
*   **`USER/HL-gaps-pub`:** in the build and coverage badges with your actual Gitlab username and project path.
*  **`PLACEHOLDER_DOI`:** with the actual DOI you get when deposit your repository.
*  **`https://zenodo.org/record/PLACEHOLDER_ZENODO_RECORD_ID`:** The link to your Zenodo deposit.
* **`YOUR NAME(S)`**: Replace this by the author names
* **`year`**: The correct Year.
* **`version`**: The correct version.

This improved README provides a good starting point. It clearly explains the project, highlights its key features, and directs users to the detailed documentation for installation and usage instructions. It also includes placeholders for important information like the graphical abstract and citation information. The use of reStructuredText directives (`.. image::`, `.. _HL-gaps-pub:`, `:ref:`) ensures proper formatting and linking within the Sphinx documentation. The inclusion of badges gives a quick overview of project health.




