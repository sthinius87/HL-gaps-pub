.. _tutorial:

Tutorial: Predicting HOMO-LUMO Gaps
====================================

This tutorial guides you through the process of predicting HOMO-LUMO gaps using the `HL-gaps-pub` package. We'll cover the following steps:

1.  **Calculating Gaps:**
    *   Using the CWL workflow for large-scale calculations on a Slurm cluster.
    *   Using a Python script for smaller datasets or individual molecules.
2.  **Calculating and Selecting Descriptors:**
    *   Using RDKit to generate molecular descriptors.
    *   Applying feature selection techniques.
3.  **Machine Learning:**
    *   Training and optimizing a Gradient Boosting Regressor (GBR) model.
    *   Training and optimizing a Multi-layer Perceptron Regressor (MLPR) model.
    *   Analyzing model results and feature importance.
    *   Visualizing prediction errors and outliers.

.. note::

    This tutorial assumes you have already installed the `HL-gaps-pub` package and its dependencies.  See the :ref:`installation` section for installation instructions.  It also assumes you have downloaded the necessary data (COCONUT database) and placed it in the `data` directory, as described in the installation guide.

1. Calculating HOMO-LUMO Gaps
-----------------------------

This section demonstrates two methods for calculating HOMO-LUMO gaps: using a CWL workflow for high-throughput calculations and a Python script for smaller-scale calculations.

1.1. CWL Workflow (Slurm)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Common Workflow Language (CWL) allows for reproducible and scalable execution of computational workflows.  This is ideal for processing large numbers of molecules on a high-performance computing (HPC) cluster using a job scheduler like Slurm.

**1.1.1. Workflow Description (`scatter_wf.cwl`)**

The main workflow, `scatter_wf.cwl`, orchestrates the gap calculation process. It takes an array of database IDs as input and scatters the `gap.cwl` tool over these IDs.

.. literalinclude:: ../cwl_workflow/scatter_wf.cwl
   :language: yaml
   :caption: scatter_wf.cwl

**1.1.2. Tool Description (`gap.cwl`)**

The `gap.cwl` tool defines the command-line execution of the `HL-gaps-pub` CLI to calculate the gap for a *single* molecule.

.. literalinclude:: ../cwl_workflow/gap.cwl
   :language: yaml
   :caption: gap.cwl

**1.1.3. Input Parameters (`scatter_inp.yml`)**

The input parameters for the workflow are defined in a YAML file (`scatter_inp.yml`).  This file specifies the database path, filenames, calculation settings (number of conformers, accuracy, electronic temperature, xTB method), and the list of database IDs to process.

.. literalinclude:: ../cwl_workflow/scatter_inp.yml
    :language: yaml
    :caption: scatter_inp.yml
    :emphasize-lines: 12-16

.. note::
    The `emphasize-lines` directive above highlights lines 12-16, which contain the `inp_id_array`.  This shows where the user would modify the file to specify different database entries.

**1.1.4. Running the Workflow (`wf_run.sh`)**

The `wf_run.sh` script executes the CWL workflow using `toil-cwl-runner`. It includes options for both normal and restart modes. The restart mode is useful if the workflow was interrupted.

.. literalinclude:: ../cwl_workflow/wf_run.sh
   :language: bash
   :caption: wf_run.sh

**1.1.5. Generating the Input YAML (`write_scatter_yaml.py`)**

The `write_scatter_yaml.py` script provides a convenient way to generate the `scatter_inp.yml` file with a specified range of database IDs.  This is useful for automating the creation of input files for different subsets of the database.

.. literalinclude:: ../cwl_workflow/write_scatter_yaml.py
   :language: python
   :caption: write_scatter_yaml.py

**Example Usage:**

To generate an input file for IDs 407261 to 407269 (exclusive):

.. code-block:: console

    python scripts/write_scatter_yaml.py

This will create a `scatter_inp.yml` file in the current directory. Then, run the workflow:

.. code-block:: console

    bash cwl_workflow/wf_run.sh

1.2. Python Script (Smaller Datasets)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For smaller datasets or individual molecules, you can use the `serial_batch.py` script directly. This script reads an SDF file, calculates the HOMO-LUMO gap for each molecule, and writes the results to a file.

.. literalinclude:: ../scripts/serial_batch.py
   :language: python
   :caption: serial_batch.py

**Example Usage:**

.. code-block:: console

    python scripts/serial_batch.py

This will process molecules from IDs 407268 up to (but not including) 407270, as defined within the script.  You can modify the `id_start` and `id_end` variables in the script to process different ranges.

2. Calculating and Selecting Descriptors
----------------------------------------

This section demonstrates how to calculate molecular descriptors using RDKit and perform feature selection. The following notebook provides a complete example:

.. nbgallery::
    notebooks/calc_descriptors.ipynb

**Key Steps:**

1.  **Import Libraries:** Import necessary libraries (RDKit, pandas, scikit-learn).
2.  **Data Loading:** Load the `id_smile.pkl` file containing SMILES strings.
3.  **Descriptor Calculation:** Use RDKit's `MoleculeDescriptors` module to calculate a wide range of molecular descriptors.
4.  **Feature Selection:**
    *   Remove highly correlated features.
    *   Remove features with low variance.
    *   Handle missing values (NaN and infinity).

The resulting descriptors are saved to `calc_descriptors_final.pkl`.

3. Machine Learning
-------------------

This section covers training, optimizing, and analyzing machine learning models to predict HOMO-LUMO gaps.

3.1. Model Training and Hyperparameter Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We use `RandomizedSearchCV` from scikit-learn to optimize the hyperparameters of two models:

*   **Gradient Boosting Regressor (GBR):** (`hyperparameter_opt_GBReg.py`)
*   **Multi-layer Perceptron Regressor (MLPR):** (`hyperparameter_opt_NN_MLPReg.py`)

**3.1.1 GBR Optimization (`hyperparameter_opt_GBReg.py`):**

.. literalinclude:: ../scripts/hyperparameter_opt_GBReg.py
   :language: python
   :caption: hyperparameter_opt_GBReg.py

**Example Usage:**

.. code-block:: console

    python scripts/hyperparameter_opt_GBReg.py

This script performs a randomized search over the specified parameter space, trains the best model, saves the trained model (as `reg_GBR.joblib` and `reg_GBR.pkl`), and generates plots of the model's performance on training and test data.

**3.1.2. MLPR Optimization (`hyperparameter_opt_NN_MLPReg.py`):**

.. literalinclude:: ../scripts/hyperparameter_opt_NN_MLPReg.py
   :language: python
   :caption: hyperparameter_opt_NN_MLPReg.py

**Example Usage:**

.. code-block:: console

    python scripts/hyperparameter_opt_NN_MLPReg.py

This script is similar to the GBR optimization script but focuses on the MLPR model. It saves the trained model as `reg_NN_MLP.joblib` and `reg_NN_MLP.pkl` and generates performance plots.

3.2. Model Analysis and Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**3.2.1 General Analysis (GBR):**

.. nbgallery::
    notebooks/general_analysis_GBR.ipynb

This notebook loads and preprocesses data, then loads a pre-trained GBR model and evaluates general metrics as well as performing a visual analysis of different parameters used.

**3.2.2 General Analysis (MLPR):**

.. nbgallery::
    notebooks/general_analysis_MLPR.ipynb

This notebook loads and preprocesses data, then loads a pre-trained MLPR model and evaluates its performance using various metrics.

The notebook calculates and prints:

*   Overall R-squared, MSE, MAE, and RMSE on the training and test sets.
*   Metrics for subsets of the data based on specific descriptors.
*   Heatmaps of actual vs. predicted values and residuals vs. predicted values.
*   Metrics for data split by HOMO-LUMO gap values (< 6 eV and >= 6 eV).

**3.2.3. Feature Permutation Importance:**

The `feature_permutation.py` script calculates and displays the permutation importance of each feature, indicating its contribution to the model's predictive power.

.. literalinclude:: ../scripts/feature_permutation.py
    :language: python
    :caption: feature_permutation.py

The following notebook visualizes these importances.

.. nbgallery::
     notebooks/feature_permutation_show.ipynb

**3.2.4. Outlier Analysis:**

.. nbgallery::
     notebooks/heatmap_range.ipynb

This notebook analyzes the occurrence of prediction outliers within different ranges of HOMO-LUMO gap values and molecular subsets.

**3.2.5 Molecule Visualization:**

.. nbgallery::
    notebooks/save_images_of_molecules.ipynb

This notebook visualizes and saves images of molecules that exhibit significant prediction errors. It includes functions for displaying the molecules with the largest and smallest prediction errors and for saving outlier structures to a PDF file.