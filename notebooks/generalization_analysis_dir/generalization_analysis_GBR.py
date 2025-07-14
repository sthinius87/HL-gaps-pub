#!/usr/bin/env python
# coding: utf-8

# # Robust Model Evaluation Script
#  This script addresses the reviewer's comment (R1, 2) by performing
#  a repeated random train-test split to evaluate model stability and generalization.

# ### import section

# In[2]:


import math
import sys
import re
import time
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.model_selection import ShuffleSplit
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import (
    r2_score,
    mean_squared_error,
    mean_absolute_error,
    root_mean_squared_error,
)
from typing import Dict

from joblib import load


# In[3]:


#get_ipython().system('pwd')


# ### loading and preprocessing

# In[4]:


print("--- Loading and Preprocessing Data ---")
try:
    # Load descriptors
    X = pd.read_pickle("../../data/processed/calc_descriptors_final.pkl")
    # Load target property (HL-Gap)
    df = pd.read_pickle("../../data/processed/gap_smile.pkl")
    y = df["GAP"].to_numpy()
    print("Data loaded successfully.")
except FileNotFoundError as e:
    print(f"Error: {e}. Make sure the data files are in the correct path.")
    exit()
    
# scale the data
scaler = MinMaxScaler()
X["Ipc"] = scaler.fit_transform(X["Ipc"].values.reshape(-1, 1))

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Load the trained model
try:
    gbr_model = load("../../data/processed/reg_GBR.joblib")
except FileNotFoundError:
    print("Error: '../../data/processed/reg_GBR.joblib' not found.")
    exit()


# ### metrics

# In[5]:


# --- Metrics ---
def get_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Calculates regression metrics."""
    metrics = {
        "R2": r2_score(y_true, y_pred),
        "MSE": mean_squared_error(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": root_mean_squared_error(y_true, y_pred),
    }
    return metrics


# ### setup splitting

# In[6]:


# Setup for 10 random splits (as requested by the reviewer)
n_splits = 10
test_size = 0.3
ss = ShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=42)


# ### evaluation

# In[ ]:


print(f"\n--- Performing evaluation over {n_splits} random splits ---")
sys.stdout.flush()
train_metrics_list = []
test_metrics_list = []
run_times = []
total_start_time = time.time()

for i, (train_index, test_index) in enumerate(ss.split(X_scaled)):
    split_start_time = time.time()
    print(f"Running split {i+1}/{n_splits}...")
    
    # Split data for the current run
    X_train, X_test = X_scaled[train_index], X_scaled[test_index]
    y_train, y_test = y[train_index], y[test_index]
    
    # Train the model from scratch on the current training fold
    gbr_model.fit(X_train, y_train)
    
    # Predict on both train and test sets
    y_pred_train = gbr_model.predict(X_train)
    y_pred_test = gbr_model.predict(X_test)
    
    # Calculate and store metrics
    train_metrics_list.append(get_metrics(y_train, y_pred_train))
    test_metrics_list.append(get_metrics(y_test, y_pred_test))
    
    # Calculate and report time for the split
    split_end_time = time.time()
    split_duration = split_end_time - split_start_time
    run_times.append(split_duration)
    print(f"--> Split {i+1} completed in {split_duration:.2f} seconds.")
sys.stdout.flush()

print("\nEvaluation complete.")


# ### aggregate and report results

# In[ ]:


print("\n--- Aggregated Performance Metrics (mean +/- std) ---")
sys.stdout.flush()
# Convert lists of dicts to DataFrames for easy calculation
df_train_metrics = pd.DataFrame(train_metrics_list)
df_test_metrics = pd.DataFrame(test_metrics_list)

# Calculate mean and std
train_mean = df_train_metrics.mean()
train_std = df_train_metrics.std()
test_mean = df_test_metrics.mean()
test_std = df_test_metrics.std()

# Print results
print("\n--- Training Set Performance ---")
for metric in train_mean.index:
    print(f"{metric}: {train_mean[metric]:.4f} +/- {train_std[metric]:.4f}")

print("\n--- Test Set Performance ---")
for metric in test_mean.index:
    print(f"{metric}: {test_mean[metric]:.4f} +/- {test_std[metric]:.4f}")

# Analyze the difference between train and test scores to discuss overfitting
print("\n--- Overfitting Analysis (Train R2 - Test R2) ---")
r2_diff = df_train_metrics['R2'] - df_test_metrics['R2']
print(f"Mean R2 Difference: {r2_diff.mean():.4f} +/- {r2_diff.std():.4f}")

print("\n--- Creating and Saving Results DataFrame ---")

# Create a dictionary to hold the final summary statistics
summary_data = {
    'Train_Mean': train_mean,
    'Train_Std': train_std,
    'Test_Mean': test_mean,
    'Test_Std': test_std
}
sys.stdout.flush()
# Create the results DataFrame
results_df = pd.DataFrame(summary_data)

# Add the R2 difference as a new row for completeness
r2_diff_mean = r2_diff.mean()
r2_diff_std = r2_diff.std()
results_df.loc['R2_Difference'] = [r2_diff_mean, r2_diff_std, np.nan, np.nan]


# Define the output path and filenames
pickle_output_filename = "../../data/processed/gbr_evaluation_metrics_summary.pkl"
csv_output_filename = "../../data/processed/gbr_evaluation_metrics_summary.csv"

try:
    # Save the DataFrame to a pickle file
    results_df.to_pickle(pickle_output_filename)
    print(f"\nResults DataFrame successfully saved to (pickle):\n{pickle_output_filename}")
    
    # Save the DataFrame to a CSV file
    results_df.to_csv(csv_output_filename)
    print(f"Results DataFrame successfully saved to (csv):\n{csv_output_filename}")
    
    # Display the final DataFrame
    print("\n--- Final Metrics Summary ---")
    print(results_df)
    
except Exception as e:
    print(f"Error saving file: {e}")
