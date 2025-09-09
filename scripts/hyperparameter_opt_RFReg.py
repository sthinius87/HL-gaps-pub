import pickle
import time
from joblib import dump, load
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import randint
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import (
    PredictionErrorDisplay as PED,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)
from sklearn.model_selection import (
    RandomizedSearchCV,
    train_test_split,
)
from sklearn.preprocessing import MinMaxScaler, StandardScaler


# --- Robust Path Resolution ---
# Get the absolute path to the directory where this script is located
SCRIPT_DIR = Path(__file__).resolve().parent

# Define the path to your data files relative to the script's directory
X_FILE_PATH = SCRIPT_DIR / "../data/processed/calc_descriptors_final.pkl"
Y_FILE_PATH = SCRIPT_DIR / "../data/processed/gap_smile.pkl"
# -----------------------------

# Load data
X = pd.read_pickle(X_FILE_PATH)
print(f"Number of features: {len(X.columns)}")

# Scale 'Ipc' feature
scaler = MinMaxScaler()
X["Ipc"] = scaler.fit_transform(X["Ipc"].values.reshape(-1, 1))

# Scale all features
scaler = StandardScaler().fit(X)
X_scaled = scaler.transform(X)

# Load target variable and split data
df = pd.read_pickle(Y_FILE_PATH)
y = df["GAP"].to_numpy()
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.3, random_state=42
)

# Define parameter grid for RandomizedSearchCV
# --- THIS IS THE KEY MODIFIED SECTION FOR RF ---
params = {
    "n_estimators": randint(200, 1000),
    "max_depth": randint(10, 50),
    "min_samples_split": randint(2, 20),
    "min_samples_leaf": randint(1, 20),
    "max_features": ['sqrt', 'log2', 1.0],
    "bootstrap": [True, False],
}

# Perform RandomizedSearchCV
t0 = time.time()
model = RandomForestRegressor(random_state=42) # Instantiate the RF model
print("Starting randomized search for Random Forest")
search = RandomizedSearchCV(
    model,
    param_distributions=params,
    random_state=42,
    n_iter=2000, 
    cv=3,
    n_jobs=-1,
    verbose=4,
    return_train_score=True,
)
search.fit(X_train, y_train)
t1 = time.time()
print(f"Time: {t1 - t0}")
print("\n")

# Train best model
best_params = search.best_params_
print("Best parameters found: ", best_params)
reg = RandomForestRegressor(**best_params, random_state=42)
reg.fit(X_train, y_train)

# Save model
with open("reg_RF.pkl", "wb") as f:
    pickle.dump(reg, f)
dump(reg, "reg_RF.joblib")

# Evaluate on training data
y_pred_train = reg.predict(X_train)
mse_train = mean_squared_error(y_train, y_pred_train)
r2_train = r2_score(y_pred=y_pred_train, y_true=y_train)

print(f"Train MSE: {mse_train:.4f}")
print(f"Train R2-score: {r2_train}")

# Plot predictions on training data
fig, axs = plt.subplots(ncols=2, figsize=(8, 4))

PED.from_predictions(
    y_train,
    y_pred=y_pred_train,
    kind="actual_vs_predicted",
    ax=axs[0],
    scatter_kwargs={"s": 4.0},
)
axs[0].set_title("Actual vs. Predicted values (Train)")
axs[0].set_ylim(0, 13)
axs[0].set_xlim(0, 13)
axs[0].set_xticks(range(0, 14, 2))
axs[0].set_yticks(range(0, 14, 2))

PED.from_predictions(
    y_train,
    y_pred=y_pred_train,
    kind="residual_vs_predicted",
    ax=axs[1],
    scatter_kwargs={"s": 4.0},
)
axs[1].set_title("Residuals vs. Predicted Values (Train)")
axs[1].set_xlim(0, 13)

plt.tight_layout()
plt.savefig("PED_RF_pred_on_trained_tl.png")
plt.clf()

# Evaluate on test data
y_pred_test = reg.predict(X_test)
mse_test = mean_squared_error(y_test, y_pred_test)
r2_test = r2_score(y_pred=y_pred_test, y_true=y_test)

print(f"Test MSE: {mse_test:.4f}")
print(f"Test R2-score: {r2_test}")

# Plot predictions on test data
fig, axs = plt.subplots(ncols=2, figsize=(8, 4))

PED.from_predictions(
    y_test,
    y_pred=y_pred_test,
    kind="actual_vs_predicted",
    ax=axs[0],
    scatter_kwargs={"s": 4.0},
)
axs[0].set_title("Actual vs. Predicted values (Test)")
axs[0].set_ylim(0, 13)
axs[0].set_xlim(0, 13)
axs[0].set_xticks(range(0, 14, 2))
axs[0].set_yticks(range(0, 14, 2))

PED.from_predictions(
    y_test,
    y_pred=y_pred_test,
    kind="residual_vs_predicted",
    ax=axs[1],
    scatter_kwargs={"s": 4.0},
)
axs[1].set_title("Residuals vs. Predicted Values (Test)")
axs[1].set_xlim(0, 13)
axs[1].set_ylim(-4, 4)

plt.tight_layout()
plt.savefig("PED_RF_pred_on_untrained_tl.png")
