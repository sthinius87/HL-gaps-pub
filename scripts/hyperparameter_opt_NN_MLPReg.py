import pickle
import time
from joblib import dump, load

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import randint, uniform
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    PredictionErrorDisplay as PED,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)
from sklearn.model_selection import (
    KFold,
    GridSearchCV,
    RandomizedSearchCV,
    cross_val_score,
    train_test_split,
)
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder, StandardScaler

# Load data
X = pd.read_pickle("../data/processed/calc_descriptors_final.pkl")
y_df = pd.read_pickle("../data/processed/gap_smile.pkl")
y = y_df["GAP"].to_numpy()

# Scale features
scaler = MinMaxScaler()
X["Ipc"] = scaler.fit_transform(X["Ipc"].values.reshape(-1, 1))
print(f"Number of features: {len(X.columns)}")

scaler = StandardScaler().fit(X)
X_scaled = scaler.transform(X)

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.3, random_state=42
)

# Define parameter grid for RandomizedSearchCV
params = {
    "hidden_layer_sizes": [
        (200, 100),
        (300, 150),
        (200, 100, 50),
        (300, 150, 50),
        (200, 100, 50, 20),
        (300, 200, 100, 50),
    ],
    "alpha": uniform(0.00006, 0.00015),
    "beta_1": uniform(0.90, 0.099),
    "n_iter_no_change": [5],
    "tol": [5e-4],
    "max_iter": [500],
}

# Perform RandomizedSearchCV
t0 = time.time()
model = MLPRegressor()
print("Starting randomized search")
search = RandomizedSearchCV(
    model,
    param_distributions=params,
    random_state=42,
    n_iter=4000,
    cv=4,
    n_jobs=12,
    verbose=4,
    return_train_score=True,
)
search.fit(X_train, y_train)
t1 = time.time()
print(f"Time: {t1 - t0}")
print("\n")

# Train best model
best_params = search.best_estimator_.get_params()
reg = MLPRegressor(**best_params)
reg.fit(X_train, y_train)

# Save model
with open("reg_NN_MLP.pkl", "wb") as f:
    pickle.dump(reg, f)
dump(reg, "reg_NN_MLP.joblib")

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
plt.savefig("PED_MLP_pred_on_trained.png")
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
plt.savefig("PED_MLP_pred_on_untrained_.png")