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
X = pd.read_pickle("../calc_descriptors_final.pkl")
print(f"Number of features: {len(X.columns)}")

# Scale 'Ipc' feature
scaler = MinMaxScaler()
X["Ipc"] = scaler.fit_transform(X["Ipc"].values.reshape(-1, 1))

# Scale all features
scaler = StandardScaler().fit(X)
X_scaled = scaler.transform(X)

# Load target variable and split data
df = pd.read_pickle("../gap_smile.pkl")
y = df["GAP"].to_numpy()
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.3, random_state=42
)

# Define parameter grid for RandomizedSearchCV
params = {
    "learning_rate": uniform(0.02, 0.05),
    "loss": ["squared_error"],
    "max_depth": randint(8, 13),
    "n_estimators": randint(800, 1200),
    "subsample": uniform(0.2, 0.5),
}

# Perform RandomizedSearchCV
t0 = time.time()
model = GradientBoostingRegressor()
print("Starting randomized search")
search = RandomizedSearchCV(
    model,
    param_distributions=params,
    random_state=42,
    n_iter=2000,
    cv=2,
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
reg = GradientBoostingRegressor(**best_params)
reg.fit(X_train, y_train)

# Save model
with open("descr_GBR.pkl", "wb") as f:
    pickle.dump(reg, f)
dump(reg, "descr_GBR.joblib")

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
plt.savefig("PED_GBR_pred_on_trained_tl.png")
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
plt.savefig("PED_GBR_pred_on_untrained_tl.png")