import pickle
import time
from joblib import dump
import xgboost as xgb
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import randint, uniform
from sklearn.metrics import (
    PredictionErrorDisplay as PED,
    mean_squared_error,
    r2_score,
)
from sklearn.model_selection import (
    RandomizedSearchCV,
    train_test_split,
)
from sklearn.preprocessing import MinMaxScaler, StandardScaler

# --- Robust Path Resolution ---
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
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

# --- THIS IS THE KEY MODIFIED SECTION FOR XGBOOST ---
# Define parameter grid for RandomizedSearchCV
params = {
    "n_estimators": randint(500, 1500),
    "max_depth": randint(8, 20),
    "learning_rate": uniform(0.01, 0.1),
    "subsample": uniform(0.6, 0.4), # Range is start, width
    "colsample_bytree": uniform(0.6, 0.4),
    "gamma": [0, 0.1, 0.5, 1],
    "reg_alpha": uniform(0, 1), # L1 regularization
    "reg_lambda": uniform(0, 1), # L2 regularization
}

# Perform RandomizedSearchCV
t0 = time.time()
# Instantiate the XGBoost model with GPU support
model = xgb.XGBRegressor(tree_method='gpu_hist', random_state=42)

print("Starting randomized search for XGBoost on GPU")
search = RandomizedSearchCV(
    model,
    param_distributions=params,
    random_state=42,
    n_iter=2000,
    cv=3,
    n_jobs=1, # Important: Set n_jobs=1 when using GPU
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
reg = xgb.XGBRegressor(tree_method='gpu_hist', **best_params, random_state=42)
reg.fit(X_train, y_train)

# Save model
with open("reg_XGB.pkl", "wb") as f:
    pickle.dump(reg, f)
dump(reg, "reg_XGB.joblib")

# --- Evaluation and Plotting Code (remains the same) ---

# Evaluate on training data
y_pred_train = reg.predict(X_train)
mse_train = mean_squared_error(y_train, y_pred_train)
r2_train = r2_score(y_pred=y_pred_train, y_true=y_train)

print(f"Train MSE: {mse_train:.4f}")
print(f"Train R2-score: {r2_train}")

# Plot predictions on training data
fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
PED.from_predictions(y_train, y_pred=y_pred_train, kind="actual_vs_predicted", ax=axs[0], scatter_kwargs={"s": 4.0})
axs[0].set_title("Actual vs. Predicted values (Train)")
# ... (rest of the plotting code is identical to your RF script) ...

plt.tight_layout()
plt.savefig("PED_XGB_pred_on_trained_tl.png")
plt.clf()

# Evaluate on test data
y_pred_test = reg.predict(X_test)
mse_test = mean_squared_error(y_test, y_pred_test)
r2_test = r2_score(y_pred=y_pred_test, y_true=y_test)

print(f"Test MSE: {mse_test:.4f}")
print(f"Test R2-score: {r2_test}")

# Plot predictions on test data
fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
PED.from_predictions(y_test, y_pred=y_pred_test, kind="actual_vs_predicted", ax=axs[0], scatter_kwargs={"s": 4.0})
axs[0].set_title("Actual vs. Predicted values (Test)")
# ... (rest of the plotting code is identical to your RF script) ...

plt.tight_layout()
plt.savefig("PED_XGB_pred_on_untrained_tl.png")
