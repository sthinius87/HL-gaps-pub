import numpy as np
import pandas as pd

from joblib import dump, load
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import MinMaxScaler, StandardScaler

def calculate_and_display_permutation_importance(
    model: MLPRegressor, X: pd.DataFrame, y: np.ndarray, dataset_name: str
) -> np.ndarray:
    """
    Calculates and displays permutation importance for a given model and dataset.

    Args:
        model: The trained machine learning model.
        X: The feature matrix.
        y: The target vector.
        dataset_name: The name of the dataset ("train" or "test").

    Returns:
        The sorted indices of feature importances.
    """

    result = permutation_importance(
        model, X, y, scoring="r2", n_repeats=12, random_state=42, n_jobs=-1
    )
    sorted_idx = result.importances_mean.argsort()[::-1]

    print(f"{dataset_name} set:\n")
    for i in sorted_idx:
        print(
            f"{X.columns[i]}: {result.importances_mean[i]:.4f} +/- {result.importances_std[i]:.4f}"
        )

    dump(result, f"PI_{dataset_name}.joblib")
    return sorted_idx


def main():
    """
    Loads data, trains a model, and calculates permutation importance.
    """
    X: pd.DataFrame = pd.read_pickle("../data/processed/calc_descriptors_final.pkl")
    df: pd.DataFrame = pd.read_pickle("../data/processed/gap_smile.pkl")
    y: np.ndarray = df["GAP"].to_numpy()

    # Scale data
    ipc_scaler = MinMaxScaler()
    X["Ipc"] = ipc_scaler.fit_transform(X["Ipc"].values.reshape(-1, 1))

    X_scaler = StandardScaler()
    X_scaled: np.ndarray = X_scaler.fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, test_size=0.3, random_state=42
    )

    # Train the model
    model_file = "../data/processed/reg_NN_MLP.joblib"
    mlpr: MLPRegressor = load(model_file)
    mlpr.fit(X_train, y_train)

    # Calculate and display permutation importance
    calculate_and_display_permutation_importance(mlpr, X_train, y_train, "train")
    calculate_and_display_permutation_importance(mlpr, X_test, y_test, "test")


if __name__ == "__main__":
    main()