#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
from sklearn.inspection import permutation_importance
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import OneHotEncoder,StandardScaler, MinMaxScaler
from joblib import dump, load

X=pd.read_pickle("../calc_descriptors_final.pkl")


df=pd.read_pickle("../gap_smile.pkl")
y=df["GAP"].to_numpy() #[0:100001:10]

# Scale data
scaler = MinMaxScaler()
X['Ipc'] = scaler.fit_transform(X['Ipc'].values.reshape(-1,1))
scaler=StandardScaler().fit(X)
scaled=scaler.transform(X)
X_train, X_test, y_train, y_test = train_test_split(scaled,y,test_size=0.3,random_state=42)

# Train your MLPR model
load_file='descr_NN_MLP.joblib'
mlpr=load(load_file)
#mlpr = MLPRegressor(random_state=42) #Example parameters
mlpr.fit(X_train, y_train)

# Calculate permutation importance train
result_tr = permutation_importance(mlpr, X_train, y_train, scoring='r2', n_repeats=12, random_state=42, n_jobs=-1)
sorted_idx = result_tr.importances_mean.argsort()[::-1]

# Display feature importances train
print("train set: \n")
for i in sorted_idx:
    print(f"{X.columns[i]}: {result_tr.importances_mean[i]:.4f} +/- {result_tr.importances_std[i]:.4f}")
#joblib
dump(result_tr, 'PI_tr.joblib') 

# Calculate permutation importance test
result_te = permutation_importance(mlpr, X_test, y_test, scoring='r2', n_repeats=12, random_state=42, n_jobs=-1)
sorted_idx = result_te.importances_mean.argsort()[::-1]

# Display feature importances test
print("test set: \n")
for i in sorted_idx:
    print(f"{X.columns[i]}: {result_te.importances_mean[i]:.4f} +/- {result_te.importances_std[i]:.4f}")

#joblib
dump(result_te, 'PI_te.joblib') 

