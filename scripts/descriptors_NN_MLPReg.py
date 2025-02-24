#!/usr/bin/env python
# coding: utf-8


import time
import selfies as sf
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import uniform, randint

from sklearn.preprocessing import OneHotEncoder,StandardScaler, MinMaxScaler
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import PredictionErrorDisplay as PED
from sklearn.metrics import r2_score,mean_absolute_error, mean_squared_error
from sklearn.model_selection import cross_val_score,train_test_split, KFold,RandomizedSearchCV,GridSearchCV
import pickle
from joblib import dump, load

X=pd.read_pickle("../calc_descriptors_final.pkl")

print(X['Ipc'].max(),X['Ipc'].min())
scaler = MinMaxScaler()
X['Ipc'] = scaler.fit_transform(X['Ipc'].values.reshape(-1,1))
print(X['Ipc'].max(),X['Ipc'].min())
print("num_features: ",len(X.columns))

scaler=StandardScaler().fit(X)
scaled=scaler.transform(X)

df=pd.read_pickle("../gap_smile.pkl")
y=df["GAP"].to_numpy() #[0:100001:10]
X_train, X_test, y_train, y_test = train_test_split(scaled,y,test_size=0.3,random_state=42)


import inspect

def get_default_args(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }

get_default_args(MLPRegressor)

params={
 'hidden_layer_sizes': (100,),
 'activation': 'relu',
 'solver': 'adam',
 'alpha': 0.0001,
 'batch_size': 'auto',
 'learning_rate': 'constant',
 'learning_rate_init': 0.001,
 'power_t': 0.5,
 'max_iter': 1000,
 'shuffle': True,
 'random_state': 42,
 'tol': 1e-3,
 'verbose': True,
 'warm_start': False,
 'momentum': 0.9,
 'nesterovs_momentum': True,
 'early_stopping': False,
 'validation_fraction': 0.1,
 'beta_1': 0.9,
 'beta_2': 0.999,
 'epsilon': 1e-08,
 'n_iter_no_change': 10,
 'max_fun': 15000
}

params={'activation': 'relu', 'alpha': 0.00015215591406453634, 'batch_size': 'auto', 'beta_1': 0.8451262295455694, 'beta_2': 0.999, 'early_stopping': False, 'epsilon': 1e-08, 'hidden_layer_sizes': (120, 60, 20), 'learning_rate': 'constant', 'learning_rate_init': 0.001, 'max_fun': 15000, 'max_iter': 300, 'momentum': 0.9, 'n_iter_no_change': 5, 'nesterovs_momentum': True, 'power_t': 0.5, 'random_state': 42, 'shuffle': True, 'solver': 'adam', 'tol': 0.0005, 'validation_fraction': 0.1, 'verbose': True, 'warm_start': False}

t0=time.time()
reg=MLPRegressor(**params)
reg.fit(X_train, y_train)
t1=time.time()
print("time",t1-t0)
y_pred=reg.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
print("The mean squared error (MSE) on test set: {:.4f}".format(mse))


params={
#    'hidden_layer_sizes': (randint(1,11),),
    'hidden_layer_sizes': [
   #     (100,)      , (200,)    , #(300)   ,
   #     (20, 20)    , (40, 20)  , (80, 40),
        (200, 100),      #300
        (300, 150),      #450
        (200, 100, 50),  #350
        (300, 150, 50),  #500 
        (200,100,50,20), #370
        (300,200,100,50) #650 <-- added more neurons or layers
    ],
    'alpha': uniform(0.00006,0.00015),
#    'power_t': uniform(0.2,0.5),
    'beta_1': uniform(0.90,0.099),
#    'beta_2': uniform(0.990,0.009),
    'n_iter_no_change': [5],
    'tol': [5e-4],
    'max_iter': [500],
#    'verbose': [True],
}

t0=time.time()
model=MLPRegressor() #**params
print("starting randomized search")
search=RandomizedSearchCV(model,param_distributions=params,random_state=42,n_iter=4000,cv=4,n_jobs=12, verbose=4,return_train_score=True)
search.fit(X_train, y_train)
t1=time.time()
print('time:',t1-t0)
print("\n")


be=search.best_estimator_
print(be.get_params())

#X_train, X_test, y_train, y_test = train_test_split(scaled,y,test_size=0.33333,random_state=42)
best_params=search.best_estimator_.get_params()
reg=MLPRegressor(**best_params)
reg.fit(X_train, y_train)

#pickle
with open('descr_NN_MLP.pkl', 'wb') as f:
    pickle.dump(reg, f)
#joblib
dump(reg, 'descr_NN_MLP.joblib') 

y_pred_train=reg.predict(X_train)
mse = mean_squared_error(y_test, y_pred_train)
print("The mean squared error (MSE) on train set: {:.4f}".format(mse))
print("\n")
print("R2-score:",r2_score(y_pred=y_pred_train,y_true=y_train))

#prediction on trained molecules
fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
PED.from_predictions(
    y_train,
    y_pred=y_pred_train,
    kind="actual_vs_predicted",
    subsample=None,
    ax=axs[0],
    random_state=0,
    scatter_kwargs={"s":4.0}
)
axs[0].set_title("Actual vs. Predicted values")
axs[0].set_ylim(0,13)
axs[0].set_xlim(0,13)
axs[0].set_xticks(range(0,14,2))
axs[0].set_yticks(range(0,14,2))

PED.from_predictions(
    y_train,
    y_pred=y_pred_train,
    kind="residual_vs_predicted",
    subsample=None,
    ax=axs[1],
    random_state=0,
    scatter_kwargs={"s":4.0}
)
axs[1].set_title("Residuals vs. Predicted Values")
axs[1].set_xlim(0,13)
#fig.suptitle("Plotting cross-validated predictions")
#plt.title("pred vs. test")
#plt.show()

plt.savefig("PED_MLP_pred_on_trained.png")
plt.tight_layout()
plt.savefig("PED_MLP_pred_on_trained_tl.png")
plt.clf()


# In[ ]:


#prediction on untrained/test molecules
y_pred=reg.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
print("The mean squared error (MSE) on test set: {:.4f}".format(mse))
print("\n")
print("R2-score:",r2_score(y_pred=y_pred,y_true=y_test))

fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
plt.tight_layout()
PED.from_predictions(
    y_test,
    y_pred=y_pred,
    kind="actual_vs_predicted",
    subsample=None,
    ax=axs[0],
    random_state=0,
    scatter_kwargs={"s":4.0}
)
axs[0].set_title("Actual vs. Predicted values")
axs[0].set_ylim(0,13)
axs[0].set_xlim(0,13)
axs[0].set_xticks(range(0,14,2))
axs[0].set_yticks(range(0,14,2))

PED.from_predictions(
    y_test,
    y_pred=y_pred,
    kind="residual_vs_predicted",
    subsample=None,
    ax=axs[1],
    random_state=0,
    scatter_kwargs={"s":4.0}
)
axs[1].set_title("Residuals vs. Predicted Values")
axs[1].set_xlim(0,13)
axs[1].set_ylim(-4,4)
#fig.suptitle("Plotting cross-validated predictions")
#plt.title("pred vs. test")
#plt.show()
plt.savefig("PED_MLP_pred_on_untrained.png")
plt.tight_layout()
plt.savefig("PED_MLP_pred_on_untrained_tl.png")


