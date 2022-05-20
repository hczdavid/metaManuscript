"""Calculate importance for PTB data."""
import pickle
import sklearn
import shap
import warnings
import argparse

import numpy as np

from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier

# Take data transformation as command-line argument
parser = argparse.ArgumentParser()
parser.add_argument('--trans',
                    default='prop',
                    choices=['prop', 'clr', 'clr4', 'log', 'rank',
                             'props', 'clrs', 'clr4s', 'logs', 'ranks'],
                    help='data transformation')
args = parser.parse_args()

# Load data
data_file = 'data/genus_data.pkl'

with open(data_file, 'rb') as f:
    genus_data = pickle.load(f)

# Get data corresponding to given transformation
abun_data = genus_data['trans'][args.trans]
meta_data = genus_data['meta']

# Use 5-fold cross-validation to calculate SHAP values for each dataset
# Repeat 10 times
n_repeat = 10
n_folds = 5
model_opts = {"n_estimators": 500,
              "max_features": "auto",
              "max_depth": None,
              "max_samples": None,
              "min_samples_split": 2}

importance_dict = dict()

for study in meta_data:
    print("Training on study {}".format(study))

    # Set train/validation data
    X_train_val = abun_data[study]
    y_train_val = meta_data[study]['Preterm']

    # List of SHAP values for every repetition
    shap_values = []

    for i in range(n_repeat):
        print("Repetition {}".format(i))

        # Split data
        kf = KFold(n_splits=n_folds, shuffle=True)

        for train_index, val_index in kf.split(X_train_val):
            X_train, X_val = (X_train_val.iloc[train_index],
                              X_train_val.iloc[val_index])
            y_train, y_val = (y_train_val.iloc[train_index],
                              y_train_val.iloc[val_index])

            # Don't train if only one label is present
            if len(set(y_train)) < 2:
                warnings.warn('Skipping fold - only one label')
                continue

            # Train model
            model = RandomForestClassifier(**model_opts)
            model.fit(X_train, y_train)

            # Calculate AUC if there are both labels in validation fold
            if len(set(y_val)) >= 2:
                y_pred = model.predict_proba(X_val)[:, 1]
                AUC = sklearn.metrics.roc_auc_score(y_val, y_pred)
            else:
                warnings.warn('Only one class in validation set. '
                              + 'AUC not calculated.')
                AUC = np.nan

            # Compute SHAP for validation data
            explainer = shap.TreeExplainer(model)

            iteration_dict = {'shap_values_val': explainer.shap_values(X_val),
                              'expected_value': explainer.expected_value,
                              'X_train': X_train,
                              'y_train': y_train,
                              'X_val': X_val,
                              'y_val': y_val,
                              'AUC': AUC}

            shap_values.append(iteration_dict)

    # Store results for this training study
    importance_dict[study] = shap_values
    print(" ")

# Save results
data_file = "importance_results_{}.pkl".format(args.trans)
with open(data_file, 'wb') as f:
    pickle.dump(importance_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
