"""Calculate AUC with and without Aerococcus for PTB data."""
import pickle

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score


def load_data(data_file, transf):
    """Load data."""
    with open(data_file, 'rb') as f:
        genus_data = pickle.load(f)

    # Get data corresponding to given transformation
    abun_data = genus_data['trans'][transf]
    meta_data = genus_data['meta']

    return abun_data, meta_data


def train_and_evaluate(abun_data, meta_data, train_study,
                       n_repeat, model_opts, taxa_to_drop=None):
    """Train n_repeat random forests and test on other datasets."""
    # List of test datasets
    test_list = [study for study in abun_data if study != train_study]

    # Create training data
    X_train = abun_data[train_study]
    if taxa_to_drop is not None and taxa_to_drop in X_train.columns:
        X_train = X_train.drop(columns=[taxa_to_drop])
    y_train = meta_data[train_study]['Preterm']

    # Create a dictionary to hold results
    AUC_test = dict()

    # Start with an empty list for each study
    for study in test_list:
        AUC_test[study] = []

    for i in range(n_repeat):

        # Train the model
        model = RandomForestClassifier(**model_opts)
        model.fit(X_train, y_train)

        # Test on other datasets
        for test_study in test_list:
            X_test = abun_data[test_study]
            if taxa_to_drop is not None and taxa_to_drop in X_test.columns:
                X_test = X_test.drop(columns=[taxa_to_drop])
            y_test = meta_data[test_study]['Preterm']

            # Calculate AUC and add to list
            y_pred_test = model.predict_proba(X_test)[:, 1]
            AUC_test[test_study].append(roc_auc_score(y_test, y_pred_test))

    avg_AUC = {study: np.mean(AUC_test[study]) for study in AUC_test}
    return avg_AUC


def main():
    """Calculate AUC with and without Aerococcus for PTB data."""
    abun_data, meta_data = load_data('data/genus_data.pkl', 'clr')

    # Set number of repetitions and random forest hyperparameters
    n_repeat = 10
    model_opts = {"n_estimators": 500,
                  "max_features": "auto",
                  "max_depth": None,
                  "max_samples": None,
                  "min_samples_split": 2}

    train_study = 'Ro'

    # Get average AUC with Aerococcus as a feature
    AUC_with_Aero = train_and_evaluate(abun_data,
                                       meta_data,
                                       train_study,
                                       n_repeat,
                                       model_opts,
                                       taxa_to_drop=None)

    # Get average AUC without Aerococcus as a feature
    AUC_without_Aero = train_and_evaluate(abun_data,
                                          meta_data,
                                          train_study,
                                          n_repeat,
                                          model_opts,
                                          taxa_to_drop='Aerococcus')

    # Create data frame with both sets of results
    AUC_compare = pd.DataFrame({'AUC_with_Aero': AUC_with_Aero,
                                'AUC_without_Aero': AUC_without_Aero})
    study_order = ['Br', 'Fe', 'Ki', 'St', 'Di', 'El',
                   'Bl', 'SC', 'Su', 'Ta', 'UC']
    AUC_compare = AUC_compare.loc[study_order, :]

    # Save results to csv
    AUC_compare.to_csv('AUC_Aerococcus.csv')


if __name__ == '__main__':
    main()
