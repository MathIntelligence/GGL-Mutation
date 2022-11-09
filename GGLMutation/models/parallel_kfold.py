#!/usr/bin/env python
"""Cross-validation experiments for Mutation Data

# Authors: Masud Rana, Duc Nguyen (@NguyenLab)

"""

import os
from sklearn import ensemble
from sklearn.model_selection import KFold, RepeatedKFold
import pandas as pd
from sklearn.metrics import mean_squared_error
from scipy import stats
import math
import numpy as np
import ntpath
import pandas as pd
import argparse
from random import randint
import warnings
import time
warnings.filterwarnings('ignore')
import catboost as ctb
import optuna
from sklearn.model_selection import cross_val_predict
import pickle


def combine_aux_feature_files(feature_file, aux_feature_file):
    """ Do kfold and save results
    Parameters
    ----------
    feature_file: str
    
    aux_feature_file: str
    
    Return
    -------
    pandas.DataFrame
    """
    
    df1 = pd.read_csv(feature_file)
    df2 = pd.read_csv(aux_feature_file, header=None)
    df1 = pd.concat([df1, df2], axis=1)
    
    return df1


def do_kfold(df, random_state=123, kfold_random_state=1, ML_type='gbt', n_splits=10):
    """
    Parameters
    ----------
    df: pandas.DataFrame
        Feature of a model in DataFrame
        
    random_state: int, default=123
        Machine learning model random state
        
    kfold_random_state: int, default=1
        Seed for generating different folds
        
    ML_type: {'gbt'}, default='gbt'
        
    n_splits : int, default=10
        Number of folds. Must be at least 2.
        

    Return
    -------
    output: dictionary
            {'y_pred': y_pred, 'y_ref': y, 'test_indices': test_indices,
              'random_state': random_state, 'kfold_random_state': kfold_random_state}
    """
    
    kf = KFold(
        n_splits=n_splits, shuffle=True,
        random_state=kfold_random_state)

    X = df.drop(['PDBID', 'ddG'], axis=1)
    X = X.values
    y = df[['ddG']]
    y = y.values.ravel()

    if ML_type == 'gbt':
        params = {
            'n_estimators': 40000,
            'max_depth': 6,
            'min_samples_split': 3,
            'learning_rate': 0.001,
            'loss': 'squared_error',
            'subsample': 0.7,
            'max_features': 'sqrt',
            'random_state': random_state
        }

        clf = ensemble.GradientBoostingRegressor(**params)
    elif ML_type == 'xgb':
        params = {
            'n_estimators': 10000,
            'max_depth': 7,
            'min_child_weight': 1,
            'min_samples_split': 3,
            'learning_rate': 0.005,
            'gamma':  0.2,
            'subsample': 0.1,
            'objective': 'reg:linear',
            'eval_metric': 'rmse',
            'colsample_bytree': 0.1,
            'n_jobs': 1,
            'max_features': 'sqrt',
            'random_state': random_state
        }
        clf = XGBRegressor(**params)
    else:
        params = {
            'n_estimators': 1000,
            'max_features': 'auto',
            'n_jobs': 5,
            'random_state': random_state
        }
        clf = ensemble.RandomForestRegressor(**params)

    tau = []
    Rp = []
    RMSE = []
    y_pred = cross_val_predict(clf, X, y, cv=kf, n_jobs=n_splits, verbose=10)
    test_indices = []
    for train_index, test_index in kf.split(y):
        test_indices.append(test_index)
        
    output = {'y_pred': y_pred, 'y_ref': y, 'test_indices': test_indices,
              'random_state': random_state, 'kfold_random_state': kfold_random_state}
    return output


def main(feature_file, aux_feature_file, output_folder,
             ML_type='gbt', random_state=123, kfold_random_state=1, n_splits=10):
    """ Do kfold and save results
    Parameters
    ----------
    feature_file: str
        feature of a model, often in *.csv
        
    aux_feature_file: str
        auxiliary features, oftein in *.csv
        
    output_folder: str
        folder to save the result
        
    ML_type: {'gt'}, default = 'gbt'       
        gradient boosting tree from sklearn: ML_type = 'gbt'
        
    random_state: int , default=123
    	machine learning random state      
        
    kfold_random_state: int, default=1
        seed for splitting folds   
        
    n_splits: int, default=10
        Number of splits in kfold
        
    Return
    -------
    save results in pickle format to <output_folder>
    """
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    file_name = ntpath.basename(feature_file).split('.')[0]
    aux_file_name = ntpath.basename(aux_feature_file).split('.')[0]
    
 
    df = combine_aux_feature_files(feature_file, aux_feature_file)
    
    output = do_kfold(df, random_state, kfold_random_state, ML_type, n_splits)
    
   
    output_predict_file = f'{output_folder}/{file_name}.{ML_type}.{n_splits}splits.{kfold_random_state}.{random_state}.pkl'
    with open(output_predict_file, 'wb') as handle:
        pickle.dump(output, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Do KFold experiments",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', '--ml-type', help='machine learning type, either gbt, or rf, or xgbt',
                        type=str, default='gbt')
    
    parser.add_argument('-o', '--output-folder', type=str,
                        help='folder to save kfold results')                       

    parser.add_argument('-f', '--feature-file', type=str,
                        help='feature file')
    parser.add_argument('-a', '--aux-feature-file', type=str,
                        help='auxiliary feature file')
    
    parser.add_argument('-r', '--random-state', type=int, 
                        help='ML random state')
    parser.add_argument('-k', '--kfold-random-state', type=int, 
                        help='Kfold random state')
    
    parser.add_argument('-n', '--n-splits', type=int, default=10,
                        help='Number of splits')
    
    args = parser.parse_args()
    
  
    feature_file = args.feature_file
    ML_type = args.ml_type
    output_folder = args.output_folder
    aux_feature_file = args.aux_feature_file

    if args.random_state:
        random_state = args.random_state
    else:
        random_state = randint(1,1000000)

    if args.kfold_random_state:
        kfold_random_state = args.kfold_random_state
    else:
        kfold_random_state = randint(1,1000000)
        
    n_splits = args.n_splits

    main(feature_file, aux_feature_file, output_folder,
         ML_type, random_state, kfold_random_state, n_splits)
