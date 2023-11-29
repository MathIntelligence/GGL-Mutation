#!/usr/bin/env python
# Author: Masud Rana <masud.rana@uky.edu>
# Last modified: Nov 28, 2023

import os
from sklearn import ensemble
from sklearn.model_selection import KFold
# from sklearn.model_selection import StratifiedKFold
import pandas as pd
from sklearn.metrics import mean_squared_error
from scipy import stats
import math
import numpy as np
import ntpath
import pandas as pd
import argparse
from random import randint


def do_Kfold(feature_df, ML_type='gbt', ml_random_state=123, kf_random_state=1, n_splits=10):
   """
    Parameters
    ----------
    df: pandas.DataFrame
        Feature of a model in DataFrame
    ML_type: str, default='gbt'
        either 'gbt' or 'rf'  
    ml_random_state: int, default=123
        Machine learning model random state    
    kf_random_state: int, default=1
        Kfold random state 
    n_splits : int, default=10
        Number of folds. Must be at least 2.
        

    Output
    -------
    Rp, tau, RMSE, Kfold predictions
    """
    
    kfold_random_state = kf_random_state
    kf = KFold(
        n_splits=n_splits, shuffle=True,
        random_state=kfold_random_state)

##############################################################
    if ML_type == 'gbt':
        # taken from ECIF
        params = {
            'n_estimators': 40000,
            'max_depth': 6,
            'min_samples_split': 3,
            'learning_rate': 0.001,
            'loss': 'ls',
            'subsample': 0.7,
            'max_features': 'sqrt',
            'random_state': ml_random_state
        }
            #old
#         params = {
#             'n_estimators': 10000,
#             'max_depth': 7,
#             'min_samples_split': 3,
#             'learning_rate': 0.01,
#             'loss': 'ls',
#             'subsample': 0.3,
#             'max_features': 'sqrt',
#             'random_state': random_state
#         }

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
#####################################################################################

    # df = feature_df.drop(['ID', 'ddG'], axis=1)
    df = feature_df

    #Separate direct and reverse mutations in two data frames
    rev_mut_size = int(df.shape[0]/2)
    df_dir = df.iloc[:rev_mut_size,:].reset_index(drop=True) # direct mutations
    df_rev = df.iloc[rev_mut_size:,:].reset_index(drop=True) # reverse mutations


    tau = []
    Rp = []
    RMSE = []

    ite = 1

    for train_index, test_index in kf.split(df_dir):
        df_dir_train, df_dir_test = df_dir.loc[train_index], df_dir.loc[test_index]
        df_rev_train, df_rev_test = df_rev.loc[train_index], df_rev.loc[test_index]

        df_train = pd.concat([df_dir_train, df_rev_train], axis=0).reset_index(drop=True)
        df_test = pd.concat([df_dir_test, df_rev_test], axis=0).reset_index(drop=True)

        X_train = df_train.drop(['ID', 'ddG'], axis=1)
        X_train = X_train.values

        X_test = df_test.drop(['ID', 'ddG'], axis=1)
        X_test = X_test.values

        y_train = df_train['ddG'].values
        y_test = df_test['ddG'].values

        clf.fit(X_train, y_train)
        y_prd = clf.predict(X_test)

        [tau_, pval] = stats.kendalltau(y_test, y_prd)
        tau.append(tau_)
        [Rp_, pval] = stats.pearsonr(y_test, y_prd)
        Rp.append(Rp_)
        RMSE.append(math.sqrt(mean_squared_error(y_test, y_prd)))

        if ite==1:
            kfold_df = pd.DataFrame({'ID':df_test['ID'], 'Predicted ddG':y_prd, 'Exp ddG':y_test})
        else:
            dummy_df = pd.DataFrame({'ID':df_test['ID'], 'Predicted ddG':y_prd, 'Exp ddG':y_test})
            kfold_df = pd.concat([kfold_df,dummy_df], axis=0).reset_index(drop=True)
        ite +=1

    return np.mean(np.array(Rp)), np.mean(np.array(tau)),\
        np.mean(np.array(RMSE)), kfold_df



def main(feature_file, ML_type, ml_random_state, kf_random_state, output_folder):
    """ Do kfold and save results
    Parameters
    ----------
    feature_file: str
        feature for dataset 1, often in *.csv
    ML_type: str
        either 'gbt' or 'rf'          
    output_folder: str
        folder to save the result
    Output
    -------
    csv file consisting of Rp, tau, RMSE
    """
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    file_name = ntpath.basename(feature_file).split('.')[0]
    # file_name2 = ntpath.basename(feature_file2).split('.')[0]
    
    feature_df = pd.read_csv(feature_file)

    # df_2 = pd.read_csv(feature_file2)
        
    Rp, tau, RMSE, df_predict = do_Kfold(feature_df, ml_random_state=ml_random_state, kf_random_state=kf_random_state) 
    
    output_file_name = f'{output_folder}/{file_name}.{ML_type}.{kf_random_state}_{ml_random_state}'

    df_predict.to_csv(f'{output_file_name}.prediction',
             index=False, float_format='%.5f')

    df_results = pd.DataFrame(data={'Rp': [Rp], 'tau': [tau], 'RMSE': [RMSE]})

    df_results.to_csv(f'{output_file_name}.results',
             index=False, float_format='%.5f')
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Do KFold experiments and save their results to a csv file",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', '--ml-type', help='machine learning type, either gbt, or rf, or xgbt',
                        type=str, default='gbt')
    parser.add_argument('-o', '--output-folder', help='folder to save kfold results')                       

    parser.add_argument('-f', '--feature-file', help='feature file')

    parser.add_argument('-k', '--kf_random_state', type=int, 
                        default=1, help='kfold random state')
    
    

    
    args = parser.parse_args()

    print(args)
    

    feature_file = args.feature_file
    # feature_file2 = args.feature_file2

    ML_type = args.ml_type
    output_folder = args.output_folder
    kf_random_state = args.kf_random_state

    ml_random_state = randint(1,1000000) 

    main(feature_file, ML_type, ml_random_state, kf_random_state, output_folder)
        
#        